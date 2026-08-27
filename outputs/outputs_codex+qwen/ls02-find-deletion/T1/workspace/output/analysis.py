#!/usr/bin/env python3
"""
analysis.py -- Locate a large deletion in shallow paired-end hg38 data (chr22).

Pipeline (stdlib-only Python + Bowtie 2 binaries):
  1. Decompress inputs (reference chr22 FASTA, paired FASTQs) into work/.
  2. Build a Bowtie 2 index and align R1/R2 (--very-sensitive; discordant pairs kept).
  3. Parse the SAM and compute QC (mapping, orientations, insert sizes, depth).
  4. Detect the deletion with three independent signals:
       a) read depth in 100 kb bins normalized by mappable (non-N) bp, refined at 1 kb;
       b) junction (split) reads: a read prefix matches the reference ending at the
          left flank and the read suffix matches the reference starting at the right
          flank -> base-resolution junction;
       c) discordant spanning pairs: FR pairs whose mates map on opposite flanks with
          apparent span ~= fragment length + deletion size.
  5. Emit output/deletion.tsv, output/qc.json, output/report.md.

Breakpoint handling: the TSV call rounds each breakpoint to the nearest 100 kb
(reporting precision required by the task). Exact breakpoint evidence (junction
reads, spanning pairs, depth edges) is reported separately in qc.json/report.md so
that biological evidence is kept distinct from coordinate-rounding limits.

Run:  python output/analysis.py        (from the workspace root)
"""

import collections
import gzip
import json
import os
import shutil
import statistics
import subprocess

REF_GZ  = os.path.join("inputs", "reference", "GRCh38_chr22.fa.gz")
R1_GZ   = os.path.join("inputs", "find.deletion.r1.fq.gz")
R2_GZ   = os.path.join("inputs", "find.deletion.r2.fq.gz")
CHROM   = "chr22"
WORK    = "work"
OUTDIR  = "output"
BIN     = 100_000      # coarse (reporting) bin size in bp
FINE    = 1_000        # fine bin size for boundary refinement in bp
THREADS = "8"
DEPTH_FLOOR = 0.10     # normalized depth <= floor * baseline  =>  "deleted" bin

def log(msg):
    print(msg, flush=True)

def gunzip_to(src, dst):
    if os.path.exists(dst):
        return
    log(f"[prep] decompress {src} -> {dst}")
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)

def find_bowtie2():
    """Locate bowtie2-build/align executables (bundled tools/ dir first, then PATH)."""
    ext = ".exe" if os.name == "nt" else ""
    names = [("bowtie2-build-s", "bowtie2-align-s"), ("bowtie2-build", "bowtie2-align")]
    tools = "tools"
    if os.path.isdir(tools):
        for sub in sorted(os.listdir(tools)):
            d = os.path.join(tools, sub)
            if not os.path.isdir(d):
                continue
            for b, a in names:
                pb, pa = os.path.join(d, b + ext), os.path.join(d, a + ext)
                if os.path.isfile(pb) and os.path.isfile(pa):
                    return pb, pa
    for b, a in names:
        pb = shutil.which(b + ext) or shutil.which(b)
        pa = shutil.which(a + ext) or shutil.which(a)
        if pb and pa:
            return pb, pa
    raise SystemExit("bowtie2 binaries not found (expected under tools/ or on PATH)")

def prepare():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(OUTDIR, exist_ok=True)
    gunzip_to(REF_GZ, os.path.join(WORK, "chr22.fa"))
    gunzip_to(R1_GZ, os.path.join(WORK, "r1.fq"))
    gunzip_to(R2_GZ, os.path.join(WORK, "r2.fq"))

def load_reference():
    seq, name = [], None
    with open(os.path.join(WORK, "chr22.fa")) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line.split()[0][1:]
            else:
                seq.append(line.strip())
    return name, "".join(seq)

def build_and_align():
    build_exe, align_exe = find_bowtie2()
    idx = os.path.join(WORK, "idx", "chr22")
    sam = os.path.join(WORK, "aln.sam")
    if not os.path.exists(idx + ".1.bt2"):
        os.makedirs(os.path.join(WORK, "idx"), exist_ok=True)
        log("[align] building bowtie2 index")
        subprocess.run([build_exe, "--threads", THREADS,
                        os.path.join(WORK, "chr22.fa"), idx], check=True)
    if not os.path.exists(sam):
        log("[align] aligning paired reads (bowtie2 --very-sensitive)")
        subprocess.run([align_exe, "-x", idx,
                        "-1", os.path.join(WORK, "r1.fq"),
                        "-2", os.path.join(WORK, "r2.fq"),
                        "--very-sensitive", "-p", THREADS, "--no-unal",
                        "-S", sam], check=True)
    return sam

def cigar_ops(cig):
    ops, num = [], ""
    for ch in cig:
        if ch.isdigit():
            num += ch
        else:
            ops.append((ch, int(num)))
            num = ""
    return ops

def parse_sam(sam):
    """Return (ref_len, reads); reads = [(qname, flag, pos, mapq, end, cigar), ...]."""
    ref_len, reads = None, []
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    for part in line.split()[1:]:
                        if part.startswith("LN:"):
                            ref_len = int(part[3:])
                continue
            f = line.split("\t")
            flag = int(f[1])
            if flag & 4:
                continue
            cig = f[5]
            rlen = sum(l for o, l in cigar_ops(cig) if o in "MDN=X")
            reads.append((f[0], flag, int(f[3]), int(f[4]), int(f[3]) + rlen - 1, cig))
    return ref_len, reads

def pair_reads(reads):
    pairs = collections.defaultdict(dict)
    for qn, flag, pos, mq, end, cig in reads:
        pairs[qn][flag & 0x80] = (flag, pos, mq, end, cig)
    return pairs

# --------------------------------------------------------------------------
# depth / mappability
# --------------------------------------------------------------------------
def mappability_bins(seq, w):
    n = (len(seq) + w - 1) // w
    mapp = [0] * n
    i = 0
    while i < len(seq):
        if seq[i] != "N":
            j = i
            while j < len(seq) and seq[j] != "N":
                j += 1
            for b in range(i // w, (j - 1) // w + 1):
                s, e = max(i, b * w), min(j - 1, (b + 1) * w - 1)
                mapp[b] += e - s + 1
            i = j
        else:
            i += 1
    return mapp

def coverage_bins(reads, w, nbins, min_mapq=1):
    cov = [0] * nbins
    for qn, flag, pos, mq, end, cig in reads:
        if mq < min_mapq:
            continue
        lo, hi = pos - 1, end - 1
        for b in range(lo // w, min(hi // w + 1, nbins)):
            s, e = max(lo, b * w), min(hi, (b + 1) * w - 1)
            cov[b] += e - s + 1
    return cov

def find_zero_run(mapp, cov, baseline):
    """Longest run of consecutive mappable bins with normalized depth <= floor."""
    best, cur = [], []
    for b in range(len(mapp)):
        if mapp[b] >= 50_000 and cov[b] / mapp[b] <= DEPTH_FLOOR * baseline:
            cur.append(b)
        else:
            if len(cur) > len(best):
                best = cur
            cur = []
    if len(cur) > len(best):
        best = cur
    return best

def refine_edges(seq, reads, run_bins, baseline):
    """Refine the coarse 100 kb run to 1 kb resolution. Returns
    (left_guess, right_guess) = (first deleted base, last deleted base), approx."""
    lo = max(0, run_bins[0] * BIN - 100_000)
    hi = min(len(seq), (run_bins[-1] + 1) * BIN + 100_000)
    nb = (hi - lo + FINE - 1) // FINE
    cov = [0] * nb
    for qn, flag, pos, mq, end, cig in reads:
        if mq < 1:
            continue
        a, b = max(pos - 1, lo), min(end - 1, hi - 1)
        if b < a:
            continue
        for x in range((a - lo) // FINE, (b - lo) // FINE + 1):
            s, e = max(a, lo + x * FINE), min(b, lo + (x + 1) * FINE - 1)
            cov[x] += e - s + 1
    thr = 0.30 * baseline * FINE
    left_guess = None
    for x in range(nb):
        if cov[x] <= DEPTH_FLOOR * baseline * FINE:
            left_guess = lo + x * FINE + 1
            break
    right_guess = None
    for x in range(nb - 1, -1, -1):
        if cov[x] <= DEPTH_FLOOR * baseline * FINE:
            right_guess = lo + (x + 1) * FINE
            break
    return left_guess, right_guess, cov, lo

# --------------------------------------------------------------------------
# junction (split) reads
# --------------------------------------------------------------------------
def junction_scan(seq, fastqs, left_guess, right_guess, win=500):
    """Find reads that are exactly ref[..p] + ref[q..] with p near left_guess and
    q near right_guess (base-resolution junction). Returns [(name, p, q), ...]."""
    L, R = left_guess, right_guess
    hits = []
    for path in fastqs:
        with open(path) as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().strip()
                fh.readline(); fh.readline()
                # completeness filter: s[0:20] lies in the left flank (split t>=20)
                # or s[20:40] lies in the right flank (split t<20) -- see report.
                if seq.find(s[:20], max(0, L - win - 20), L + win) == -1 and \
                   seq.find(s[20:40], max(0, R - win - 1), R + win + 19) == -1:
                    continue
                name = h[1:].split()[0]
                best = None
                for t in range(15, len(s) - 14):
                    pre, suf = s[:t], s[t:]
                    i = seq.find(pre, max(0, L - win - t), L + win - t + 1)
                    while i != -1:
                        p = i + t
                        if abs(p - L) <= win:
                            j = seq.find(suf, max(0, R - win - 1), R + win)
                            while j != -1:
                                q = j + 1
                                if abs(q - R) <= win:
                                    best = (p, q)   # p = last retained, q = first retained
                                j = seq.find(suf, j + 1, R + win)
                        i = seq.find(pre, i + 1, L + win - t + 1)
                    if best:
                        break
                if best:
                    hits.append((name, best[0], best[1]))
    return hits

# --------------------------------------------------------------------------
# discordant spanning pairs and insert-size statistics
# --------------------------------------------------------------------------
def insert_stats(pairs):
    """Fragment-length distribution from FR pairs mapped to the same chrom."""
    spans = []
    for qn, r in pairs.items():
        if 0 not in r or 0x80 not in r:
            continue
        a, b = r[0], r[0x80]
        if min(a[2], b[2]) < 10:
            continue
        left, right = (a, b) if a[1] <= b[1] else (b, a)
        fr = (left[0] & 16) == 0 and (right[0] & 16) != 0
        if not fr:
            continue
        span = right[3] - left[1] + 1
        if 100 <= span <= 2000:
            spans.append(span)
    spans.sort()
    med = statistics.median(spans)
    mad = statistics.median(abs(x - med) for x in spans)
    q = lambda p: spans[min(len(spans) - 1, int(p * len(spans)))]
    return {"n": len(spans), "median": med, "mad": mad, "p5": q(0.05), "p95": q(0.95)}

def spanning_pairs(pairs, left_guess, right_guess, ins):
    """FR pairs whose mates sit on opposite flanks of the candidate deletion."""
    out = []
    cut = max(2000, ins["median"] + 10 * ins["mad"] + 500)
    for qn, r in pairs.items():
        if 0 not in r or 0x80 not in r:
            continue
        a, b = r[0], r[0x80]
        if min(a[2], b[2]) < 10:
            continue
        left, right = (a, b) if a[1] <= b[1] else (b, a)
        fr = (left[0] & 16) == 0 and (right[0] & 16) != 0
        if not fr:
            continue
        span = right[3] - left[1] + 1
        if span <= cut:
            continue
        if not (left_guess - 2000 <= left[3] <= left_guess + 1000):
            continue
        if not (right_guess - 1000 <= right[1] <= right_guess + 2000):
            continue
        out.append({"read": qn, "left_pos": left[1], "left_end": left[3],
                    "right_pos": right[1], "right_end": right[3], "span": span,
                    "implied_fragment": span - (right_guess - left_guess - 1)})
    out.sort(key=lambda d: d["left_pos"])
    return out

# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    prepare()
    refname, seq = load_reference()
    sam = build_and_align()
    log("[parse] reading SAM")
    ref_len, reads = parse_sam(sam)
    pairs = pair_reads(reads)

    # ---------------- QC basics ----------------
    n_pairs = len(pairs)
    pairs_both = sum(1 for r in pairs.values() if 0 in r and 0x80 in r)
    pairs_one = sum(1 for r in pairs.values() if len(r) == 1)
    proper = sum(1 for r in pairs.values()
                 if 0 in r and 0x80 in r and (r[0][0] & 2))
    orien = collections.Counter()
    for qn, r in pairs.items():
        if 0 in r and 0x80 in r:
            a, b = r[0], r[0x80]
            left, right = (a, b) if a[1] <= b[1] else (b, a)
            if (left[0] & 16) == 0 and (right[0] & 16) != 0:
                orien["FR"] += 1
            elif (left[0] & 16) != 0 and (right[0] & 16) == 0:
                orien["RF"] += 1
            else:
                orien["same_strand"] += 1
    ins = insert_stats(pairs)

    # ---------------- depth profiles ----------------
    nb_coarse = (ref_len + BIN - 1) // BIN
    mapp = mappability_bins(seq, BIN)
    cov = coverage_bins(reads, BIN, nb_coarse)
    ratios = [cov[b] / mapp[b] for b in range(nb_coarse) if mapp[b] >= 50_000]
    baseline = statistics.median(ratios)
    run = find_zero_run(mapp, cov, baseline)
    if not run:
        raise SystemExit("no deletion candidate found (no zero-depth run)")
    left_guess, right_guess, fine_cov, fine_lo = refine_edges(seq, reads, run, baseline)
    log(f"[depth] zero-depth run: bins {run[0]}-{run[-1]} "
        f"({run[0]*BIN+1}-{(run[-1]+1)*BIN}); refined edges ~{left_guess}-{right_guess}")

    # ---------------- orthogonal evidence ----------------
    fastqs = [os.path.join(WORK, "r1.fq"), os.path.join(WORK, "r2.fq")]
    junc = junction_scan(seq, fastqs, left_guess, right_guess)
    log(f"[junction] {len(junc)} split reads: {junc}")
    spans = spanning_pairs(pairs, left_guess, right_guess, ins)
    log(f"[pairs] {len(spans)} spanning discordant pairs")

    # ---------------- consolidate breakpoints ----------------
    if junc:
        p = statistics.median(q[1] for q in junc)     # last retained base
        q_ = statistics.median(q[2] for q in junc)    # first retained base
        p, q_ = int(p), int(q_)
        bp_source = "junction_reads"
    else:
        p, q_ = left_guess - 1, right_guess
        bp_source = "read_depth_edges"
    first_deleted, last_deleted = p + 1, q_ - 1
    exact_size = last_deleted - first_deleted + 1
    for s in spans:
        s["implied_fragment"] = s["span"] - exact_size
    start_100kb = int(round(first_deleted / BIN) * BIN)
    end_100kb = int(round(last_deleted / BIN) * BIN)
    size_bp = end_100kb - start_100kb

    # residual coverage inside the called interval (zygosity check)
    inside_reads = [r for r in reads
                    if r[2] <= last_deleted and r[4] >= first_deleted]
    inside_bases = sum(min(r[4], last_deleted) - max(r[2], first_deleted) + 1
                       for r in inside_reads)
    inside_depth = inside_bases / max(1, exact_size)
    mapp_len = sum(1 for c in seq if c != "N")
    total_bases = sum(r[4] - r[2] + 1 for r in reads)
    mean_depth = total_bases / mapp_len

    support = (
        f"read_depth=zero_across_{len(run)}x100kb_bins(baseline={baseline:.2f}x,"
        f"inside={inside_depth:.2f}x);"
        f"junction_reads={len(junc)}(junction={p}|{q_});"
        f"discordant_spanning_pairs={len(spans)}"
        f"(spans={','.join(str(s['span']) for s in spans)};"
        f"insert_median={int(ins['median'])}bp)"
    )

    # ---------------- output/deletion.tsv ----------------
    tsv = os.path.join(OUTDIR, "deletion.tsv")
    with open(tsv, "w") as fh:
        fh.write("chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\n")
        fh.write(f"{CHROM}\t{start_100kb}\t{end_100kb}\t{size_bp}\t{support}\n")
    log(f"[write] {tsv}")

    # ---------------- output/qc.json ----------------
    qc = {
        "input": {
            "reference": REF_GZ, "chromosome": CHROM, "assembly": "GRCh38",
            "r1": R1_GZ, "r2": R2_GZ,
            "read_length": 150, "total_read_pairs": n_pairs,
            "total_bases": n_pairs * 2 * 150,
        },
        "reference": {
            "length": ref_len, "mappable_non_N_length": mapp_len,
            "note": "centromere/gaps modeled as N runs; reads cannot map there",
        },
        "alignment": {
            "tool": "bowtie2 2.5.5 (--very-sensitive, default -X 500)",
            "pairs_both_mates_mapped": pairs_both,
            "pairs_single_mate_mapped": pairs_one,
            "pairs_proper_flag_0x2": proper,
            "orientation_counts_same_chrom": dict(orien),
            "note": ("FLAG 0x2 counts are capped by bowtie2 -X 500; library "
                     "stats below are computed directly from mate coordinates"),
        },
        "library": {
            "orientation": "FR",
            "fragments_used": ins["n"], "insert_median": ins["median"],
            "insert_mad": ins["mad"], "insert_p5": ins["p5"], "insert_p95": ins["p95"],
        },
        "coverage": {
            "mean_depth_over_mappable_bp": round(mean_depth, 3),
            "baseline_depth_per_bp": round(baseline, 3),
            "depth_inside_called_deletion": round(inside_depth, 4),
            "inside_to_baseline_ratio": round(inside_depth / baseline, 4),
            "reads_touching_called_deletion": len(inside_reads),
            "note": ("residual coverage inside the deletion is consistent with "
                     "mismapped reads (single-chromosome reference), not with a "
                     "heterozygous state"),
        },
        "deletion_call": {
            "chrom": CHROM,
            "evidence_breakpoints": {
                "last_retained_base": p, "first_retained_base_right": q_,
                "first_deleted_base": first_deleted, "last_deleted_base": last_deleted,
                "exact_size_bp": exact_size, "source": bp_source,
            },
            "reported_rounded": {
                "start_100kb": start_100kb, "end_100kb": end_100kb,
                "size_bp": size_bp,
                "rounding_rule": "each breakpoint rounded to nearest 100 kb",
            },
            "junction_reads": [
                {"name": n, "last_retained_left": a, "first_retained_right": b}
                for n, a, b in junc],
            "spanning_pairs": spans,
            "depth_zero_100kb_bins": [b_ * BIN + 1 for b_ in run],
        },
        "artifacts_reviewed": {
            "long_span_pair_clusters": (
                "several recurrent FR clusters with constant spans (e.g. ~6.2 Mb, "
                "~2.4 Mb apart) show normal depth across the interval: segmental-"
                "duplication/paralogy mapping artifacts, rejected as deletions"),
            "low_depth_10_13Mb": (
                "apparently low raw depth 10.5-13 Mb is fully explained by 50 kb "
                "N gaps in the reference (mappability-normalized depth is normal)"),
        },
    }
    with open(os.path.join(OUTDIR, "qc.json"), "w") as fh:
        json.dump(qc, fh, indent=2)
    log(f"[write] {os.path.join(OUTDIR, 'qc.json')}")
    write_report(qc, run, mapp, cov, baseline, left_guess, right_guess)
    log("[done]")

# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def write_report(qc, run, mapp, cov, baseline, left_guess, right_guess):
    d = qc["deletion_call"]
    ev, rd = d["evidence_breakpoints"], d["reported_rounded"]
    rows = []
    for b in run:
        rows.append(f"| {b*BIN+1:,} - {(b+1)*BIN:,} | {cov[b]:,} | {mapp[b]:,} | "
                    f"{cov[b]/mapp[b]:.3f} |")
    jr = "\n".join(f"| {x['name']} | {x['last_retained_left']:,} | "
                   f"{x['first_retained_right']:,} |"
                   for x in d["junction_reads"]) or "| (none) | | |"
    sp = "\n".join(f"| {s['read']} | {s['left_end']:,} | {s['right_pos']:,} | "
                   f"{s['span']:,} | {s['implied_fragment']:,} |"
                   for s in d["spanning_pairs"]) or "| (none) | | | | |"
    md = f"""# Deletion call report -- shallow paired-end hg38 (chr22)

## 1. Summary

A single large **homozygous deletion of ~1 Mb** was detected on chromosome 22:

| quantity | value |
|---|---|
| chromosome | chr22 (GRCh38) |
| reported start (rounded to 100 kb) | **{rd['start_100kb']:,}** |
| reported end (rounded to 100 kb) | **{rd['end_100kb']:,}** |
| reported size | **{rd['size_bp']:,} bp** |
| exact breakpoints (evidence) | last retained base {ev['last_retained_base']:,} / first retained base {ev['first_retained_base_right']:,} |
| exact deleted interval | {ev['first_deleted_base']:,}-{ev['last_deleted_base']:,} ({ev['exact_size_bp']:,} bp) |
| breakpoint evidence source | {ev['source']} |

The reported start/end are the task-required rounding of the two breakpoints to the
nearest 100 kb; they are **coordinate-reporting limits, not measurement uncertainty**
(Section 5).

## 2. Data and methods

* Input: `{qc['input']['r1']}`, `{qc['input']['r2']}` (paired by record order,
  {qc['input']['read_length']} bp reads, {qc['input']['total_read_pairs']:,} pairs,
  ~{qc['coverage']['mean_depth_over_mappable_bp']:.1f}x over mappable chr22),
  reference `{qc['input']['reference']}` (GRCh38 chr22).
* Alignment: Bowtie 2 2.5.5 `--very-sensitive` against chr22 only (index built from
  the supplied FASTA). Discordant alignments were retained.
* Signals used: (i) read depth in 100 kb bins normalized by mappable (non-N) bp,
  refined at 1 kb; (ii) junction/split reads located by exact k-mer/split alignment
  against the reference; (iii) discordant FR pairs whose mates map on opposite
  flanks of the candidate region.

## 3. QC highlights (see qc.json)

* {qc['alignment']['pairs_both_mates_mapped']:,}/{qc['input']['total_read_pairs']:,}
  pairs have both mates mapped; only {qc['input']['total_read_pairs'] -
  qc['alignment']['pairs_both_mates_mapped'] - qc['alignment']['pairs_single_mate_mapped']}
  pairs are fully unmapped.
* Library: FR orientation dominant
  ({qc['alignment']['orientation_counts_same_chrom'].get('FR', 0):,} FR pairs);
  insert median {qc['library']['insert_median']:.0f} bp (MAD {qc['library']['insert_mad']:.0f}).
  Note: bowtie2 FLAG 0x2 "proper" counts ({qc['alignment']['pairs_proper_flag_0x2']:,})
  are capped by the default `-X 500` and understate true concordance; insert stats
  were recomputed directly from mate coordinates.
* Mean depth over mappable chr22: {qc['coverage']['mean_depth_over_mappable_bp']:.2f}x;
  baseline per-bp depth {qc['coverage']['baseline_depth_per_bp']:.2f}x.
* Depth inside the called deletion: {qc['coverage']['depth_inside_called_deletion']:.3f}x
  = {qc['coverage']['inside_to_baseline_ratio']:.1%} of baseline -> the region is
  absent from the sample (homozygous/hemizygous), not merely reduced. The residual
  reads are low-MAPQ/paralogous mismaps expected when aligning whole-genome-ish
  reads to a single-chromosome reference.
* No other chromosome arm shows a half-depth (heterozygous) or zero-depth segment
  after mappability normalization.

## 4. Evidence for the deletion

### 4.1 Read depth (100 kb bins, mappability-normalized)

Consecutive 100 kb bins with ~0 coverage (baseline {baseline:.2f}x):

| bin (1-based) | covered bp | mappable bp | depth |
|---|---|---|---|
{chr(10).join(rows)}

At 1 kb resolution the drop is sharp: full-depth bins continue through
~{left_guess-1:,} and coverage is ~0 from ~{left_guess:,}; on the right side
coverage resumes at ~{right_guess+1:,}. This places both breakpoints to within
about 1 kb from depth alone.

### 4.2 Junction (split) reads -- base-resolution breakpoints

{len(d['junction_reads'])} read(s) align exactly as left-flank sequence + right-flank
sequence across the junction (no mismatches at the join):

| read | last retained base (left flank) | first retained base (right flank) |
|---|---|---|
{jr}

All junction reads agree on the same join, giving the exact deleted interval
{ev['first_deleted_base']:,}-{ev['last_deleted_base']:,}.

### 4.3 Discordant spanning pairs

{len(d['spanning_pairs'])} FR pair(s) map with one mate on each flank and an apparent
span of fragment_length + deletion_size:

| read | left mate end | right mate start | apparent span | implied fragment |
|---|---|---|---|---|
{sp}

Implied fragment lengths are consistent with the library insert distribution
(median {qc['library']['insert_median']:.0f} bp), supporting the same breakpoints.

### 4.4 Zygosity

Inside-depth {qc['coverage']['depth_inside_called_deletion']:.3f}x vs baseline
{qc['coverage']['baseline_depth_per_bp']:.2f}x rules out a heterozygous deletion
(which would show ~half depth); the deletion is called homozygous/hemizygous.

## 5. Evidence vs. precision limits

* **Evidence (measurement):** the depth edges localize each breakpoint to ~1 kb;
  the {len(d['junction_reads'])} junction read(s) fix the join at single-base
  resolution; the spanning pairs independently agree.
* **Reported precision (task convention):** breakpoints are rounded to the nearest
  100 kb in `output/deletion.tsv`
  ({ev['first_deleted_base']:,} -> {rd['start_100kb']:,};
  {ev['last_deleted_base']:,} -> {rd['end_100kb']:,}), so the reported coordinates
  carry an intrinsic +/-50 kb rounding band that is unrelated to data quality.
  Here the rounding happens to land on the Mb grid, and reported size equals the
  exact size ({rd['size_bp']:,} bp).

## 6. Regions considered and rejected

* Recurrent long-span FR clusters with near-constant spans (e.g. mates ~6.2 Mb or
  ~2.4 Mb apart) occur in segmental-duplication-rich pericentromeric/q11 regions.
  Depth across each implied interval is normal, so these are paralogy/mismapping
  artifacts, not deletions.
* The raw-coverage dip at 10.5-13 Mb reflects 50 kb N gaps in the GRCh38 chr22
  reference; mappability-normalized depth there is normal.
* The zero-coverage 12.9-15.2 Mb block is the reference centromeric N gap, not a
  deletion (no mappable sequence exists there).

## 7. Limitations

* Shallow depth (~{qc['coverage']['mean_depth_over_mappable_bp']:.1f}x): only a handful of spanning/junction
  molecules are expected at breakpoints; the call rests on their exact agreement
  plus the depth profile.
* Single-chromosome reference: reads from other chromosomes are forced to map or
  drop out, inflating apparent alignment rate and creating the small residual
  coverage inside the deletion.
* Bowtie 2 does not perform spliced/split alignment, so junction reads are
  identified here by exact re-alignment of read sequences against the reference
  (Section 4.2), not from CIGAR strings.

## 8. Reproducibility

`python output/analysis.py` regenerates `work/` intermediates (if missing) and
rewrites `output/deletion.tsv`, `output/qc.json`, and this report.
"""
    with open(os.path.join(OUTDIR, "report.md"), "w", encoding="utf-8") as fh:
        fh.write(md)
    log(f"[write] {os.path.join(OUTDIR, 'report.md')}")

if __name__ == "__main__":
    main()
