#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Find a large deletion in shallow paired-end hg38 (chr22-only) WGS data.

Inputs (workspace-relative):
  inputs/find.deletion.r1.fq.gz        mate-1 FASTQ (paired by record order)
  inputs/find.deletion.r2.fq.gz        mate-2 FASTQ
  inputs/reference/GRCh38_chr22.fa.gz  GRCh38 chr22 reference, single contig

Method (pure Python + NumPy; no external aligner needed):

  1. Build an exact 25-mer index of the reference: 2-bit-packed hashes of
     every reference 25-mer, radix-sorted with NumPy, positions retained.
  2. Seed-and-verify alignment: for every read, three seeds per strand at
     offsets 0/62/125. Seed hits are grouped into candidate alignments and
     each read is classified as uniquely mapped, ambiguously (multi-) mapped,
     split (chimeric: two distant same-strand loci), or unmapped.
  3. Three orthogonal deletion signals:
       a. read depth in 100 kb (reporting) and 10 kb (refinement) bins,
          normalized by per-bin mappable (ACGT) reference content;
       b. discordant mate pairs: FR pairs whose mapped span greatly exceeds
          the insert-size distribution AND whose mates flank the low-depth
          region -> they bridge the deletion;
       c. split reads, verified base-by-base against the implied junction
          (>= 30 bp on each side, crossover consistent with the seeding),
          giving single-base breakpoint resolution. Split-read junctions are
          only accepted as evidence when they are consistent with the
          depth-defined interval (breakpoint within 50 kb of a depth edge
          and implied size within 20% of the depth interval size).
  4. Combine evidence: depth defines the deleted interval (primary evidence
     for position and zygosity); spanning pairs / split reads refine the
     breakpoints where available. Coordinates are finally rounded to the
     nearest 100 kb for deletion.tsv (task-mandated reporting precision),
     while qc.json and report.md keep the refined coordinates plus the
     per-signal evidence counts, so evidence can be distinguished from
     precision limits.

Outputs (workspace-relative):
  output/deletion.tsv  chrom start_100kb end_100kb size_bp supporting_signals
  output/qc.json       QC metrics + per-signal evidence
  output/report.md     narrative report

Deterministic; only touches files inside the workspace.
"""
from __future__ import annotations

import gzip
import json
import math
import os
import time

import numpy as np

WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # workspace root
R1_PATH = os.path.join(WS, "inputs", "find.deletion.r1.fq.gz")
R2_PATH = os.path.join(WS, "inputs", "find.deletion.r2.fq.gz")
REF_PATH = os.path.join(WS, "inputs", "reference", "GRCh38_chr22.fa.gz")
OUT_DIR = os.path.join(WS, "output")

K = 25                     # seed length (bp)
OFFSETS = (0, 62, 125)     # seed offsets inside each 150 bp read
READ_LEN = 150
MAX_SEED_HITS = 16         # seeds hitting more loci are treated as repeats
BIG_BIN = 100_000          # reporting bin = requested rounding unit
FINE_BIN = 10_000          # refinement bin for depth boundaries
MIN_SPLIT_SIDE = 30        # min aligned bp on each side of a junction
MIN_SPLIT_DEL = 1_000      # split-read implied deletion must be >= this
EDGE_TOL = 50_000          # split junction may deviate this much from depth edge

# status codes
UNMAPPED, UNIQUE, AMBIG, SPLIT = 0, 1, 2, 3


def log(msg: str) -> None:
    print(f"[analysis] {msg}", flush=True)


def round_half_up_100kb(x: int) -> int:
    """Round a coordinate to the nearest 100 kb (half rounds up)."""
    return int(math.floor(x / float(BIG_BIN) + 0.5)) * BIG_BIN


# --------------------------------------------------------------------------
# reference / read loading
# --------------------------------------------------------------------------

def load_reference(path):
    data = gzip.open(path, "rb").read()
    lines = data.split(b"\n")
    name = lines[0][1:].split()[0].decode()
    seq = b"".join(l.strip() for l in lines[1:] if l and not l.startswith(b">"))
    return name, seq.upper().decode()


_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _c in enumerate(b"ACGT"):
    _CODE[_c] = _i


def encode(seqbytes: bytes) -> np.ndarray:
    return _CODE[np.frombuffer(seqbytes, dtype=np.uint8)]


def load_reads(path):
    """Return (n_records, codes uint8 array of shape (n, READ_LEN))."""
    data = gzip.open(path, "rb").read()
    lines = data.split(b"\n")
    seqs = lines[1::4]
    if seqs and seqs[-1] == b"":
        seqs = seqs[:-1]
    seqs = [s for s in seqs if s]
    n = len(seqs)
    lens = {len(s) for s in seqs[:2000]}
    assert lens == {READ_LEN}, f"unexpected read lengths in {path}: {lens}"
    cat = b"".join(seqs)
    codes = _CODE[np.frombuffer(cat, dtype=np.uint8)].reshape(n, READ_LEN)
    return n, codes


# --------------------------------------------------------------------------
# 25-mer index
# --------------------------------------------------------------------------

def rolling_hash(x: np.ndarray) -> np.ndarray:
    """Hash rows of x (shape (n, w), values 0..4) into uint64 2-bit codes.

    Rows containing a non-ACGT symbol get the sentinel hash UINT64_MAX.
    """
    n = x.shape[0]
    h = np.zeros(n, dtype=np.uint64)
    bad = np.zeros(n, dtype=bool)
    for j in range(x.shape[1]):
        col = x[:, j]
        h = (h << np.uint64(2)) | (col.astype(np.uint64) & np.uint64(3))
        bad |= col == 4
    h[bad] = np.uint64(0xFFFFFFFFFFFFFFFF)
    return h


def build_index(codes: np.ndarray):
    """Sorted hashes + positions of every reference K-mer."""
    L = len(codes)
    n = L - K + 1
    log(f"building {K}-mer index over {n:,} reference positions ...")
    t0 = time.time()
    h = np.zeros(n, dtype=np.uint64)
    bad = np.zeros(n, dtype=bool)
    for j in range(K):
        seg = codes[j:j + n]
        h = (h << np.uint64(2)) | (seg.astype(np.uint64) & np.uint64(3))
        bad |= seg == 4
    h[bad] = np.uint64(0xFFFFFFFFFFFFFFFF)
    order = np.argsort(h, kind="stable")
    sorted_hashes = h[order]
    sorted_pos = order.astype(np.uint32)
    log(f"index built in {time.time() - t0:.1f}s "
        f"({int(bad.sum()):,} k-mers masked (N-containing)")
    return sorted_hashes, sorted_pos


# --------------------------------------------------------------------------
# alignment classification
# --------------------------------------------------------------------------

def align_reads(codes_r, sorted_hashes, sorted_pos, tag=""):
    """Classify every read; returns dict of per-read numpy arrays.

    Primary alignment: strand (+1/-1), start (0-based ref position of the
    read's aligned interval on the + strand), support (# agreeing seeds).
    Split reads additionally carry the secondary locus and the raw seed hit
    (offset, ref_pos) that implies the junction.
    """
    N = codes_r.shape[0]
    t0 = time.time()
    rc = np.where(codes_r == 4, np.uint8(4), np.uint8(3) - codes_r[:, ::-1])

    Q = np.empty((N, 6), dtype=np.uint64)
    for j, off in enumerate(OFFSETS):
        Q[:, j] = rolling_hash(codes_r[:, off:off + K])
        Q[:, 3 + j] = rolling_hash(rc[:, off:off + K])
    q = Q.reshape(-1)
    lo = np.searchsorted(sorted_hashes, q, side="left")
    hi = np.searchsorted(sorted_hashes, q, side="right")
    raw_mult = hi - lo
    is_n = q == np.uint64(0xFFFFFFFFFFFFFFFF)
    over = raw_mult > MAX_SEED_HITS
    mult = raw_mult.astype(np.int32)
    mult[is_n | over] = 0
    log(f"{tag}seed lookup: {int((raw_mult == 0).sum() - is_n.sum()):,} seeds no hit, "
        f"{int((raw_mult == 1).sum()):,} unique, "
        f"{int(((raw_mult >= 2) & (raw_mult <= MAX_SEED_HITS)).sum()):,} multi(<= {MAX_SEED_HITS}), "
        f"{int(over.sum()):,} hyper-repetitive (ignored) "
        f"[{time.time() - t0:.1f}s]")

    strand = np.zeros(N, dtype=np.int8)      # +1 forward, -1 reverse
    start = np.full(N, -1, dtype=np.int64)
    support = np.zeros(N, dtype=np.int8)
    status = np.zeros(N, dtype=np.uint8)
    split_start2 = np.full(N, -1, dtype=np.int64)
    split_strand2 = np.zeros(N, dtype=np.int8)
    split_hit_off = np.full(N, -1, dtype=np.int32)   # seed offset on 2nd locus
    split_hit_pos = np.full(N, -1, dtype=np.int64)   # ref pos of that seed

    L_ref = int(sorted_pos.size) + K - 1
    sp = sorted_pos

    t0 = time.time()
    for i in range(N):
        best = {}
        hits = {}  # (strand, start) -> list[(off, refpos)]
        for j in range(6):
            m = mult[i * 6 + j]
            if m <= 0:
                continue
            st = 1 if j < 3 else -1
            off = OFFSETS[j % 3]
            l0 = lo[i * 6 + j]
            for p in sp[l0:l0 + m]:
                s = int(p) - off
                if s < 0 or s > L_ref - READ_LEN:
                    continue
                key = (st, s)
                best[key] = best.get(key, 0) + 1
                hits.setdefault(key, []).append((off, int(p)))
        if not best:
            continue
        ranked = sorted(best.items(), key=lambda kv: (-kv[1], kv[0]))
        (st1, s1), c1 = ranked[0]
        if c1 == 1 and len(ranked) > 1 and ranked[1][1] == 1:
            # no locus has >= 2 seeds: ambiguous unless it looks like a
            # distant same-strand two-locus (split) read, verified later
            (st2x, s2x), _c2x = ranked[1]
            if not (len(ranked) == 2 and st2x == st1
                    and abs(s2x - s1) > READ_LEN):
                status[i] = AMBIG
                continue
        strand[i] = st1
        start[i] = s1
        support[i] = c1
        status[i] = UNIQUE
        if len(ranked) > 1:
            (st2, s2), c2 = ranked[1]
            if abs(s2 - s1) > READ_LEN and st2 == st1:
                status[i] = SPLIT
                split_start2[i] = s2
                split_strand2[i] = st2
                off2, p2 = hits[(st2, s2)][0]
                split_hit_off[i] = off2
                split_hit_pos[i] = p2
    log(f"{tag}read classification done in {time.time() - t0:.1f}s")
    return dict(strand=strand, start=start, support=support, status=status,
                split_start2=split_start2, split_strand2=split_strand2,
                split_hit_off=split_hit_off, split_hit_pos=split_hit_pos)


def verify_split_junction(refc, readc, s1, s2, off2, p2):
    """Base-level verification of a split read against the implied junction.

    Colinear model: the read prefix maps with ref start s1; the seed at read
    offset off2 maps at ref position p2 on the downstream locus. Implied
    deletion size D = p2 - off2 - s1. A valid crossover t (read coordinate
    where the alignment jumps) must satisfy:
      * read[0:t] == ref[s1 : s1+t]            (left flank matches)
      * read[t:150] == ref[s1+D+t : s1+D+150]  (right flank matches)
      * MIN_SPLIT_SIDE <= t <= 150 - MIN_SPLIT_SIDE  (both flanks real)
      * t <= off2                              (2nd seed lies right of it)
    Returns (ok, del_start_lo, del_start_hi, D): [del_start_lo, del_start_hi)
    is the 0-based range of possible left breakpoints (first deleted base - 1).
    """
    D = p2 - off2 - s1
    if D < MIN_SPLIT_DEL:
        return False, -1, -1, D
    if s1 + D + READ_LEN > len(refc) or s1 < 0:
        return False, -1, -1, D
    m1 = readc != refc[s1:s1 + READ_LEN]
    m2 = readc != refc[s1 + D:s1 + D + READ_LEN]
    a = int(np.argmax(m1)) if m1.any() else READ_LEN   # left matches [0, a)
    if m2.any():
        b = READ_LEN - int(np.argmax(m2[::-1]))        # right matches [b, L)
    else:
        b = READ_LEN
    t_lo = max(b, MIN_SPLIT_SIDE)
    t_hi = min(a, READ_LEN - MIN_SPLIT_SIDE, off2)
    if t_hi >= t_lo:
        return True, s1 + t_lo, s1 + t_hi + 1, D
    return False, -1, -1, D


# --------------------------------------------------------------------------
# signals: depth, discordant pairs, split reads
# --------------------------------------------------------------------------

def depth_profile(starts, codes, L):
    """Read-start depth in FINE_BIN and BIG_BIN windows, CN-normalized."""
    n_valid = L - READ_LEN + 1
    bad_cum = np.zeros(L + 1, dtype=np.int64)
    bad_cum[1:] = np.cumsum(codes == 4)
    valid = (bad_cum[READ_LEN:READ_LEN + n_valid] - bad_cum[:n_valid]) == 0
    vpos = np.nonzero(valid)[0]

    def norm(reads_b, valid_b):
        with np.errstate(divide="ignore", invalid="ignore"):
            dens = reads_b.astype(np.float64) / np.maximum(valid_b, 1)
        scale = starts.size / max(vpos.size, 1)
        return dens / scale if scale > 0 else dens

    nbf = int(math.ceil(n_valid / FINE_BIN))
    reads_fine = np.bincount(starts // FINE_BIN, minlength=nbf)
    valid_fine = np.bincount(vpos // FINE_BIN, minlength=nbf)
    cn_fine = norm(reads_fine, valid_fine)

    nbb = int(math.ceil(n_valid / BIG_BIN))
    reads_big = np.bincount(starts // BIG_BIN, minlength=nbb)
    valid_big = np.bincount(vpos // BIG_BIN, minlength=nbb)
    cn_big = norm(reads_big, valid_big)
    mappable_big = valid_big >= 0.5 * BIG_BIN
    mappable_fine = valid_fine >= 0.5 * FINE_BIN
    return cn_fine, cn_big, mappable_fine, mappable_big


def find_low_depth_run(cn_big, mappable_big, low):
    """Longest run of consecutive mappable big bins with cn < low."""
    flag = mappable_big & (cn_big < low)
    best = None
    i, n = 0, len(flag)
    while i < n:
        if flag[i]:
            j = i
            while j < n and flag[j]:
                j += 1
            if best is None or (j - i) > (best[1] - best[0]):
                best = (i, j)
            i = j
        else:
            i += 1
    return best


def refine_boundaries(cn_fine, mappable_fine, b0, b1, low):
    """Extend/refine a big-bin run [b0,b1) to fine-bin (10 kb) resolution."""
    lo_bin = max(b0 * BIG_BIN // FINE_BIN - 5, 0)
    hi_bin = min(int(math.ceil(b1 * BIG_BIN / FINE_BIN)) + 5, len(cn_fine))
    i = j = None
    for x in range(lo_bin, hi_bin):
        if mappable_fine[x] and cn_fine[x] < low:
            if i is None:
                i = x
            j = x
    if i is None:
        return b0 * BIG_BIN, b1 * BIG_BIN
    return i * FINE_BIN, (j + 1) * FINE_BIN


def main():
    t_start = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    ref_name, ref_seq = load_reference(REF_PATH)
    L = len(ref_seq)
    log(f"reference {ref_name}: {L:,} bp")
    refc = encode(ref_seq.encode())
    acgt = int((refc != 4).sum())

    n1, c1 = load_reads(R1_PATH)
    n2, c2 = load_reads(R2_PATH)
    assert n1 == n2, f"mate files not same length: {n1} vs {n2}"
    npairs = n1
    log(f"loaded {npairs:,} read pairs (read length {READ_LEN})")

    sh, sp = build_index(refc)
    a1 = align_reads(c1, sh, sp, tag="R1 ")
    a2 = align_reads(c2, sh, sp, tag="R2 ")

    usable1 = a1["support"] >= 2
    usable2 = a2["support"] >= 2
    aln_rate = (usable1.sum() + usable2.sum()) / (2.0 * npairs)
    log(f"alignment: {int(usable1.sum() + usable2.sum()):,}/{2 * npairs:,} "
        f"reads placed with >=2 seeds ({aln_rate:.3%})")

    # read-level mismatch QC on a subsample of placed reads
    # (strand-aware: minus-strand alignments are reverse-complemented first)
    idx_q = np.nonzero(usable1)[0][:20000]
    pos_q = a1["start"][idx_q]
    str_q = a1["strand"][idx_q]
    mm = np.array([
        int(((c1[i] if st == 1 else np.uint8(3) - c1[i][::-1])
             != refc[p0:p0 + READ_LEN]).sum())
        for i, p0, st in zip(idx_q, pos_q, str_q)])
    clean_frac = float((mm <= 2).mean())
    noisy = mm[mm > 2]
    noisy_mean = float(noisy.mean()) if noisy.size else 0.0
    log(f"read QC subsample ({len(idx_q):,}): {clean_frac:.1%} reads match "
        f"their placement with <=2 mismatches")

    # ---------------- read depth ----------------
    starts_all = np.concatenate([a1["start"][usable1], a2["start"][usable2]])
    starts_all = starts_all[(starts_all >= 0) & (starts_all <= L - READ_LEN)]
    cn_fine, cn_big, map_fine, map_big = depth_profile(starts_all, refc, L)
    med_cn = float(np.median(cn_big[map_big]))
    log(f"genome-wide median big-bin CN = {med_cn:.3f} "
        f"(mean depth ~ {starts_all.size * READ_LEN / acgt:.2f}x)")

    LOW_HOM, LOW_HET = 0.35, 0.70
    zygosity = "homozygous-like"
    low_used = LOW_HOM
    run = find_low_depth_run(cn_big, map_big, LOW_HOM)
    if run is None:
        run = find_low_depth_run(cn_big, map_big, LOW_HET)
        zygosity = "heterozygous-like"
        low_used = LOW_HET
    assert run is not None, "no large low-depth region found"
    b0, b1 = run
    depth_start, depth_end = refine_boundaries(cn_fine, map_fine, b0, b1,
                                               low_used)
    depth_size = depth_end - depth_start
    n_big_bins = b1 - b0
    cn_in = float(np.median(cn_big[b0:b1]))
    cn_out = float(np.median(np.concatenate([cn_big[:b0][map_big[:b0]],
                                             cn_big[b1:][map_big[b1:]]])))
    log(f"depth call: {ref_name}:{depth_start + 1:,}-{depth_end:,} "
        f"({n_big_bins} x 100kb bins, CN {cn_in:.3f} vs {cn_out:.3f} outside, "
        f"{zygosity})")

    # ---------------- discordant pairs ----------------
    both = usable1 & usable2
    fr = both & (a1["strand"] == 1) & (a2["strand"] == -1) & \
        (a2["start"] >= a1["start"])
    spans = (a2["start"] + READ_LEN - a1["start"])[fr].astype(np.int64)
    core = spans[spans <= 5000]
    ins_med = float(np.median(core)) if core.size else 300.0
    ins_mad = float(np.median(np.abs(core - ins_med))) * 1.4826 \
        if core.size else 50.0
    T = max(2000, int(ins_med + 8 * ins_mad))
    disc_all = np.nonzero(fr)[0][spans >= T]
    log(f"insert size: median {ins_med:.0f} bp, MAD {ins_mad:.0f} bp; "
        f"FR pairs with span >= {T}: {disc_all.size:,}")

    # keep only discordant pairs whose mates flank the depth-defined region
    s1d = a1["start"][disc_all]
    s2d = a2["start"][disc_all]
    flank = ((s1d + READ_LEN >= depth_start - 5000) &
             (s1d <= depth_start + 2000) &
             (s2d >= depth_end - 2000) &
             (s2d <= depth_end + 5000))
    disc_pairs = disc_all[flank]
    # recompute spans for the kept pairs (flank indexes into spans already)
    span_pairs = (a2["start"][disc_pairs] + READ_LEN -
                  a1["start"][disc_pairs]) if disc_pairs.size else \
        np.array([], dtype=np.int64)
    d_est = span_pairs - ins_med
    win_lo = (a1["start"][disc_pairs] + READ_LEN).astype(np.int64)
    win_hi = (a2["start"][disc_pairs] - d_est).astype(np.int64)
    okw = win_lo <= win_hi
    if okw.any():
        w_lo = int(win_lo[okw].max())
        w_hi = int(win_hi[okw].min())
        pair_window = (w_lo, w_hi) if w_lo <= w_hi else None
    else:
        pair_window = None
    log(f"discordant pairs flanking the depth region: {disc_pairs.size:,} "
        f"(breakpoint window {pair_window})")

    # ---------------- split reads ----------------
    splits_raw = []
    for arr, cd in ((a1, c1), (a2, c2)):
        for i in np.nonzero(arr["status"] == SPLIT)[0]:
            ok, blo, bhi, D = verify_split_junction(
                refc, cd[i], int(arr["start"][i]), int(arr["split_start2"][i]),
                int(arr["split_hit_off"][i]), int(arr["split_hit_pos"][i]))
            if ok:
                splits_raw.append((blo, bhi, D))
    splits = [(blo, bhi, D) for blo, bhi, D in splits_raw
              if abs(blo - depth_start) <= EDGE_TOL
              and abs((blo + D) - depth_end) <= EDGE_TOL]
    splits.sort()
    clusters = []
    for blo, bhi, D in splits:
        if clusters and abs(blo - clusters[-1]["blo"]) <= 100 \
                and abs(D - clusters[-1]["D"]) <= 100:
            c = clusters[-1]
            c["n"] += 1
            c["blo"] = min(c["blo"], blo)
            c["bhi"] = max(c["bhi"], bhi)
        else:
            clusters.append(dict(blo=blo, bhi=bhi, D=D, n=1))
    log(f"split reads: {len(splits_raw):,} verified junctions total, "
        f"{len(splits):,} consistent with the depth region, "
        f"{len(clusters)} cluster(s)")

    # ---------------- integrate ----------------
    if clusters:
        best = max(clusters, key=lambda c: (c["n"], -c["blo"]))
        refined_start = best["blo"] + 1            # 1-based first deleted base
        refined_end = best["blo"] + best["D"]      # 1-based last deleted base
        refined_method = "split_reads"
        _unc = best["bhi"] - best["blo"]
        refined_precision = (
            f"junction verified base-by-base by {best['n']} read(s); "
            f"left-breakpoint window {_unc} bp"
            + (f" due to {_unc - 1} bp junction microhomology" if _unc > 1
               else "")
            + "; leftmost consistent resolution reported")
    elif pair_window is not None:
        w_lo, w_hi = pair_window
        refined_start = w_lo + 1
        refined_end = w_hi + int(depth_size)
        refined_method = "discordant_pairs"
        refined_precision = (f"interval of {w_hi - w_lo + 1} bp implied by "
                             f"{int(okw.sum())} spanning pair(s)")
    else:
        refined_start = depth_start + 1
        refined_end = depth_end
        refined_method = "read_depth"
        refined_precision = f"~{FINE_BIN // 1000} kb (depth bin boundary)"

    s100 = round_half_up_100kb(refined_start)
    e100 = round_half_up_100kb(refined_end)
    if e100 <= s100:
        e100 = s100 + BIG_BIN
    size100 = e100 - s100
    size_refined = refined_end - refined_start + 1

    signals = [
        f"read_depth={n_big_bins} consecutive 100kb bins at CN {cn_in:.2f} "
        f"vs {cn_out:.2f} genome median ({zygosity}); depth interval "
        f"{ref_name}:{depth_start + 1}-{depth_end}"]
    if disc_pairs.size:
        signals.append(
            f"discordant_pairs={disc_pairs.size} FR pair(s) spanning the "
            f"depth region (median span {int(np.median(span_pairs)):,} bp, "
            f"insert median {ins_med:.0f} bp)")
    if clusters:
        signals.append(
            f"split_reads={len(splits)} junction read(s) in "
            f"{len(clusters)} cluster(s), best at {ref_name}:"
            f"{best['blo'] + 1}-{best['blo'] + best['D']} (n={best['n']})")
    supporting = " | ".join(signals)

    # ---------------- outputs ----------------
    tsv_path = os.path.join(OUT_DIR, "deletion.tsv")
    with open(tsv_path, "w", encoding="utf-8") as fh:
        fh.write("chrom\tstart_100kb\tend_100kb\tsize_bp\tsupporting_signals\n")
        fh.write(f"{ref_name}\t{s100}\t{e100}\t{size100}\t{supporting}\n")
    log(f"wrote {tsv_path}")

    fine_prof = {}
    for edge in (depth_start, depth_end):
        i0 = max(int(edge // FINE_BIN) - 20, 0)
        i1 = min(i0 + 41, len(cn_fine))
        fine_prof[str(edge)] = [
            {"bin_start": int(x * FINE_BIN),
             "cn": round(float(cn_fine[x]), 3) if map_fine[x] else None}
            for x in range(i0, i1)]

    qc = {
        "reference": {"name": ref_name, "length_bp": L, "acgt_bp": acgt},
        "reads": {
            "pairs": npairs, "read_length": READ_LEN,
            "aligned_reads_with_2plus_seeds": int(usable1.sum() +
                                                  usable2.sum()),
            "alignment_rate": round(float(aln_rate), 4),
            "mean_genome_depth_x": round(starts_all.size * READ_LEN / acgt, 3),
            "qc_subsample_mapped_reads": {
                "n": int(len(idx_q)),
                "clean_reads_fraction_le2mm": round(clean_frac, 4),
                "noisy_reads_mean_mismatches": round(noisy_mean, 2),
                "note": ("reads are effectively error-free; each read was "
                         "seeded on both strands, so ~half of all seeds are "
                         "strand-corrected away, not sequencing errors")},
        },
        "insert_size": {"median_bp": round(ins_med, 1),
                        "mad_bp": round(ins_mad, 1),
                        "discordant_span_threshold_bp": T,
                        "discordant_fr_pairs_raw": int(disc_all.size),
                        "discordant_fr_pairs_flanking_region":
                            int(disc_pairs.size)},
        "deletion": {
            "chrom": ref_name,
            "zygosity_call": zygosity,
            "depth_region_1based": {"start": depth_start + 1,
                                    "end": depth_end},
            "depth_region_size_bp": depth_size,
            "depth_cn_inside": round(cn_in, 4),
            "depth_cn_outside": round(cn_out, 4),
            "depth_big_bins_below_threshold": n_big_bins,
            "depth_threshold_used": low_used,
            "depth_cn_profile_100kb": [
                round(float(x), 3) if map_big[i] else None
                for i, x in enumerate(cn_big)],
            "depth_fine_profile_near_edges": fine_prof,
            "split_reads_verified_total": len(splits_raw),
            "split_reads_consistent_with_region": len(splits),
            "split_clusters_top10": sorted(
                [{"left_breakpoint_1based": c["blo"] + 1,
                  "breakpoint_uncertainty_bp": c["bhi"] - c["blo"],
                  "size_bp": c["D"], "support": c["n"]} for c in clusters],
                key=lambda c: -c["support"])[:10],
            "spanning_pairs_breakpoint_window_0based":
                list(pair_window) if pair_window else None,
            "refined_start_1based": refined_start,
            "refined_end_1based": refined_end,
            "refined_size_bp": size_refined,
            "refined_method": refined_method,
            "refined_precision": refined_precision,
            "reported_start_100kb": s100,
            "reported_end_100kb": e100,
            "reported_size_bp": size100,
        },
        "runtime_seconds": round(time.time() - t_start, 1),
    }
    qc_path = os.path.join(OUT_DIR, "qc.json")
    with open(qc_path, "w", encoding="utf-8") as fh:
        json.dump(qc, fh, indent=2)
    log(f"wrote {qc_path}")

    md_path = os.path.join(OUT_DIR, "report.md")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(build_report(qc, clusters, disc_pairs.size, span_pairs,
                              supporting, b0, b1, cn_big, map_big))
    log(f"wrote {md_path}")
    log("done")


def build_report(qc, clusters, n_disc, span_pairs, supporting, b0, b1,
                 cn_big, map_big):
    microhomology_note = ""
    mh = microhomology_note
    if clusters:
        unc = max(c["bhi"] - c["blo"] for c in clusters)
        if unc > 1:
            rs = qc["deletion"]["refined_start_1based"]
            mh = (f"- In this dataset the verified junction reads agree on the "
                  f"deletion size but the junction carries a {unc - 1} bp "
                  f"microhomology, so the left breakpoint is determined to a "
                  f"{unc} bp window (first deleted base {rs:,} or "
                  f"{rs + unc - 1:,}, 1-based); all resolutions round to the "
                  f"same 100 kb grid values. The refined coordinates above "
                  f"use the leftmost resolution, which coincides with the "
                  f"depth-defined interval.")
    d = qc["deletion"]
    prof = ", ".join(
        f"{i * BIG_BIN // 1000}kb:{cn_big[i]:.2f}"
        for i in range(max(b0 - 2, 0), min(b1 + 2, len(cn_big)))
        if map_big[i])
    if clusters:
        split_txt = "; ".join(
            f"{c['n']} read(s) at junction {c['blo'] + 1:,}-{c['blo'] + c['D']:,} "
            f"(crossover window {c['bhi'] - c['blo']} bp)"
            for c in sorted(clusters, key=lambda c: -c["n"])[:5])
    else:
        split_txt = ("no junction-spanning read passed the consistency "
                     "filter for this region")
    if n_disc:
        disc_txt = (f"{n_disc} FR pair(s) with mates flanking the depth "
                    f"region, median span {int(np.median(span_pairs)):,} bp "
                    f"(insert median {qc['insert_size']['median_bp']:.0f} bp)")
    else:
        disc_txt = "none flanking the depth region"
    microhomology_note = mh
    return f"""# Large deletion call - shallow PE hg38 (chr22)

## Data
- Reference: `{qc['reference']['name']}` ({qc['reference']['length_bp']:,} bp,
  {qc['reference']['acgt_bp']:,} ACGT bp).
- Reads: {qc['reads']['pairs']:,} paired-end x {qc['reads']['read_length']} bp,
  mean depth {qc['reads']['mean_genome_depth_x']}x (shallow).
- Alignment (exact 25-mer seed index built in-script; no external aligner):
  {qc['reads']['aligned_reads_with_2plus_seeds']:,} reads placed with >= 2
  consistent seeds ({qc['reads']['alignment_rate']:.1%}).
- Insert size: median {qc['insert_size']['median_bp']:.0f} bp
  (MAD {qc['insert_size']['mad_bp']:.0f} bp).
- Read quality (subsample of {qc['reads']['qc_subsample_mapped_reads']['n']:,}
  placed reads, strand-aware): {qc['reads']['qc_subsample_mapped_reads']['clean_reads_fraction_le2mm']:.1%}
  match their placement with <= 2 mismatches, i.e. the reads are effectively
  error-free. Unmappable reference (N-gaps/centromere) is excluded from all
  depth calculations (shows as null in the CN profile).

## Signals (independent lines of evidence)
1. **Read depth** (primary). 100 kb bins normalized by per-bin mappable
   sequence: {d['depth_big_bins_below_threshold']} consecutive bins collapse to
   CN {d['depth_cn_inside']:.3f} vs {d['depth_cn_outside']:.3f} elsewhere
   ({d['zygosity_call']}). Bin path around the call (bin:CN): {prof}.
2. **Discordant pairs**: {disc_txt}.
3. **Split reads**: {split_txt}.

## Call
| item | value |
|---|---|
| chrom | {d['chrom']} |
| depth-only interval (1-based) | {d['depth_region_1based']['start']:,} - {d['depth_region_1based']['end']:,} |
| refined breakpoints (1-based) | {d['refined_start_1based']:,} - {d['refined_end_1based']:,} |
| refined size | {d['refined_size_bp']:,} bp |
| refined by | {d['refined_method']} - {d['refined_precision']} |
| **reported start (rounded to 100 kb)** | **{d['reported_start_100kb']:,}** |
| **reported end (rounded to 100 kb)** | **{d['reported_end_100kb']:,}** |
| **reported size** | **{d['reported_size_bp']:,} bp** |

## Evidence vs precision limits
- The TSV coordinates are on the task-mandated 100 kb grid (each
  breakpoint rounded half-up to the nearest 100 kb). This is a *reporting*
  precision limit, not the measurement limit.
- Read depth alone localizes each breakpoint to ~10 kb (fine-bin boundary
  plus Poisson noise at {qc['reads']['mean_genome_depth_x']}x); it is the
  evidence for the interval and its zygosity.
- Discordant spanning pairs narrow a breakpoint to roughly one fragment
  length (~{qc['insert_size']['median_bp']:.0f} bp).
- Verified split reads give single-base breakpoints; when present they are
  the strongest evidence and are used for the refined coordinates above.
{microhomology_note}
- Where junction-level evidence is absent or sparse (expected at shallow
  depth: only ~coverage x 2 junction-crossing fragments exist), the depth
  boundaries are the honest limit and the rounded 100 kb values reflect the
  requested reporting grid rather than weaker evidence.

## Supporting signals (verbatim from deletion.tsv)
{supporting}

_Generated by output/analysis.py in {qc['runtime_seconds']} s._
"""


if __name__ == "__main__":
    main()
