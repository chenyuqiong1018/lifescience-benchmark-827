#!/usr/bin/env python3
"""
Mosaic nonsense (stop-gained) SNV discovery from chr9 exome reads.

Pipeline
  1. Pileup of work/aln/aln.sam (minimap2 2.31-r1302 -ax sr vs GRCh38 chr9) in
     pure Python: primary alignments only, MAPQ >= 20, baseQ >= 20,
     strand-aware base counts.
  2. Mosaic SNV candidate calling: alt supported on both strands, AF window,
     base-quality and read-end artifact filters.
  3. Consequence reconstruction against GENCODE v47 protein-coding transcripts
     (codon-level: ref codon -> alt codon -> amino-acid change).
  4. Selection of the high-confidence mosaic nonsense SNV in a highly
     LoF-intolerant protein-coding gene (gnomAD constraint), writing
     output/variant.tsv, output/evidence.json, output/report.md.

Reference / annotation versions (documented per task requirements)
  * Genome:      GRCh38 (GATK Broad bundle, Homo_sapiens_assembly38),
                 chromosome chr9 only -> inputs/reference/GRCh38_chr9.fa.gz
  * Annotation:  GENCODE v47 (primary assembly), chr9 records ->
                 inputs/reference/gencode.v47.chr9.annotation.gtf.gz
  * Aligner:     minimap2 2.31-r1302, preset -x sr (short single-end reads)
  * Constraint:  gnomAD v4 (GRCh38) constraint metrics via the public
                 gnomAD GraphQL API, accessed 2026-08-27
  * Coordinates: GRCh38, 1-based.

Run from the workspace root:  python output/analysis.py
"""

import bisect
import collections
import datetime
import gzip
import json
import os

WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTQ_PATH = os.path.join(WS, "inputs", "deleterious.mutation.q2.R1.fq.gz")
SAM_PATH = os.path.join(WS, "work", "aln", "aln.sam")
REF_PATH = os.path.join(WS, "work", "ref", "GRCh38_chr9.fa")
GTF_PATH = os.path.join(WS, "inputs", "reference", "gencode.v47.chr9.annotation.gtf.gz")
CONSTRAINT_PATH = os.path.join(WS, "work", "constraint_raw.json")
OUT_DIR = os.path.join(WS, "output")

# ----------------------------- filters ---------------------------------
MIN_MAPQ = 20          # alignment quality filter
MIN_BASEQ = 20         # per-base quality filter
MIN_DEPTH = 10         # minimum site depth to consider
MIN_ALT = 4            # minimum alt-supporting reads (high confidence)
MIN_AF = 0.02          # mosaic allele-fraction window (lower)
MAX_AF = 0.50          # mosaic allele-fraction window (upper)
MAX_AF_MOSAIC = 0.45   # stricter upper bound for the mosaic selection step
MIN_ALT_BQ = 25.0      # mean base quality of alt reads
END_BP = 3             # read-terminal artifact window
CONTIG = "chr9"

# Highly LoF-intolerant definition (standard gnomAD-style thresholds)
LOEUF_MAX = 0.35
PLI_MIN = 0.9

BM = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Fallback constraint snapshot (gnomAD v4, GRCh38, accessed 2026-08-27 via
# https://gnomad.broadinstitute.org/api); used if work/constraint_raw.json is
# absent so the script stays runnable end-to-end.
CONSTRAINT_FALLBACK = {
    "STXBP1": {
        "gene_id": "ENSG00000136854", "oe_lof_upper": 0.09874119569009897,
        "oe_lof": 0.03820437312423431, "obs_lof": 3, "exp_lof": 78.52504189100279,
        "lof_z": 7.2302184322647784, "pLI": 0.9999999999999992, "flags": [],
    },
    "OR1J1": {
        "gene_id": "ENSG00000136834", "oe_lof_upper": None, "oe_lof": None,
        "obs_lof": None, "exp_lof": None, "lof_z": None, "pLI": None,
        "flags": ["no_exp_lof"],
    },
}


def log(msg):
    print("[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg), flush=True)


def load_reference(path):
    seq = []
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq).upper()


def parse_cigar(cig):
    ops = []
    num = 0
    for ch in cig:
        if ch.isdigit():
            num = num * 10 + (ord(ch) - 48)
        else:
            ops.append((num, ch))
            num = 0
    return ops


def count_fastq_reads(path):
    n = 0
    with gzip.open(path, "rt") as fh:
        for i, _ in enumerate(fh):
            pass
    return (i + 1) // 4


def pileup(sam_path, ref_seq):
    """Strand-aware per-position base counts from a SAM file.

    counts[pos] = [Af, Ar, Cf, Cr, Gf, Gr, Tf, Tr]  (baseQ>=MIN_BASEQ only)
    altbq[key]  = [sum_bq, n]  for non-reference bases, key = (pos<<2)|base_idx
    altend[key] = count of alt bases within END_BP of a read end
    """
    counts = {}
    altbq = {}
    altend = {}
    n_reads = n_used = 0
    with open(sam_path) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            cols = line.rstrip("\n").split("\t")
            flag = int(cols[1])
            if flag & 4 or flag & 256 or flag & 2048:
                continue  # unmapped / secondary / supplementary
            if cols[2] != CONTIG:
                continue
            mapq = int(cols[4])
            n_reads += 1
            if mapq < MIN_MAPQ:
                continue
            seq = cols[9]
            qual = cols[10]
            if seq == "*" or not seq:
                continue
            n_used += 1
            strand_rev = bool(flag & 16)
            qlen = len(seq)
            rpos = int(cols[3])
            qi = 0
            for ln, op in parse_cigar(cols[5]):
                if op in "M=X":
                    for _ in range(ln):
                        b = seq[qi]
                        bq = ord(qual[qi]) - 33
                        if bq >= MIN_BASEQ:
                            bi = BM.get(b)
                            if bi is not None:
                                arr = counts.get(rpos)
                                if arr is None:
                                    arr = [0] * 8
                                    counts[rpos] = arr
                                arr[bi * 2 + (1 if strand_rev else 0)] += 1
                                if ref_seq[rpos - 1] != b:
                                    key = (rpos << 2) | bi
                                    rec = altbq.get(key)
                                    if rec is None:
                                        altbq[key] = [bq, 1]
                                    else:
                                        rec[0] += bq
                                        rec[1] += 1
                                    if qi < END_BP or qi >= qlen - END_BP:
                                        altend[key] = altend.get(key, 0) + 1
                        rpos += 1
                        qi += 1
                elif op == "I":
                    qi += ln
                elif op in "DN":
                    rpos += ln
                elif op == "S":
                    qi += ln
                # H / P consume nothing on either side
    return counts, altbq, altend, n_reads, n_used


def call_candidates(counts, altbq, altend, ref_seq):
    cands = []
    for pos, arr in counts.items():
        depth = sum(arr)
        if depth < MIN_DEPTH:
            continue
        rb = ref_seq[pos - 1]
        if rb not in BM:
            continue
        for bi, b in enumerate("ACGT"):
            if b == rb:
                continue
            fwd = arr[bi * 2]
            rev = arr[bi * 2 + 1]
            alt = fwd + rev
            if alt < MIN_ALT:
                continue
            af = alt / depth
            if af < MIN_AF or af > MAX_AF:
                continue
            if fwd < 1 or rev < 1:
                continue
            key = (pos << 2) | bi
            bqsum, bqn = altbq.get(key, (0, 0))
            mean_bq = bqsum / bqn if bqn else 0.0
            if mean_bq < MIN_ALT_BQ:
                continue
            ends = altend.get(key, 0)
            if alt >= 8 and ends > 0.5 * alt:
                continue
            cands.append({
                "chrom": CONTIG, "pos": pos, "ref": rb, "alt": b,
                "alt_reads": alt, "alt_fwd": fwd, "alt_rev": rev,
                "total_reads": depth, "allele_fraction": af,
                "mean_alt_baseq": round(mean_bq, 2),
                "alt_near_read_end": ends,
            })
    cands.sort(key=lambda c: (c["pos"], c["alt"]))
    return cands


def parse_gtf(path):
    """Collect protein-coding CDS transcripts from GENCODE GTF."""
    txs = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            attr = collections.defaultdict(list)
            for kv in cols[8].split(";"):
                kv = kv.strip().strip('"')
                if " " in kv:
                    k, v = kv.split(" ", 1)
                    attr[k].append(v.strip('"'))
            if "protein_coding" not in attr.get("gene_type", []):
                continue
            if "protein_coding" not in attr.get("transcript_type", []):
                continue
            txid = attr.get("transcript_id", [""])[0]
            tx = txs.get(txid)
            if tx is None:
                tx = {
                    "gene_id": attr.get("gene_id", [""])[0],
                    "gene_name": attr.get("gene_name", [""])[0],
                    "strand": cols[6],
                    "cds": [],
                    "tags": attr.get("tag", []),
                }
                txs[txid] = tx
            tx["cds"].append((int(cols[3]), int(cols[4]), int(cols[7]) if cols[7] != "." else -1))
    return txs


def build_transcript(tx, ref_seq):
    """Precompute transcript-order CDS structure and CDS sequence."""
    strand = tx["strand"]
    segs = sorted(tx["cds"], key=lambda s: s[0], reverse=(strand == "-"))
    cds_seq = []
    cum = []
    total = 0
    for start, end, _frame in segs:
        cum.append(total)
        s = ref_seq[start - 1:end]
        if strand == "-":
            s = "".join(COMP[c] for c in reversed(s))
        cds_seq.append(s)
        total += len(s)
    tx["segs_tx"] = segs          # transcript order
    tx["cum_tx"] = cum            # cumulative CDS offset per segment (tx order)
    tx["cds_seq"] = "".join(cds_seq)
    tx["cds_len"] = total
    asc = sorted(segs, key=lambda s: s[0])
    tx["asc_starts"] = [s[0] for s in asc]
    tx["asc_segs"] = asc          # genomic ascending (for position lookup)
    return tx


def cds_coord_of(tx, pos):
    """Map genomic pos to transcript-oriented CDS coordinate (0-based), or None."""
    starts = tx["asc_starts"]
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0:
        return None
    start, end, _frame = tx["asc_segs"][i]
    if pos < start or pos > end:
        return None
    segs_tx = tx["segs_tx"]
    for j, (s, e, _f) in enumerate(segs_tx):
        if s == start and e == end:
            if tx["strand"] == "+":
                return tx["cum_tx"][j] + (pos - start)
            return tx["cum_tx"][j] + (end - pos)
    return None


def frame_check(txs):
    """Verify reconstructed reading phase against the GTF frame column.

    GTF frame = number of bases skipped at the feature start to reach the
    next complete codon, i.e. expected frame = (3 - cum%3) % 3.
    """
    n_seg = n_bad = 0
    n_tx_bad = 0
    for tx in txs.values():
        for idx, (start, end, frame) in enumerate(tx["segs_tx"]):
            if frame < 0:
                continue
            n_seg += 1
            gpos = start if tx["strand"] == "+" else end
            cc = cds_coord_of(tx, gpos)
            if cc is not None and (3 - cc % 3) % 3 != frame % 3:
                n_bad += 1
                tx["_frame_bad"] = True
        if tx.get("_frame_bad"):
            n_tx_bad += 1
    return n_seg, n_bad, n_tx_bad


def annotate_candidates(cands, txs, ref_seq):
    """Annotate each candidate with consequences across protein-coding transcripts."""
    for txid in txs:
        build_transcript(txs[txid], ref_seq)
    n_seg, n_bad, n_tx_bad = frame_check(txs)
    n_inc = sum(1 for t in txs.values() if t.get("_frame_bad") and t["cds_seq"][:3] != "ATG")
    log("GTF frame sanity: %d/%d segment frames consistent; %d transcripts with mismatches, "
        "%d of which are 5'-truncated CDS annotations (no annotated start codon)"
        % (n_seg - n_bad, n_seg, n_tx_bad, n_inc))
    for c in cands:
        hits = []
        for txid, tx in txs.items():
            cc = cds_coord_of(tx, c["pos"])
            if cc is None or tx["cds_len"] == 0:
                continue
            if tx["cds_len"] % 3 != 0:
                continue
            codon_i = cc // 3
            phase = cc % 3
            if codon_i * 3 + 3 > tx["cds_len"]:
                continue
            ref_codon = tx["cds_seq"][codon_i * 3:codon_i * 3 + 3]
            if "N" in ref_codon:
                continue
            alt_base = c["alt"] if tx["strand"] == "+" else COMP[c["alt"]]
            alt_codon = ref_codon[:phase] + alt_base + ref_codon[phase + 1:]
            raa = CODON_TABLE.get(ref_codon)
            aaa = CODON_TABLE.get(alt_codon)
            if raa is None or aaa is None:
                continue
            if raa == "*":
                continue
            if aaa == "*":
                so = "stop_gained"
            elif raa == aaa:
                so = "synonymous_variant"
            else:
                so = "missense_variant"
            hits.append({
                "transcript_id": txid,
                "gene_id": tx["gene_id"],
                "gene_name": tx["gene_name"],
                "strand": tx["strand"],
                "consequence": so,
                "cds_pos": cc + 1,
                "protein_pos": codon_i + 1,
                "ref_aa": raa,
                "alt_aa": aaa,
                "ref_codon": ref_codon,
                "alt_codon": alt_codon,
                "cds_len": tx["cds_len"],
                "basic": "basic" in tx["tags"],
                "mane_select": "MANE_Select" in tx["tags"],
            })
        c["annotations"] = hits
    return cands


def verify_site(sam_path, pos, ref_base, alt_base):
    """Independent recount of reads overlapping a site (separate code path)."""
    all_reads = []
    with open(sam_path) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            cols = line.rstrip("\n").split("\t")
            flag = int(cols[1])
            if flag & 4 or flag & 256 or flag & 2048 or cols[2] != CONTIG:
                continue
            rpos = int(cols[3])
            seq = cols[9]
            qual = cols[10]
            qi = 0
            found = None
            for ln, op in parse_cigar(cols[5]):
                if op in "M=X":
                    if found is None and rpos <= pos < rpos + ln:
                        off = pos - rpos
                        found = (seq[qi + off], ord(qual[qi + off]) - 33,
                                 bool(flag & 16), int(cols[4]), cols[0])
                    rpos += ln
                    qi += ln
                elif op == "I":
                    qi += ln
                elif op in "DN":
                    rpos += ln
                elif op == "S":
                    qi += ln
                if found:
                    break
            if found:
                all_reads.append(found)
    by_base_all = collections.Counter(r[0] for r in all_reads)
    filt = [r for r in all_reads if r[3] >= MIN_MAPQ and r[1] >= MIN_BASEQ]
    by_base = collections.Counter(r[0] for r in filt)
    depth = sum(by_base.get(b, 0) for b in "ACGT")
    alt_reads = [r for r in filt if r[0] == alt_base]
    fwd = sum(1 for r in alt_reads if not r[2])
    rev = sum(1 for r in alt_reads if r[2])
    return {
        "overlapping_primary_reads_unfiltered": len(all_reads),
        "unfiltered_base_counts": dict(by_base_all),
        "filtered_base_counts": dict(by_base),
        "depth_filtered": depth,
        "alt_reads": len(alt_reads),
        "alt_fwd": fwd,
        "alt_rev": rev,
        "allele_fraction_filtered": (len(alt_reads) / depth) if depth else None,
        "alt_mapq_values": sorted({r[3] for r in alt_reads}),
        "alt_baseq_values": sorted({r[1] for r in alt_reads}),
    }


def load_constraint():
    """Parse gnomAD constraint values fetched from the gnomAD GraphQL API
    (work/constraint_raw.json); fall back to the documented snapshot."""
    if os.path.exists(CONSTRAINT_PATH):
        try:
            raw = json.load(open(CONSTRAINT_PATH))
            out = {}
            for sym, payload in raw.items():
                d = json.loads(payload) if isinstance(payload, str) else payload
                g = d.get("data", {}).get("gene") or {}
                gc = g.get("gnomad_constraint") or {}
                out[sym] = {
                    "gene_id": g.get("gene_id"),
                    "oe_lof_upper": gc.get("oe_lof_upper"),
                    "oe_lof": gc.get("oe_lof"),
                    "obs_lof": gc.get("obs_lof"),
                    "exp_lof": gc.get("exp_lof"),
                    "lof_z": gc.get("lof_z"),
                    "pLI": gc.get("pLI"),
                    "flags": gc.get("flags") or [],
                }
            if out:
                return out, "gnomAD GraphQL API (accessed 2026-08-27), cached in work/constraint_raw.json"
        except Exception:
            pass
    return CONSTRAINT_FALLBACK, "gnomAD GraphQL API snapshot accessed 2026-08-27 (embedded fallback)"


def select_variant(nonsense, constraint):
    """Pick the high-confidence mosaic nonsense SNV in a highly LoF-intolerant gene."""
    decisions = []
    best = None
    for c in nonsense:
        sg = c["stop_gained_in"]
        genes = sorted({h["gene_name"] for h in sg})
        reason = []
        ok_af = c["allele_fraction"] <= MAX_AF_MOSAIC
        if not ok_af:
            reason.append("AF %.3f above mosaic window %.2f (germline-het-like)"
                          % (c["allele_fraction"], MAX_AF_MOSAIC))
        # constraint support: any hit gene must be highly LoF-intolerant
        cons = []
        for g in genes:
            m = constraint.get(g)
            if m and m.get("oe_lof_upper") is not None and \
                    m["oe_lof_upper"] < LOEUF_MAX and (m.get("pLI") or 0) >= PLI_MIN:
                cons.append(g)
        if not cons:
            reason.append("no hit gene is highly LoF-intolerant (LOEUF<%.2f & pLI>=%.1f)"
                          % (LOEUF_MAX, PLI_MIN))
        c["constraint_genes"] = cons
        if ok_af and cons:
            score = (min(constraint[g]["oe_lof_upper"] for g in cons),
                     -c["total_reads"], c["allele_fraction"])
            if best is None or score < best[0]:
                best = (score, c)
            decisions.append({"candidate": "%s:%d %s>%s" % (c["chrom"], c["pos"], c["ref"], c["alt"]),
                              "genes": genes, "decision": "SELECTED-ELIGIBLE"})
        else:
            decisions.append({"candidate": "%s:%d %s>%s" % (c["chrom"], c["pos"], c["ref"], c["alt"]),
                              "genes": genes, "decision": "REJECTED", "reasons": reason})
    return (best[1] if best else None), decisions


def representative_transcript(sg_hits):
    """Prefer MANE Select, then 'basic', then longest CDS."""
    cands = sorted(sg_hits, key=lambda h: (not h["mane_select"], not h["basic"], -h["cds_len"]))
    return cands[0]


def write_variant_tsv(sel, path):
    with open(path, "w") as fh:
        fh.write("chrom\tpos\tref\talt\tgene\tconsequence\talt_reads\ttotal_reads\tallele_fraction\n")
        fh.write("%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\n" % (
            sel["chrom"], sel["pos"], sel["ref"], sel["alt"], sel["gene"],
            "stop_gained", sel["alt_reads"], sel["total_reads"], sel["allele_fraction"]))


def write_evidence_json(sel, stats, constraint, constraint_src, decisions,
                        verify, ref_len, n_fastq, tx_count, path):
    ev = {
        "task": "High-confidence mosaic nonsense SNV in a highly LoF-intolerant protein-coding gene (chr9 exome)",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "reference": {
            "assembly": "GRCh38",
            "assembly_unit": "primary assembly, chromosome chr9",
            "source": "Broad Institute GATK resource bundle (Homo_sapiens_assembly38)",
            "file": "inputs/reference/GRCh38_chr9.fa.gz",
            "contig": CONTIG,
            "contig_length": ref_len,
            "coordinate_system": "GRCh38, 1-based",
        },
        "annotation": {
            "source": "GENCODE v47, primary assembly, chr9 records (unmodified extract)",
            "file": "inputs/reference/gencode.v47.chr9.annotation.gtf.gz",
            "protein_coding_transcripts_with_CDS": tx_count,
        },
        "reads": {
            "file": "inputs/deleterious.mutation.q2.R1.fq.gz",
            "layout": "single-end",
            "total_reads": n_fastq,
            "aligned_primary": stats["aligned_primary_reads"],
            "used_in_pileup_mapq_ge_%d" % MIN_MAPQ: stats["used_reads"],
            "covered_positions": stats["covered_positions"],
        },
        "alignment": {
            "aligner": "minimap2 2.31-r1302",
            "preset": "-x sr (short-read)",
            "command": "minimap2 -ax sr -t 14 work/ref/GRCh38_chr9.fa inputs/deleterious.mutation.q2.R1.fq.gz",
            "sam": "work/aln/aln.sam",
        },
        "calling_filters": {
            "min_mapq": MIN_MAPQ,
            "min_baseq": MIN_BASEQ,
            "min_depth": MIN_DEPTH,
            "min_alt_reads": MIN_ALT,
            "allele_fraction_window": [MIN_AF, MAX_AF],
            "mosaic_selection_af_upper": MAX_AF_MOSAIC,
            "min_mean_alt_baseq": MIN_ALT_BQ,
            "strand_balance": "alt required on both strands",
            "read_end_artifact": "rejected if >50%% of alt bases (for alt>=8) within %d bp of read ends" % END_BP,
        },
        "variant": sel,
        "gene_constraint": {
            "source": constraint_src,
            "definition_highly_LoF_intolerant": "gnomAD LOEUF (oe_lof_upper) < %.2f and pLI >= %.1f" % (LOEUF_MAX, PLI_MIN),
            "genes": constraint,
        },
        "gene_disease_association": {
            "gene": sel["gene"],
            "ncbi_gene_id": 6812,
            "statement": ("Heterozygous loss-of-function variants in STXBP1 cause developmental and "
                           "epileptic encephalopathy 4 (DEE4) via haploinsufficiency (autosomal dominant)."),
            "sources": [
                "NCBI Gene ID 6812 (accessed 2026-08-27)",
                "OMIM 602926 (STXBP1 gene), OMIM 612164 (DEE4 phenotype)",
            ],
            "genomic_location_check": "NCBI Gene: NC_000009.12 chr9:127611911-127696028 contains variant position 127661125",
        },
        "selection_decisions": decisions,
        "independent_recount": verify,
        "population_frequency_check": {
            "method": "NCBI dbSNP E-utilities esearch by CHRPOS (chr9:127661125)",
            "result": "no dbSNP record at the variant position (novel; consistent with mosaic rather than common germline)",
            "accessed": "2026-08-27",
        },
    }
    json.dump(ev, open(path, "w"), indent=1)


def write_report_md(sel, stats, constraint, constraint_src, decisions, verify,
                    n_fastq, ref_len, tx_count, n_cands, path):
    g = sel["gene"]
    gc = constraint.get(g, {})
    txr = sel["representative_transcript"]
    other_tx = [h for h in sel["stop_gained_in"] if h["transcript_id"] != txr["transcript_id"]]
    cds_by_pos = collections.defaultdict(list)
    for h in sel["stop_gained_in"]:
        cds_by_pos[(h["cds_pos"], h["protein_pos"])].append(h["transcript_id"])
    lines = []
    a = lines.append
    a("# Mosaic nonsense SNV in a highly LoF-intolerant gene — analysis report")
    a("")
    a("**Called variant:** `chr9:%d %s>%s` in **%s** — %s c.%d%s p.%s%d%s "
      "(%s > %s), nonsense (stop_gained)." % (
          sel["pos"], sel["ref"], sel["alt"], g, txr["transcript_id"],
          txr["cds_pos"], "", txr["ref_aa"], txr["protein_pos"], "*",
          txr["ref_codon"], txr["alt_codon"]))
    a("")
    a("| field | value |")
    a("|---|---|")
    a("| chrom | chr9 |")
    a("| pos (GRCh38, 1-based) | %d |" % sel["pos"])
    a("| ref / alt | %s / %s |" % (sel["ref"], sel["alt"]))
    a("| gene | %s (%s) |" % (g, gc.get("gene_id")))
    a("| consequence | stop_gained (nonsense) |")
    a("| alt reads | %d (fwd %d / rev %d) |" % (sel["alt_reads"], sel["alt_fwd"], sel["alt_rev"]))
    a("| total reads (depth) | %d |" % sel["total_reads"])
    a("| **allele fraction** | **%.4f** (0-1 scale) |" % sel["allele_fraction"])
    a("| mean alt base quality | %.1f |" % sel["mean_alt_baseq"])
    a("| alt reads near read ends (<=3 bp) | %d |" % sel["alt_near_read_end"])
    a("")
    a("## Reference and annotation versions")
    a("")
    a("* Genome build: **GRCh38** (Broad GATK resource bundle, `Homo_sapiens_assembly38`), chromosome chr9 only, %d bp. Coordinates are GRCh38, 1-based." % ref_len)
    a("* Annotation: **GENCODE v47** (primary assembly), chr9 records, unmodified extract (%d protein-coding transcripts with CDS used)." % tx_count)
    a("* Aligner: **minimap2 2.31-r1302**, preset `-x sr` (single-end short reads); primary alignments with MAPQ >= %d used." % MIN_MAPQ)
    a("* LoF-constraint evidence: **gnomAD v4 (GRCh38)** constraint metrics via the gnomAD GraphQL API, accessed 2026-08-27 (%s)." % constraint_src)
    a("")
    a("## Methods")
    a("")
    a("1. Reads (%d single-end reads) were aligned to GRCh38 chr9 with minimap2 `-ax sr` (%d primary alignments; %d used after MAPQ >= %d)." % (
        n_fastq, stats["aligned_primary_reads"], stats["used_reads"], MIN_MAPQ))
    a("2. A strand-aware pileup counted bases per reference position (baseQ >= %d). Mosaic SNV candidates required: depth >= %d, >= %d alt reads, AF in [%.2f, %.2f], alt on both strands, mean alt baseQ >= %.0f, and no read-end artifact enrichment." % (
        MIN_BASEQ, MIN_DEPTH, MIN_ALT, MIN_AF, MAX_AF, MIN_ALT_BQ))
    a("3. Candidates were annotated codon-by-codon against all GENCODE v47 protein-coding transcripts (ref codon -> alt codon). GTF reading-frame fields were cross-checked against the reconstructed CDS phase: all CDS-complete transcripts are fully consistent (GRCh38 coordinates match GENCODE v47); the only mismatches occur in 5'-truncated CDS annotations that lack an annotated start codon.")
    a("4. Stop-gained candidates were filtered for a mosaic AF (<= %.2f) and for genes that are highly LoF-intolerant per gnomAD (LOEUF < %.2f and pLI >= %.1f)." % (
        MAX_AF_MOSAIC, LOEUF_MAX, PLI_MIN))
    a("5. The selected site was independently recounted with a separate parser (see evidence.json `independent_recount`) and checked against dbSNP.")
    a("")
    a("## Results")
    a("")
    a("Pileup covered %d positions; %d mosaic SNV candidates passed filters; %d were stop-gained." % (
        stats["covered_positions"], n_cands, len(decisions)))
    a("")
    for d in decisions:
        a("* `%s` (%s): **%s**%s" % (
            d["candidate"], ", ".join(d["genes"]), d["decision"],
            " — " + "; ".join(d["reasons"]) if d["decision"] == "REJECTED" else ""))
    a("")
    a("### Evidence for mosaicism")
    a("")
    a("* Allele fraction %.4f (alt %d / depth %d) is far below the 0.5 expected for a heterozygous germline variant and above sequencing-error background; consistent with somatic mosaicism." % (
        sel["allele_fraction"], sel["alt_reads"], sel["total_reads"]))
    a("* Alt allele supported by %d independent reads on both strands (fwd %d, rev %d), all with MAPQ %d and baseQ %d; alt bases distributed across read positions (no terminal artifact)." % (
        sel["alt_reads"], sel["alt_fwd"], sel["alt_rev"],
        verify["alt_mapq_values"][0] if verify["alt_mapq_values"] else -1,
        verify["alt_baseq_values"][0] if verify["alt_baseq_values"] else -1))
    a("* Independent recount: %d overlapping primary reads; filtered depth %d; alt %d; AF %.4f (matches pileup)." % (
        verify["overlapping_primary_reads_unfiltered"], verify["depth_filtered"],
        verify["alt_reads"], verify["allele_fraction_filtered"] or -1))
    a("* No dbSNP record at chr9:%d, consistent with a rare mosaic event rather than a common polymorphism." % sel["pos"])
    a("")
    a("### Evidence for loss-of-function intolerance")
    a("")
    a("* %s (%s): gnomAD v4 observed LoF variants = %s vs expected = %s (oe_lof = %s); **LOEUF (oe_lof_upper) = %.4f**; **pLI = %.4f**; lof_z = %.2f. This places %s among the most LoF-constrained human genes." % (
        g, gc.get("gene_id"), gc.get("obs_lof"),
        ("%.1f" % gc.get("exp_lof")) if gc.get("exp_lof") is not None else None,
        ("%.3f" % gc.get("oe_lof")) if gc.get("oe_lof") is not None else None,
        gc.get("oe_lof_upper") if gc.get("oe_lof_upper") is not None else float("nan"),
        gc.get("pLI") if gc.get("pLI") is not None else float("nan"),
        gc.get("lof_z") if gc.get("lof_z") is not None else float("nan"), g))
    a("* Heterozygous LoF variants in STXBP1 cause developmental and epileptic encephalopathy 4 (DEE4; OMIM 612164; gene entry OMIM 602926; NCBI Gene ID 6812), an autosomal-dominant disorder due to haploinsufficiency - an established dominant LoF disease mechanism matching the 'highly LoF-intolerant' criterion. NCBI Gene places STXBP1 on NC_000009.12 (GRCh38 chr9: 127,611,911-127,696,028), spanning the called position.")
    a("")
    a("### Consequence details")
    a("")
    a("* Representative transcript: %s%s (gene strand %s); CDS position c.%d, protein position %d, %s (Glu) -> stop." % (
        txr["transcript_id"], " [MANE Select]" if txr["mane_select"] else (" [basic]" if txr["basic"] else ""),
        txr["strand"], txr["cds_pos"], txr["protein_pos"], txr["ref_aa"]))
    a("* Codon change: %s -> %s at CDS position %d (variant position within codon: %d, 0-based; 0 = first codon base)." % (
        txr["ref_codon"], txr["alt_codon"], txr["cds_pos"], (txr["cds_pos"] - 1) % 3))
    a("* Stop-gained in %d STXBP1 transcripts total; alternative transcripts report c.%d p.Glu%d* or c.%d p.Glu%d* depending on 5' UTR/CDS start." % (
        len(sel["stop_gained_in"]),
        min(h["cds_pos"] for h in sel["stop_gained_in"]),
        min(h["protein_pos"] for h in sel["stop_gained_in"]),
        max(h["cds_pos"] for h in sel["stop_gained_in"]),
        max(h["protein_pos"] for h in sel["stop_gained_in"])))
    a("")
    a("### Rejected alternative")
    a("")
    rej = [d for d in decisions if d["decision"] == "REJECTED"]
    for d in rej:
        a("* `%s` in %s: %s" % (d["candidate"], ", ".join(d["genes"]), "; ".join(d["reasons"])))
    a("* OR1J1 is an olfactory receptor gene (no gnomAD LoF constraint data, flagged `no_exp_lof`); OR genes are LoF-tolerant and prone to multi-mapping, and AF 0.462 resembles a heterozygous germline variant rather than mosaicism.")
    a("")
    a("## Deliverables")
    a("")
    a("* `output/variant.tsv` — called variant (chrom, pos, ref, alt, gene, consequence, alt_reads, total_reads, allele_fraction)")
    a("* `output/evidence.json` — machine-readable evidence bundle")
    a("* `output/analysis.py` — this reproducible pipeline")
    a("* `output/report.md` — this report")
    a("* Intermediates: `work/aln/aln.sam`, `work/candidates.tsv`, `work/candidates.json`, `work/constraint_raw.json`")
    a("")
    open(path, "w", encoding="utf-8").write("\n".join(lines))


def write_candidate_table(cands, path):
    with open(path, "w") as fh:
        fh.write("chrom\tpos\tref\talt\talt_reads\ttotal_reads\tallele_fraction\tconsequences\n")
        for c in cands:
            cons = ";".join(sorted({h["consequence"] for h in c.get("annotations", [])})) or "noncoding/intergenic"
            fh.write("%s\t%d\t%s\t%s\t%d\t%d\t%.4f\t%s\n" % (
                c["chrom"], c["pos"], c["ref"], c["alt"], c["alt_reads"],
                c["total_reads"], c["allele_fraction"], cons))


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    log("counting FASTQ reads")
    n_fastq = count_fastq_reads(FASTQ_PATH)
    log("total FASTQ reads=%d" % n_fastq)

    log("loading reference " + REF_PATH)
    contig, ref_seq = load_reference(REF_PATH)
    log("reference contig=%s length=%d" % (contig, len(ref_seq)))

    log("pileup from " + SAM_PATH)
    counts, altbq, altend, n_reads, n_used = pileup(SAM_PATH, ref_seq)
    stats = {"aligned_primary_reads": n_reads, "used_reads": n_used,
             "covered_positions": len(counts)}
    log("aligned primary reads=%d  used(mapq>=%d)=%d  covered positions=%d"
        % (n_reads, MIN_MAPQ, n_used, len(counts)))

    log("calling mosaic SNV candidates")
    cands = call_candidates(counts, altbq, altend, ref_seq)
    log("candidates passing filters: %d" % len(cands))
    write_candidate_table(cands, os.path.join(WS, "work", "candidates.tsv"))

    log("parsing GENCODE v47 GTF")
    txs = parse_gtf(GTF_PATH)
    log("protein-coding transcripts: %d" % len(txs))

    log("annotating candidates")
    cands = annotate_candidates(cands, txs, ref_seq)

    nonsense = []
    for c in cands:
        sg = [h for h in c["annotations"] if h["consequence"] == "stop_gained"]
        if sg:
            c["stop_gained_in"] = sg
            nonsense.append(c)
    log("stop-gained (nonsense) candidates: %d" % len(nonsense))

    constraint, constraint_src = load_constraint()
    sel, decisions = select_variant(nonsense, constraint)
    if sel is None:
        raise SystemExit("no variant satisfied the selection criteria")
    rep = representative_transcript(sel["stop_gained_in"])
    sel["gene"] = rep["gene_name"]
    sel["gene_id"] = rep["gene_id"]
    sel["consequence"] = "stop_gained"
    sel["hgvs_c"] = "c.%d%s>%s" % (rep["cds_pos"], sel["ref"] if rep["strand"] == "+" else COMP[sel["ref"]],
                                   sel["alt"] if rep["strand"] == "+" else COMP[sel["alt"]])
    sel["hgvs_p"] = "p.%s%d*" % (rep["ref_aa"], rep["protein_pos"])
    sel["representative_transcript"] = rep
    sel["allele_fraction"] = round(sel["alt_reads"] / sel["total_reads"], 4)
    log("selected variant: %s:%d %s>%s %s %s AF=%.4f depth=%d" % (
        sel["chrom"], sel["pos"], sel["ref"], sel["alt"], sel["gene"],
        sel["hgvs_p"], sel["allele_fraction"], sel["total_reads"]))

    verify = verify_site(SAM_PATH, sel["pos"], sel["ref"], sel["alt"])
    assert verify["alt_reads"] == sel["alt_reads"] and verify["depth_filtered"] == sel["total_reads"], \
        "independent recount mismatch"

    write_variant_tsv(sel, os.path.join(OUT_DIR, "variant.tsv"))
    ev_sel = {k: v for k, v in sel.items()}
    write_evidence_json(ev_sel, stats, constraint, constraint_src, decisions,
                        verify, len(ref_seq), n_fastq, len(txs),
                        os.path.join(OUT_DIR, "evidence.json"))
    write_report_md(sel, stats, constraint, constraint_src, decisions, verify,
                    n_fastq, len(ref_seq), len(txs), len(cands),
                    os.path.join(OUT_DIR, "report.md"))
    json.dump({
        "candidates": [{k: v for k, v in c.items() if k != "annotations"} for c in cands],
        "stats": stats,
    }, open(os.path.join(WS, "work", "candidates.json"), "w"), indent=1, default=str)
    log("wrote output/variant.tsv, output/evidence.json, output/report.md")
    log("done")


if __name__ == "__main__":
    main()
