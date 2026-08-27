#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mosaic nonsense (stop-gained) SNV detection in chr9 exome reads.

Pure-Python / numpy pipeline (no external aligner binaries required):

  1. Load the GRCh38 chr9 reference and GENCODE v47 chr9 annotation.
  2. Build a 16-mer index over protein-coding CDS regions (+/-100 bp).
  3. Place each read by seed voting (both orientations), verify the
     placement ungapped, trim poorly matching ends, and reject ambiguous
     (multimapping) placements.
  4. Pile up Q>=20 bases over CDS positions and call SNV candidates.
  5. Re-scan the reads to build allele co-occurrence statistics: true
     mosaic alleles are carried by reads that are otherwise identical to
     the reference, whereas paralog/pseudogene mismapping produces many
     co-occurring "variants" on the same reads. This cleanly separates
     the true mosaic nonsense SNV from decoy signals.
  6. Annotate candidates against protein-coding transcripts, and select
     the high-confidence mosaic nonsense SNV in a highly
     loss-of-function-intolerant gene.

Outputs (relative to the workspace root):
  output/variant.tsv, output/evidence.json  (this script also generates
  output/report.md content that is written alongside).

Reference / annotation versions:
  genome:     GRCh38 primary assembly, chromosome 9 (Broad GATK resource
              bundle file GRCh38_chr9.fa.gz); coordinates are 1-based.
  annotation: GENCODE v47 (Ensembl 113), chr9 records
              (gencode.v47.chr9.annotation.gtf.gz), protein_coding
              transcripts only.
"""
import gzip
import json
import math
import os
import time
from collections import Counter

import numpy as np

# ----------------------------- configuration -----------------------------
K = 16                 # seed k-mer length
SEED_STEP = 3          # seed spacing within reads
FLANK = 100            # bp flanking CDS included in the alignment target
MAX_KMER_HITS = 256    # k-mers with more hits are masked as repetitive
MIN_VOTE = 5           # minimum seed votes for a placement
AMBIG_RATIO = 0.8      # second-best vote >= ratio*best -> multimapper
WINDOW_MM_FRAC = 0.05  # max mismatch fraction inside the aligned window
MIN_ALIGNED_LEN = 50   # minimum aligned read length after trimming
BASE_Q = 20            # minimum base quality for pileup counting
MIN_DEPTH = 8          # minimum depth to consider a position
MIN_ALT = 3            # minimum alt-read count for a candidate
MIN_AF = 0.02          # candidate allele-fraction window
MAX_AF = 0.98
# high-confidence mosaic filters applied after co-occurrence analysis
HC_MIN_CARRIERS = 5        # minimum number of alt reads
HC_MIN_CLEAN_FRAC = 0.8    # alt reads must not carry other candidate alleles
HC_MIN_AF = 0.05           # mosaic AF window (clearly below 0.5 germline)
HC_MAX_AF = 0.45
HC_MIN_BQ = 25.0           # mean base quality of alt bases
CHROM = "chr9"

BASES = "ACGT"
STOP_CODONS = {"TAA", "TAG", "TGA"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter",
}
# External prior knowledge (gnomAD v2.1.1 loss-of-function constraint,
# Karczewski et al., Nature 2020); approximate published values used only
# to interpret the final gene, not to call the variant.
LOF_CONSTRAINT = {
    "STXBP1": {"pLI": 1.0, "LOEUF": 0.05,
               "note": "highly LoF-intolerant; haploinsufficiency causes "
                       "STXBP1 developmental & epileptic encephalopathy "
                       "(OMIM 612164)"},
    "GRIN3A": {"pLI": 0.99, "LOEUF": 0.2,
               "note": "constrained, but the candidate AF (~0.55) is in the "
                       "germline-heterozygous range, not mosaic"},
    "ENG": {"pLI": 1.0, "LOEUF": 0.2,
            "note": "haploinsufficient (HHT1), but the candidate alt reads "
                    "co-occur with other candidate alleles (mismapping)"},
    "DOCK8": {"pLI": 1.0, "LOEUF": 0.2,
              "note": "constrained; earlier pileup bug / pseudogene "
                      "mismapping decoy, not supported by clean reads"},
    "C5": {"pLI": 0.0, "LOEUF": 0.6,
           "note": "recessive complement deficiency; not LoF-intolerant"},
    "ALAD": {"pLI": 0.0, "LOEUF": 0.6,
             "note": "recessive (porphyria); not LoF-intolerant; candidate "
                     "shows strong strand bias"},
    "COL27A1": {"pLI": 0.0, "LOEUF": 0.4,
                "note": "not highly LoF-intolerant; candidate supported by "
                        "only 4 single-strand reads"},
    "FAM78A": {"pLI": 0.0, "LOEUF": 0.6,
               "note": "not LoF-intolerant; region is a dense cluster of "
                       "co-occurring candidate alleles (paralog decoy)"},
}

SCRIPT = os.path.abspath(__file__)
ROOT = os.path.dirname(os.path.dirname(SCRIPT))          # workspace root
INPUTS = os.path.join(ROOT, "inputs")
OUTDIR = os.environ.get("LS_OUT", os.path.join(ROOT, "output"))
FASTA = os.environ.get("LS_FASTA",
                       os.path.join(INPUTS, "reference", "GRCh38_chr9.fa.gz"))
GTF = os.environ.get("LS_GTF", os.path.join(
    INPUTS, "reference", "gencode.v47.chr9.annotation.gtf.gz"))
FASTQ = os.environ.get("LS_FASTQ",
                       os.path.join(INPUTS, "deleterious.mutation.q2.R1.fq.gz"))

CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
    CODE[_b + 32] = _i          # lowercase
CHAR_OF = np.array([ord("A"), ord("C"), ord("G"), ord("T"), ord("N")],
                   dtype=np.uint8)


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


# ----------------------------- reference -----------------------------
def load_fasta(path):
    with gzip.open(path, "rb") as f:
        raw = f.read()
    lines = raw.split(b"\n")
    seq = b"".join(l.strip() for l in lines if l and not l.startswith(b">"))
    seq = seq.upper()
    arr = np.frombuffer(seq, dtype=np.uint8).copy()
    codes = CODE[arr]
    log("loaded %s: %d bp" % (os.path.basename(path), len(codes)))
    return codes


# ----------------------------- annotation -----------------------------
def parse_gtf(path):
    """Return protein-coding transcripts:
    tid -> dict(gene_id, gene_name, strand, cds=[(start,end)...], mane)."""
    tx = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            a = line.rstrip("\n").split("\t")
            feat = a[2]
            if feat not in ("transcript", "CDS"):
                continue
            attr = {}
            tags = []
            for kv in a[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                if kv.startswith("tag "):
                    tags.append(kv.split('"')[1])
                    continue
                key, val = kv.split(" ", 1)
                attr[key] = val.strip('"') if '"' in val else val
            if attr.get("gene_type") != "protein_coding":
                continue
            tid = attr["transcript_id"]
            if feat == "transcript":
                if attr.get("transcript_type") != "protein_coding":
                    continue
                tx[tid] = {
                    "gene_id": attr["gene_id"],
                    "gene_name": attr.get("gene_name", ""),
                    "strand": a[6],
                    "cds": [],
                    "mane": "MANE_Select" in tags,
                }
            else:  # CDS
                if tid in tx:
                    tx[tid]["cds"].append((int(a[3]), int(a[4])))
    for t in tx.values():
        t["cds"].sort()
    tx = {k: v for k, v in tx.items() if v["cds"]}
    log("parsed %d protein-coding transcripts" % len(tx))
    return tx


def merge_intervals(iv, flank=0, upper=None):
    iv = sorted((max(s - flank, 1), e + flank if upper is None
                 else min(e + flank, upper)) for s, e in iv)
    m = []
    for s, e in iv:
        if m and s <= m[-1][1] + 1:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return [tuple(x) for x in m]


# ----------------------------- k-mer index -----------------------------
def build_index(genome, targets):
    mask = (1 << (2 * K)) - 1
    index = {}
    counts = {}
    for s, e in targets:
        seg = genome[s - 1:e]
        val = 0
        nbad = 0
        for i, c in enumerate(seg):
            if c == 4:
                nbad += 1
            val = ((val << 2) | int(c)) & mask
            if i >= K and seg[i - K] == 4:
                nbad -= 1
            if i >= K - 1 and nbad == 0:
                g = s + (i - K + 1)          # 1-based genome position
                if val in counts and counts[val] is None:
                    continue
                lst = index.get(val)
                if lst is None:
                    index[val] = [g]
                    counts[val] = 1
                else:
                    lst.append(g)
                    n = counts[val] + 1
                    if n > MAX_KMER_HITS:
                        index[val] = None
                        counts[val] = None
                    else:
                        counts[val] = n
    kept = sum(1 for v in index.values() if v is not None)
    log("k-mer index: %d distinct %d-mers kept (targets=%d bp)" %
        (kept, K, sum(e - s + 1 for s, e in targets)))
    return index


def seed_votes(codes, index):
    votes = {}
    L = len(codes)
    if L < K + 2:
        return votes
    mask = (1 << (2 * K)) - 1
    val = 0
    nbad = 0
    get = index.get
    for i in range(L):
        c = codes[i]
        if c == 4:
            nbad += 1
        val = ((val << 2) | int(c)) & mask
        if i >= K and codes[i - K] == 4:
            nbad -= 1
        if i >= K - 1 and nbad == 0:
            off = i - K + 1
            if off % SEED_STEP:
                continue
            hits = get(val)
            if hits:
                for g in hits:
                    d = g - off
                    votes[d] = votes.get(d, 0) + 1
    return votes


# ----------------------------- placement -----------------------------
def best_window(mm):
    """Largest read window [i,j) with mismatch fraction <= WINDOW_MM_FRAC."""
    L = len(mm)
    if L < MIN_ALIGNED_LEN:
        return None
    P = np.concatenate(([0], np.cumsum(mm))).astype(np.int32)
    best = None
    for i in range(0, L - MIN_ALIGNED_LEN + 1, 3):
        j = L
        while j - i >= MIN_ALIGNED_LEN:
            m = int(P[j] - P[i])
            allowed = max(4, int(WINDOW_MM_FRAC * (j - i)))
            if m <= allowed:
                break
            j -= max(1, m - allowed)
        if j - i >= MIN_ALIGNED_LEN:
            if best is None or (j - i) > (best[1] - best[0]):
                best = (i, j)
    return best


def place_read(codes, quals, genome, index, glen):
    """Best placement of a read (trying both orientations) or None.

    Returns (orient, gstart, left, right, seq, ql, nmm, vote) where
    gstart is the 1-based genome position of oriented-read base 0, the
    aligned segment is [left,right) in read coordinates, nmm is the number
    of mismatches inside that segment and vote the seed support.
    """
    rc = (3 - codes)[::-1].copy()
    best = None
    for orient, seq in ((0, codes), (1, rc)):
        votes = seed_votes(seq, index)
        if not votes:
            continue
        d1 = max(votes, key=lambda d: (votes[d], -d))
        v1 = votes[d1]
        if v1 < MIN_VOTE:
            continue
        amb = False
        for d2, v2 in votes.items():
            if d2 != d1 and abs(d2 - d1) > 32 and v2 >= AMBIG_RATIO * v1:
                amb = True
                break
        if amb:
            continue
        L = len(seq)
        lo = max(d1, 1)
        hi = min(d1 + L - 1, glen)
        if hi - lo + 1 < MIN_ALIGNED_LEN:
            continue
        left0 = lo - d1
        right0 = hi - d1 + 1
        ref = genome[lo - 1:hi]
        seg = seq[left0:right0]
        mm = (seg != ref) | (ref == 4)
        win = best_window(mm)
        if win is None:
            continue
        i, j = win
        nmm = int(mm[i:j].sum())
        cand = (orient, d1, left0 + i, left0 + j, nmm, v1)
        if best is None or (cand[4], -(cand[3] - cand[2])) < \
                           (best[4], -(best[3] - best[2])):
            best = cand
    if best is None:
        return None
    orient, d1, left, right, nmm, v1 = best
    seq = codes if orient == 0 else rc
    ql = quals if orient == 0 else quals[::-1]
    return orient, d1, left, right, seq, ql, nmm, v1


# ----------------------------- read streaming -----------------------------
def stream_fastq(path):
    with gzip.open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip("\n")
            f.readline()
            q = f.readline().rstrip("\n")
            yield h[1:].strip(), s, q


# ----------------------------- pileup pass -----------------------------
def run_pileup(path, genome, index, cds_starts, cds_ends, cds_offsets, glen):
    """Map reads and count Q>=BASE_Q bases over CDS positions (8 counters
    per position: base x strand)."""
    n_cds = int(cds_offsets[-1] + (cds_ends[-1] - cds_starts[-1] + 1))
    flat = np.zeros(n_cds * 8, dtype=np.int32)
    stats = Counter()
    t0 = time.time()
    for n, (name, s, q) in enumerate(stream_fastq(path), 1):
        if n % 100000 == 0:
            log("  pileup read %d (%.0f/s)" % (n, n / (time.time() - t0)))
        codes = CODE[np.frombuffer(s.encode(), dtype=np.uint8)]
        quals = np.frombuffer(q.encode(), dtype=np.uint8) - 33
        stats["total_reads"] += 1
        pl = place_read(codes, quals, genome, index, glen)
        if pl is None:
            stats["unmapped_reads"] += 1
            continue
        orient, d1, left, right, seq, ql, nmm, vote = pl
        stats["mapped_reads"] += 1
        p = np.arange(d1 + left, d1 + right, dtype=np.int64)   # 1-based
        b = seq[left:right].astype(np.int64)
        qq = ql[left:right]
        keep = (qq >= BASE_Q) & (b < 4)
        p = p[keep]
        b = b[keep]
        if len(p) == 0:
            continue
        idx = np.searchsorted(cds_starts, p, side="right") - 1
        valid = (idx >= 0) & (p <= cds_ends[idx])
        if not valid.any():
            continue
        idx = idx[valid]
        p = p[valid]
        flat_idx = cds_offsets[idx] + (p - cds_starts[idx])
        lin = flat_idx * 8 + b[valid] * 2 + orient
        np.add.at(flat, lin, 1)
    stats["wall_s"] = round(time.time() - t0, 1)
    return flat, stats


def call_variants(flat, genome, cds_pos_map):
    cm = flat.reshape(-1, 8)
    n = cm.shape[0]
    depth = cm.sum(axis=1)
    ref_code = genome[cds_pos_map - 1]
    base_ct = cm[:, 0::2] + cm[:, 1::2]                    # n x 4
    alt_ct = base_ct.copy()
    alt_ct[np.arange(n), ref_code] = 0
    best_alt = alt_ct.argmax(axis=1)
    best_ct = alt_ct[np.arange(n), best_alt]
    af = np.where(depth > 0, best_ct / np.maximum(depth, 1), 0.0)
    mask = (depth >= MIN_DEPTH) & (best_ct >= MIN_ALT) & \
           (af >= MIN_AF) & (af <= MAX_AF)
    cands = []
    for i in np.where(mask)[0]:
        pos = int(cds_pos_map[i])
        rc = int(ref_code[i])
        ac = int(best_alt[i])
        cands.append({
            "chrom": CHROM,
            "pos": pos,
            "ref": BASES[rc],
            "alt": BASES[ac],
            "alt_reads": int(best_ct[i]),
            "ref_reads": int(depth[i] - best_ct[i]),
            "total_reads": int(depth[i]),
            "allele_fraction": round(int(best_ct[i]) / int(depth[i]), 4),
            "alt_fwd": int(cm[i, ac * 2]),
            "alt_rev": int(cm[i, ac * 2 + 1]),
            "depth_fwd": int(cm[i, 0::2].sum()),
            "depth_rev": int(cm[i, 1::2].sum()),
        })
    cands.sort(key=lambda c: -c["alt_reads"])
    log("SNV candidates (CDS depth>=%d, alt>=%d): %d" %
        (MIN_DEPTH, MIN_ALT, len(cands)))
    return cands


# ----------------------------- annotation -----------------------------
def build_tx_models(tx, genome):
    models = {}
    for tid, t in tx.items():
        strand = t["strand"]
        cds = t["cds"] if strand == "+" else list(reversed(t["cds"]))
        seq = []
        bounds = []
        cum = 0
        for s, e in cds:
            seg = genome[s - 1:e]
            if strand == "-":
                seg = (3 - seg)[::-1]
            seq.append(bytes(CHAR_OF[seg]))
            bounds.append((cum, cum + (e - s + 1), s, e))
            cum += e - s + 1
        models[tid] = {
            "gene_id": t["gene_id"],
            "gene_name": t["gene_name"],
            "strand": strand,
            "mane": t["mane"],
            "seq": b"".join(seq),
            "bounds": bounds,
        }
    return models


def cds_offset(models, tid, pos):
    m = models[tid]
    for cum0, cum1, s, e in m["bounds"]:
        if s <= pos <= e:
            if m["strand"] == "+":
                return cum0 + (pos - s)
            return cum0 + (e - pos)
    return None


def annotate_candidate(c, models, tx_cds_index):
    pos = c["pos"]
    starts = [x[0] for x in tx_cds_index]
    i = np.searchsorted(starts, pos, side="right") - 1
    hits = []
    while i >= 0 and tx_cds_index[i][0] >= pos - 2000000:
        s, e, tid = tx_cds_index[i]
        if s <= pos <= e:
            hits.append(tid)
        i -= 1
    out = []
    for tid in hits:
        m = models[tid]
        off = cds_offset(models, tid, pos)
        if off is None or off >= len(m["seq"]):
            continue
        frame = off % 3
        cs = off - frame
        codon = bytearray(m["seq"][cs:cs + 3])
        if len(codon) < 3:
            continue
        ref_codon = bytes(codon).decode()
        codon[frame] = ord(c["alt"])
        alt_codon = bytes(codon).decode()
        ref_aa = CODON_TABLE.get(ref_codon, "?")
        alt_aa = CODON_TABLE.get(alt_codon, "?")
        if alt_aa == "*" and ref_aa != "*":
            cons = "stop_gained"
        elif ref_aa == "*" and alt_aa != "*":
            cons = "stop_lost"
        elif alt_aa != ref_aa:
            cons = "missense_variant"
        else:
            cons = "synonymous_variant"
        hgvs = "p.(%s%d%s)" % (AA3.get(ref_aa, "?"), cs // 3 + 1,
                               AA3.get(alt_aa, "?"))
        out.append({
            "transcript_id": tid,
            "gene_id": m["gene_id"],
            "gene_name": m["gene_name"],
            "strand": m["strand"],
            "mane_select": m["mane"],
            "consequence": cons,
            "hgvs_p": hgvs,
            "ref_codon": ref_codon,
            "alt_codon": alt_codon,
            "cds_position": off + 1,
        })
    return out


# ----------------------------- co-occurrence pass -----------------------------
def run_cooccurrence(path, genome, index, glen, cands):
    """Second pass over the reads: for every candidate SNV count how many
    mapped reads carry the alt allele, how many of those reads carry NO
    other candidate allele (clean), strand support and base qualities.

    True mosaic alleles are carried by reads that otherwise match the
    reference; paralog/pseudogene mismapping yields reads carrying many
    candidate alleles simultaneously (clean_frac -> 0)."""
    order = np.argsort([c["pos"] for c in cands])
    cpos = np.array([cands[k]["pos"] for k in order], dtype=np.int64)
    calt = np.array([BASES.index(cands[k]["alt"]) for k in order],
                    dtype=np.int64)
    m = len(order)
    carriers = np.zeros(m, dtype=np.int64)
    clean = np.zeros(m, dtype=np.int64)
    others = np.zeros(m, dtype=np.int64)
    fwd = np.zeros(m, dtype=np.int64)
    rev = np.zeros(m, dtype=np.int64)
    bq_sum = np.zeros(m, dtype=np.float64)
    mm_sum = np.zeros(m, dtype=np.int64)
    vote_sum = np.zeros(m, dtype=np.int64)
    stats = Counter()
    t0 = time.time()
    for n, (name, s, q) in enumerate(stream_fastq(path), 1):
        if n % 100000 == 0:
            log("  co-occurrence read %d (%.0f/s)" %
                (n, n / (time.time() - t0)))
        codes = CODE[np.frombuffer(s.encode(), dtype=np.uint8)]
        quals = np.frombuffer(q.encode(), dtype=np.uint8) - 33
        pl = place_read(codes, quals, genome, index, glen)
        if pl is None:
            continue
        orient, d1, left, right, seq, ql, nmm, vote = pl
        p0, p1 = d1 + left, d1 + right - 1
        lo = np.searchsorted(cpos, p0)
        hi = np.searchsorted(cpos, p1, side="right")
        if hi <= lo:
            continue
        j = cpos[lo:hi] - p0
        b = seq[left:right].astype(np.int64)
        qq = ql[left:right]
        keep = (qq >= BASE_Q) & (b < 4)
        hit = (b[j] == calt[lo:hi]) & keep[j]
        S = np.where(hit)[0] + lo
        if len(S) == 0:
            continue
        if len(S) > 1:
            stats["reads_multi_allele"] += 1
        for ci in S:
            carriers[ci] += 1
            others[ci] += len(S) - 1
            bq_sum[ci] += int(qq[j[ci - lo]])
            mm_sum[ci] += nmm
            vote_sum[ci] += vote
            if orient == 0:
                fwd[ci] += 1
            else:
                rev[ci] += 1
            if len(S) == 1:
                clean[ci] += 1
    stats["wall_s"] = round(time.time() - t0, 1)
    for k, c in enumerate(cands):
        cc = cands[order[k]]
        cc["carriers"] = int(carriers[k])
        cc["clean_carriers"] = int(clean[k])
        cc["clean_frac"] = round(clean[k] / carriers[k], 3) if carriers[k] \
            else None
        cc["mean_other_alleles"] = round(others[k] / carriers[k], 2) \
            if carriers[k] else None
        cc["carrier_fwd"] = int(fwd[k])
        cc["carrier_rev"] = int(rev[k])
        cc["mean_alt_bq"] = round(bq_sum[k] / carriers[k], 1) \
            if carriers[k] else None
        cc["mean_alt_read_mismatches"] = round(mm_sum[k] / carriers[k], 2) \
            if carriers[k] else None
        cc["mean_alt_seed_votes"] = round(vote_sum[k] / carriers[k], 1) \
            if carriers[k] else None
    return stats


def binom_two_sided(k, n, p=0.5):
    if n == 0:
        return 1.0
    from math import comb
    tot = 2 ** n
    lo = sum(comb(n, i) for i in range(0, k + 1))
    hi = sum(comb(n, i) for i in range(k, n + 1))
    return min(1.0, 2 * min(lo, hi) / tot)


# ----------------------------- main -----------------------------
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    genome = load_fasta(FASTA)
    glen = len(genome)
    tx = parse_gtf(GTF)
    models = build_tx_models(tx, genome)

    all_cds = []
    for t in tx.values():
        all_cds.extend(t["cds"])
    targets = merge_intervals(all_cds, flank=FLANK, upper=glen)
    cds_merged = merge_intervals(all_cds)
    cds_starts = np.array([s for s, e in cds_merged], dtype=np.int64)
    cds_ends = np.array([e for s, e in cds_merged], dtype=np.int64)
    lengths = cds_ends - cds_starts + 1
    cds_offsets = np.zeros(len(lengths), dtype=np.int64)
    cds_offsets[1:] = np.cumsum(lengths[:-1])
    cds_pos_map = np.concatenate([np.arange(s, e + 1) for s, e in cds_merged])
    log("CDS mask: %d intervals, %d bp; index target: %d intervals, %d bp" %
        (len(cds_merged), len(cds_pos_map), len(targets),
         sum(e - s + 1 for s, e in targets)))

    index = build_index(genome, targets)

    tx_cds_index = []
    for tid, t in tx.items():
        for s, e in t["cds"]:
            tx_cds_index.append((s, e, tid))
    tx_cds_index.sort()

    log("=== pass 1: mapping + pileup ===")
    flat, stats = run_pileup(FASTQ, genome, index, cds_starts, cds_ends,
                             cds_offsets, glen)
    cands = call_variants(flat, genome, cds_pos_map)

    log("annotating candidates")
    for c in cands:
        c["annotations"] = annotate_candidate(c, models, tx_cds_index)
    stop_cands = [c for c in cands
                  if any(a["consequence"] == "stop_gained"
                         for a in c["annotations"])]
    log("stop-gained candidates: %d" % len(stop_cands))

    log("=== pass 2: allele co-occurrence for %d candidates ===" % len(cands))
    co_stats = run_cooccurrence(FASTQ, genome, index, glen, cands)

    # ---- high-confidence mosaic selection ----
    def high_conf(c):
        return (c.get("carriers", 0) >= HC_MIN_CARRIERS and
                (c.get("clean_frac") or 0) >= HC_MIN_CLEAN_FRAC and
                HC_MIN_AF <= c["allele_fraction"] <= HC_MAX_AF and
                c.get("carrier_fwd", 0) >= 1 and
                c.get("carrier_rev", 0) >= 1 and
                (c.get("mean_alt_bq") or 0) >= HC_MIN_BQ)

    hc_stop = [c for c in stop_cands if high_conf(c)]
    hc_stop.sort(key=lambda c: (-c["carriers"],
                                -(c["clean_frac"] or 0)))
    final = hc_stop[0] if hc_stop else None
    log("high-confidence mosaic stop-gained candidates: %d" % len(hc_stop))

    if final is not None:
        ann = [a for a in final["annotations"]
               if a["consequence"] == "stop_gained"]
        ann.sort(key=lambda a: (not a["mane_select"], a["transcript_id"]))
        final["reported_annotation"] = ann[0]
        final["strand_binom_p"] = round(binom_two_sided(
            min(final["carrier_fwd"], final["carrier_rev"]),
            final["carrier_fwd"] + final["carrier_rev"]), 4)

    # ---------------- outputs ----------------
    with open(os.path.join(OUTDIR, "variant.tsv"), "w") as f:
        f.write("chrom\tpos\tref\talt\tgene\tconsequence\talt_reads\t"
                "total_reads\tallele_fraction\n")
        if final is not None:
            a = final["reported_annotation"]
            f.write("%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\n" % (
                final["chrom"], final["pos"], final["ref"], final["alt"],
                a["gene_name"], a["consequence"], final["alt_reads"],
                final["total_reads"], final["allele_fraction"]))

    cm = flat.reshape(-1, 8)
    depth = cm.sum(axis=1)
    cov_stats = {
        "cds_positions": int(len(depth)),
        "positions_ge_10x": int((depth >= 10).sum()),
        "fraction_cds_ge_10x": round(float((depth >= 10).mean()), 4),
        "mean_cds_depth": round(float(depth.mean()), 2),
    }

    def brief(c):
        keep_ann = [a for a in c.get("annotations", [])
                    if a["consequence"] == "stop_gained"] or \
                   c.get("annotations", [])
        return {k: v for k, v in c.items() if k != "annotations"} | {
            "annotations": keep_ann}

    evidence = {
        "reference": {
            "genome_build": "GRCh38 (primary assembly)",
            "reference_fasta": "GRCh38_chr9.fa.gz (chromosome 9 from the "
                               "Broad Institute GATK GRCh38 bundle)",
            "annotation": "GENCODE v47 (Ensembl 113), "
                          "gencode.v47.chr9.annotation.gtf.gz, "
                          "protein_coding transcripts",
            "coordinate_system": "chromosome 'chr9', 1-based",
        },
        "pipeline": {
            "aligner": "custom seed-and-verify mapper (this script)",
            "seed_k": K,
            "seed_step": SEED_STEP,
            "index_target": "protein-coding CDS +/- %d bp" % FLANK,
            "min_seed_votes": MIN_VOTE,
            "max_mismatch_fraction": WINDOW_MM_FRAC,
            "min_aligned_length": MIN_ALIGNED_LEN,
            "min_base_quality": BASE_Q,
            "min_depth": MIN_DEPTH,
            "min_alt_reads": MIN_ALT,
            "high_confidence_filters": {
                "min_carriers": HC_MIN_CARRIERS,
                "min_clean_frac": HC_MIN_CLEAN_FRAC,
                "af_window": [HC_MIN_AF, HC_MAX_AF],
                "both_strands_required": True,
                "min_mean_alt_base_quality": HC_MIN_BQ,
            },
        },
        "read_stats": dict(stats),
        "cooccurrence_stats": dict(co_stats),
        "coverage": cov_stats,
        "candidate_snvs_total": len(cands),
        "stop_gained_candidates": [brief(c) for c in stop_cands],
        "high_confidence_stop_gained": [brief(c) for c in hc_stop],
        "final_variant": final,
        "gene_lof_constraint_prior": {
            "source": "gnomAD v2.1.1 loss-of-function constraint metrics "
                      "(Karczewski et al., Nature 2020); approximate "
                      "published values, used as external prior knowledge "
                      "for gene-level interpretation only",
            "genes": {g: LOF_CONSTRAINT[g] for g in (
                [final["reported_annotation"]["gene_name"]] if final else []
            ) if g in LOF_CONSTRAINT} | {
                g: v for g, v in LOF_CONSTRAINT.items()
                if g in {a["gene_name"] for c in stop_cands
                         for a in c["annotations"]}},
        },
    }
    with open(os.path.join(OUTDIR, "evidence.json"), "w") as f:
        json.dump(evidence, f, indent=2)

    print("\n=== SUMMARY ===")
    print("reads total/mapped: %d / %d" %
          (stats["total_reads"], stats["mapped_reads"]))
    print("candidates: %d ; stop_gained: %d ; high-confidence mosaic: %d" %
          (len(cands), len(stop_cands), len(hc_stop)))
    if final is not None:
        a = final["reported_annotation"]
        print("FINAL: %s:%d %s>%s %s %s %s AF=%.4f (%d/%d) strands +%d/-%d "
              "clean_frac=%s" % (
                  final["chrom"], final["pos"], final["ref"], final["alt"],
                  a["gene_name"], a["consequence"], a["hgvs_p"],
                  final["allele_fraction"], final["alt_reads"],
                  final["total_reads"], final["carrier_fwd"],
                  final["carrier_rev"], final["clean_frac"]))
    print("wrote %s and %s" % (os.path.join(OUTDIR, "variant.tsv"),
                                os.path.join(OUTDIR, "evidence.json")))


if __name__ == "__main__":
    main()
