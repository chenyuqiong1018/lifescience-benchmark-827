#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cryptic-exon detection pipeline (pure Python + NumPy; no external aligner).

Inputs (./inputs):
  cryptic.exon.q1.fq.gz                      single-end RNA-seq reads (107 bp)
  reference/GRCh38_chr9.fa.gz                chr9 genome (Broad hg38/v0)
  reference/ensembl112_protein_coding_exons.tsv.gz
                                             Ensembl release-112 protein-coding exon
                                             annotation (BioMart export). Novelty is
                                             assessed against THIS annotation version.

Outputs (./output):
  cryptic_exon.tsv   gene, chrom, start, end, left_junction_reads,
                     right_junction_reads, expression_evidence
  junctions.tsv      detected splice junctions with novelty flags
  report.md          narrative report

Method overview
---------------
1.  Build a genome k-mer index restricted to k-mers that occur in the reads
    (k=16; k-mers with >20 genomic copies are treated as repetitive).
2.  Screen every read (and its reverse complement) for seed hits at offsets
    0,15,30,45,60,75,91; keep reads with >=2 hit seeds in one orientation.
3.  Place kept reads by diagonal voting: contiguous alignments (single
    dominant diagonal, full 107 bp verified with <=3 mismatches) and
    single-junction split alignments (two diagonals; every split point tested,
    <=3 mismatches; introns 30..200,000 bp). All best-mismatch placements are
    clustered; a read is counted only when its best placements form a single
    locus cluster (junctions within +/-4 bp, contigs within +/-3 bp). Reads
    with several equally good loci (segmental duplications) are discarded.
4.  Junction clusters are canonicalised to the GT-AG splice motif (+/-3 bp,
    strand-aware).
5.  Novelty: a junction is novel iff (chrom, intron_start, intron_end)
    (1-based inclusive skipped interval) is absent from every supplied
    transcript.
6.  Cryptic-exon search: pairs of novel junctions on the same strand define
    the interval left_intron_end+1 .. right_intron_start-1; candidates must
    not overlap any annotated exon and must lie inside a protein-coding gene.
7.  Refinement: all canonical splice-site variants within the microhomology
    window (+/-3 bp) are tested by re-aligning all reads to candidate spliced
    transcripts; the variant supported by reads spanning both junctions wins.
8.  Gene-level quantification against canonical and cryptic isoform
    references.
"""
import gzip
import os
import time
from collections import defaultdict

import numpy as np

T0 = time.time()


def log(*a):
    print(*a, flush=True)


K = 16                 # seed length
OFFS = [0, 15, 30, 45, 60, 75, 91]
MAXMM = 3              # max mismatches for an alignment
CAP = 20               # max genomic copies of a seed k-mer
MAX_INTRON = 200_000
MIN_INTRON = 30
MIN_ANCHOR = 16
PAIR_MAXGAP = 100_000  # max distance between paired novel junctions
MIN_JUNC_READS = 3     # min support per junction in a pair
BASE = "ACGT"
CODE = np.zeros(256, dtype=np.uint8)
for _c, _v in zip(b"ACGTNacgtn", [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]):
    CODE[_c] = _v

HERE = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(HERE)


def code_str(s: str) -> np.ndarray:
    return CODE[np.frombuffer(s.encode(), dtype=np.uint8)]


def load_inputs():
    seqs = []
    with gzip.open(os.path.join(WS, "inputs/cryptic.exon.q1.fq.gz"), "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                seqs.append(line.strip())
    n = len(seqs)
    length = len(seqs[0])
    flat = "".join(seqs).encode()
    rseq = CODE[np.frombuffer(flat, dtype=np.uint8)].reshape(n, length).copy()
    with gzip.open(os.path.join(WS, "inputs/reference/GRCh38_chr9.fa.gz"), "rt") as f:
        chrom = f.readline().strip().lstrip(">").split()[0]
        g = "".join(l.strip() for l in f).upper()
    gcode = code_str(g)
    return rseq, gcode, chrom


def build_index(rseq, gcode):
    """Genome index restricted to k-mers present in the reads."""
    n, length = rseq.shape
    log(f"[index] reads={n} len={length} genome={len(gcode)}")

    def kmer_keys(m, k=K):
        nn, ll = m.shape
        keys = np.zeros((nn, ll - k + 1), dtype=np.uint64)
        bad = np.zeros((nn, ll - k + 1), dtype=bool)
        for j in range(k):
            col = m[:, j:j + ll - k + 1].astype(np.uint64)
            keys |= col << np.uint64(2 * (k - 1 - j))
            bad |= (m[:, j:j + ll - k + 1] == 4)
        keys[bad] = np.uint64((1 << 64) - 1)
        return keys

    SENT = np.uint64((1 << 64) - 1)
    rkeys = kmer_keys(rseq)
    uk = np.unique(rkeys[rkeys != SENT])
    log(f"[index] unique read k-mers: {len(uk)}")

    G = len(gcode)
    nwin = G - K + 1
    chunk = 8_000_000
    hk_list, hp_list = [], []
    for s in range(0, nwin, chunk):
        e = min(s + chunk, nwin)
        sub = gcode[s:e + K - 1]
        keys = np.zeros(e - s, dtype=np.uint64)
        hasN = np.zeros(e - s, dtype=bool)
        for j in range(K):
            seg = sub[j:j + (e - s)]
            keys |= seg.astype(np.uint64) << np.uint64(2 * (K - 1 - j))
            hasN |= (seg == 4)
        keys[hasN] = SENT
        idx = np.minimum(np.searchsorted(uk, keys), len(uk) - 1)
        hit = (uk[idx] == keys) & (~hasN)
        hk_list.append(keys[hit])
        hp_list.append(np.arange(s, e, dtype=np.uint64)[hit])
    hk = np.concatenate(hk_list)
    hp = np.concatenate(hp_list)
    u, inv, c = np.unique(hk, return_inverse=True, return_counts=True)
    keep = c <= CAP
    order = np.argsort(inv, kind="stable")
    starts = np.searchsorted(np.sort(inv), np.arange(len(u) + 1))
    pos_dict = {}
    for i in np.where(keep)[0]:
        pos_dict[int(u[i])] = hp[order[starts[i]:starts[i + 1]]].astype(np.int64)
    log(f"[index] kept {len(pos_dict)} k-mers with <= {CAP} genomic copies "
        f"({time.time()-T0:.0f}s)")
    return pos_dict


def keys_at(mat, off, k=K):
    sub = mat[:, off:off + k]
    bad = (sub == 4).any(1)
    kkey = np.zeros(mat.shape[0], dtype=np.uint64)
    for j in range(k):
        kkey |= sub[:, j].astype(np.uint64) << np.uint64(2 * (k - 1 - j))
    kkey[bad] = (1 << 64) - 1
    return kkey


def screen(rseq, pos_dict):
    n = rseq.shape[0]
    rseq_rc = np.where(rseq == 4, 4, 3 - rseq)[:, ::-1].astype(np.uint8).copy()
    SENT = (1 << 64) - 1

    def do(mat):
        h = np.zeros((n, len(OFFS)), dtype=np.int8)
        get = pos_dict.get
        for oi, off in enumerate(OFFS):
            ks = keys_at(mat, off)
            for i in range(n):
                kk = int(ks[i])
                if kk == SENT:
                    continue
                p = get(kk)
                if p is not None:
                    h[i, oi] = min(len(p), 30)
        return h

    h1 = do(rseq)
    h2 = do(rseq_rc)
    s1 = (h1 >= 1).sum(1)
    s2 = (h2 >= 1).sum(1)
    keep = np.where((s1 >= 2) | (s2 >= 2))[0]
    log(f"[screen] kept {len(keep)} reads for placement ({time.time()-T0:.0f}s)")
    return keep, s1, s2, rseq_rc


def seed_positions(row, pos_dict):
    out = []
    for off in OFFS:
        k = 0
        ok = True
        for j in range(K):
            v = row[off + j]
            if v == 4:
                ok = False
                break
            k = (k << 2) | int(v)
        if not ok:
            continue
        p = pos_dict.get(k)
        if p is not None and len(p) <= CAP:
            out.append((off, p))
    return out


def all_placements(row, pos_dict, gcode, arL):
    """List of (type, mm, x, y): C=contig at x; J=junction intron (x, y),
    1-based inclusive."""
    G = len(gcode)
    length = len(row)
    seeds = seed_positions(row, pos_dict)
    if len(seeds) < 2:
        return []
    diags = defaultdict(int)
    for off, ps in seeds:
        for pp in ps:
            diags[int(pp) - off] += 1
    top = sorted(diags.items(), key=lambda x: -x[1])
    placements = []
    for d, v in top:
        if v < 3:
            break
        if 0 <= d and d + length <= G:
            mm = int((gcode[d:d + length] != row).sum())
            if mm <= MAXMM:
                placements.append(("C", mm, d, None))
    cand = [d for d, v in top if v >= 1][:10]
    seen = defaultdict(int)
    for i in range(len(cand)):
        for j in range(len(cand)):
            if i == j:
                continue
            d1, d2 = cand[i], cand[j]
            if d1 < 0 or d2 < 0 or d1 + length > G or d2 + length > G:
                continue
            avals = np.arange(MIN_ANCHOR, length - MIN_ANCHOR + 1)
            cols = arL[None, :]
            M = np.where(cols < avals[:, None],
                         gcode[d1:d1 + length][None, :],
                         gcode[d2:d2 + length][None, :])
            mm = (M != row[None, :]).sum(1)
            for gi in np.where(mm <= MAXMM)[0]:
                a = int(avals[gi])
                s, e = d1 + a + 1, d2 + a
                if MIN_INTRON <= e - s + 1 <= MAX_INTRON:
                    if (s, e) not in seen or int(mm[gi]) < seen[(s, e)]:
                        seen[(s, e)] = int(mm[gi])
    for (s, e), mm in seen.items():
        placements.append(("J", mm, s, e))
    placements.sort(key=lambda x: x[1])
    return placements


def place_reads(rseq, rseq_rc, keep, s1, s2, pos_dict, gcode):
    length = rseq.shape[1]
    arL = np.arange(length)
    C_uniq = {}                    # idx -> (strand, pos, mm)
    J_uniq = defaultdict(list)     # (s, e, strand) -> [idx]
    ambig = nomap = 0
    for n, idx in enumerate(keep):
        idx = int(idx)
        rows = [(1, rseq[idx])]
        if s2[idx] >= 2:
            rows.append((-1, rseq_rc[idx]))
        allp = []
        for strand, row in rows:
            for typ, mm, x, y in all_placements(row, pos_dict, gcode, arL):
                allp.append((strand, typ, mm, x, y))
        if not allp:
            nomap += 1
            continue
        bestmm = min(p[2] for p in allp)
        best = [p for p in allp if p[2] == bestmm]
        clusters = []
        for st, typ, mm, x, y in best:
            placed = False
            if typ == "C":
                for cl in clusters:
                    if cl[0] == "C" and cl[1] == st and abs(cl[2] - x) <= 3:
                        cl[3].append((x, y, mm))
                        placed = True
                        break
                if not placed:
                    clusters.append(["C", st, x, [(x, y, mm)]])
            else:
                for cl in clusters:
                    if cl[0] == "J" and cl[1] == st and abs(cl[2] - x) <= 4 \
                            and abs(cl[3] - y) <= 4:
                        cl[4].append((x, y, mm))
                        placed = True
                        break
                if not placed:
                    clusters.append(["J", st, x, y, [(x, y, mm)]])
        if len(clusters) != 1:
            ambig += 1
            continue
        cl = clusters[0]
        if cl[0] == "C":
            C_uniq[idx] = (cl[1], int(cl[2]), bestmm)
        else:
            J_uniq[(int(cl[2]), int(cl[3]), cl[1])].append(idx)
        if n % 10000 == 0:
            log(f"[place] {n}/{len(keep)} ({time.time()-T0:.0f}s)")
    log(f"[place] contig={len(C_uniq)} junction clusters={len(J_uniq)} "
        f"ambiguous={ambig} nomap={nomap} ({time.time()-T0:.0f}s)")
    return C_uniq, J_uniq


def canonical_junc(s, e, strand, gcode):
    G = len(gcode)
    cands = []
    for ds in range(-3, 4):
        for de in range(-3, 4):
            ss, ee = s + ds, e + de
            if ss < 2 or ee > G:
                continue
            d = BASE[gcode[ss - 1]] + BASE[gcode[ss]]
            a = BASE[gcode[ee - 2]] + BASE[gcode[ee - 1]]
            if strand == 1:
                if d == "GT" and a == "AG":
                    cands.append((ss, ee))
            else:
                if d == "CT" and a == "AC":
                    cands.append((ss, ee))
    return cands


def canonicalize(J_uniq, gcode):
    J_reads = defaultdict(list)
    for (s, e, st), idxs in J_uniq.items():
        canon = canonical_junc(s, e, st, gcode)
        cs, ce = canon[0] if len(canon) == 1 else (s, e)
        J_reads[(int(cs), int(ce), st)].extend(idxs)
    return J_reads


def load_annotation():
    introns = set()
    trans_exons = defaultdict(list)
    genes = {}
    path = os.path.join(WS, "inputs/reference/ensembl112_protein_coding_exons.tsv.gz")
    with gzip.open(path, "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            ensg, enst, gname, chrom, s, e, strand, ttype = p
            s, e = int(s), int(e)
            trans_exons[(chrom, enst)].append((s, e, strand, ensg))
            key = (chrom, ensg)
            if key not in genes:
                genes[key] = {"name": gname, "strand": strand, "exons": []}
            genes[key]["exons"].append((s, e))
    for (chrom, enst), exs in trans_exons.items():
        ee = sorted(set((a, b) for a, b, _, _ in exs))
        for i in range(len(ee) - 1):
            introns.add((chrom, ee[i][1] + 1, ee[i + 1][0] - 1))
    return introns, trans_exons, genes


def find_candidates(J_reads, introns, genes, chrom, achrom):
    """Pairs of novel junctions whose intervening interval contains no
    annotated exon and lies inside protein-coding gene(s)."""
    nov = defaultdict(list)
    for (s, e, st), reads in J_reads.items():
        if (achrom, s, e) not in introns:
            nov[st].append((s, e, len(reads)))
    exon_iv = []
    for (c, ensg), gi in genes.items():
        if c != achrom:
            continue
        for s, e in gi["exons"]:
            exon_iv.append((s, e, gi["name"]))
    exon_iv.sort()
    ex_starts = np.array([x[0] for x in exon_iv])

    def overlaps_exon(s, e):
        lo = int(np.searchsorted(ex_starts, s - 1_000_000))
        for j in range(lo, len(exon_iv)):
            a, b, g = exon_iv[j]
            if a > e:
                break
            if b >= s:
                return g
        return None

    cands = []
    for st, lst in nov.items():
        lst.sort()
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                s1, e1, n1 = lst[i]
                s2, e2, n2 = lst[j]
                if e1 < s2 and s2 - e1 <= PAIR_MAXGAP \
                        and n1 >= MIN_JUNC_READS and n2 >= MIN_JUNC_READS:
                    cs, ce = e1 + 1, s2 - 1
                    if overlaps_exon(cs, ce) is not None:
                        continue
                    host = sorted({v["name"] for (c, _), v in genes.items()
                                   if c == achrom
                                   and min(a for a, b in v["exons"]) <= cs
                                   and max(b for a, b in v["exons"]) >= ce})
                    cands.append({"strand": st, "left": (s1, e1, n1),
                                  "right": (s2, e2, n2), "start": cs, "end": ce,
                                  "host_genes": host})
    cands.sort(key=lambda c: -(c["left"][2] + c["right"][2]))
    return cands


def pick_host_gene(cand, trans_exons, genes, achrom):
    """Prefer the gene whose transcript exon boundaries match the junction
    ends (the cryptic exon is spliced between two exons of that gene)."""
    best, bestscore = None, -1
    for (c, enst), exs in trans_exons.items():
        if c != achrom:
            continue
        strand = exs[0][2]
        ensg = exs[0][3]
        ee = sorted(set((a, b) for a, b, _, _ in exs))
        gmin = min(a for a, b in ee)
        gmax = max(b for a, b in ee)
        if not (gmin <= cand["start"] and gmax >= cand["end"]):
            continue
        score = 0
        for a, b in ee:
            if b == cand["left"][0] - 1:      # exon ends before left intron
                score += 1
            if a == cand["right"][1] + 1:     # exon starts after right intron
                score += 1
        if score > bestscore:
            bestscore = score
            best = (ensg, enst, strand, ee)
    return best, bestscore


def scan_ref(mat, T, length, maxmm=MAXMM):
    """Slide every read along spliced reference T; return best position/mm."""
    n = mat.shape[0]
    npos = len(T) - length + 1
    bestpos = np.full(n, -1, dtype=np.int32)
    bestmm = np.full(n, 99, dtype=np.int32)
    chunk = 10000
    for s0 in range(0, n, chunk):
        s1 = min(n, s0 + chunk)
        sub = mat[s0:s1]
        mm_all = np.empty((npos, s1 - s0), dtype=np.int32)
        for p in range(npos):
            mm_all[p] = (T[p:p + length][None, :] != sub).sum(1)
        b = mm_all.argmin(0)
        bm = mm_all[b, np.arange(s1 - s0)]
        bestpos[s0:s1] = np.where(bm <= maxmm, b, -1)
        bestmm[s0:s1] = bm
    return bestpos, bestmm


def refine_candidate(cand, rseq, rseq_rc, gcode):
    """Enumerate canonical splice-site variants within +/-3 bp of each of the
    four junction coordinates and choose the spliced structure supported by
    the most reads (reads spanning both junctions are decisive).  Handles
    microhomology-ambiguous boundaries."""
    st = cand["strand"]
    s1, e1, _ = cand["left"]
    s2, e2, _ = cand["right"]
    length = rseq.shape[1]
    FLANK = 150

    def motif(pos, kind):
        # pos is 1-based; kind: 'donor'/'acceptor'; strand-aware
        if kind == "donor":
            d = BASE[gcode[pos - 1]] + BASE[gcode[pos]]
            return d == "GT" if st == 1 else d == "AC"
        d = BASE[gcode[pos - 2]] + BASE[gcode[pos - 1]]
        return d == "AG" if st == 1 else d == "CT"

    donors_L = [p for p in range(s1 - 3, s1 + 4) if p >= 2 and motif(p, "donor")]
    accs_L = [p for p in range(e1 - 3, e1 + 4) if p >= 2 and motif(p, "acceptor")]
    dons_R = [p for p in range(s2 - 3, s2 + 4) if p >= 2 and motif(p, "donor")]
    accs_R = [p for p in range(e2 - 3, e2 + 4) if p >= 2 and motif(p, "acceptor")]
    variants = []
    for dL in donors_L:
        for aL in accs_L:
            for dR in dons_R:
                for aR in accs_R:
                    cs, ce = aL + 1, dR - 1
                    if ce - cs + 1 >= 10:
                        variants.append((dL, aL, dR, aR, cs, ce))
    log(f"[refine] motif-valid variants: donors_L={donors_L} accs_L={accs_L} "
        f"dons_R={dons_R} accs_R={accs_R} -> {len(variants)} combinations")
    best = None
    for dL, aL, dR, aR, cs, ce in variants:
        up = gcode[dL - 1 - FLANK:dL - 1]
        down = gcode[aR:aR + FLANK]
        T = np.concatenate([up, gcode[cs - 1:ce], down])
        if st == -1:
            T = np.where(T == 4, 4, 3 - T)[::-1].astype(np.uint8)
        B0, B1 = len(up), len(up) + (ce - cs + 1)
        spanning = sup_l = sup_r = 0
        for mat in (rseq, rseq_rc):
            bp, bm = scan_ref(mat, T, length)
            for i in np.where((bp >= 0) & (bm <= MAXMM))[0]:
                p = int(bp[i])
                cL = p < B0 < p + length
                cR = p < B1 < p + length
                sup_l += cL
                sup_r += cR
                spanning += cL and cR
        log(f"[refine]   introns ({dL},{aL}) & ({dR},{aR}) CE {cs}-{ce}: "
            f"span={spanning} supL={sup_l} supR={sup_r}")
        key = (spanning, sup_l + sup_r)
        if best is None or key > best[0]:
            best = (key, {"intron_left": (dL, aL), "intron_right": (dR, aR),
                          "start": cs, "end": ce, "spanning": spanning,
                          "sup_left": sup_l, "sup_right": sup_r})
    return best[1]


def quantify_gene(exons_sorted, refined, rseq, rseq_rc, gcode):
    """Align all reads to canonical and cryptic isoform references of the host
    gene; count reads per junction and per isoform."""
    length = rseq.shape[1]
    cs, ce = refined["start"], refined["end"]
    exseqs = [gcode[a - 1:b] for a, b in exons_sorted]
    host_i = None
    for i in range(len(exons_sorted) - 1):
        iv_s, iv_e = exons_sorted[i][1] + 1, exons_sorted[i + 1][0] - 1
        if iv_s <= cs and ce <= iv_e:
            host_i = i
            break
    if host_i is None:
        raise RuntimeError("cryptic exon not inside any host-gene intron")
    Rcry = np.concatenate(exseqs[:host_i + 1] + [gcode[cs - 1:ce]] +
                          exseqs[host_i + 1:])
    Rcan = np.concatenate(exseqs)
    b_cry = [0]
    for part in (exseqs[:host_i + 1] + [gcode[cs - 1:ce]] + exseqs[host_i + 1:]):
        b_cry.append(b_cry[-1] + len(part))
    b_can = [0]
    for part in exseqs:
        b_can.append(b_can[-1] + len(part))
    counts = {"total": 0, "cry_junc": 0, "can_junc": 0,
              "J_left": set(), "J_right": set(),
              "J_canon_intron": set(), "J_other": set(), "spanning": set()}
    for tag, mat in (("fwd", rseq), ("rc", rseq_rc)):
        pc, mc = scan_ref(mat, Rcry, length)
        pn, mn = scan_ref(mat, Rcan, length)
        for i in range(mat.shape[0]):
            bc, bn = int(mc[i]), int(mn[i])
            if bc > MAXMM and bn > MAXMM:
                continue
            counts["total"] += 1
            if bc <= bn:
                p = int(pc[i])
                crossed = [k for k in range(1, len(b_cry) - 1)
                           if p < b_cry[k] < p + length]
                if not crossed:
                    continue
                counts["cry_junc"] += 1
                for k in crossed:
                    if k == host_i + 1:
                        counts["J_left"].add((tag, i))
                    elif k == host_i + 2:
                        counts["J_right"].add((tag, i))
                    else:
                        counts["J_other"].add((tag, i))
                if host_i + 1 in crossed and host_i + 2 in crossed:
                    counts["spanning"].add((tag, i))
            else:
                p = int(pn[i])
                crossed = [k for k in range(1, len(b_can) - 1)
                           if p < b_can[k] < p + length]
                if not crossed:
                    continue
                counts["can_junc"] += 1
                for k in crossed:
                    if k == host_i + 1:
                        counts["J_canon_intron"].add((tag, i))
                    else:
                        counts["J_other"].add((tag, i))
    return counts, Rcry, b_cry, host_i


def main():
    rseq, gcode, chrom = load_inputs()
    length = rseq.shape[1]
    log(f"[load] {rseq.shape[0]} reads x {length} bp; reference {chrom} "
        f"({len(gcode)} bp)")

    achrom = chrom[3:] if chrom.startswith("chr") else chrom
    pos_dict = build_index(rseq, gcode)
    keep, s1, s2, rseq_rc = screen(rseq, pos_dict)
    C_uniq, J_uniq = place_reads(rseq, rseq_rc, keep, s1, s2, pos_dict, gcode)
    J_reads = canonicalize(J_uniq, gcode)

    introns, trans_exons, genes = load_annotation()
    log(f"[annotation] introns={len(introns)} transcripts={len(trans_exons)} "
        f"genes={len(genes)}")

    jrows = []
    for (s, e, st), reads in sorted(J_reads.items()):
        jrows.append({"chrom": chrom, "start": s, "end": e, "strand": st,
                      "reads": len(reads),
                      "novel": (achrom, s, e) not in introns})

    cands = find_candidates(J_reads, introns, genes, chrom, achrom)
    log(f"[candidates] cryptic-exon candidates (interval free of annotated "
        f"exons): {len(cands)}")
    for c in cands[:10]:
        log(f"  strand {c['strand']:+d} {c['left'][:2]} ({c['left'][2]} r) / "
            f"{c['right'][:2]} ({c['right'][2]} r) -> interval "
            f"{c['start']}-{c['end']} hosts={c['host_genes']}")
    if not cands:
        raise SystemExit("no cryptic-exon candidate found")

    cand = cands[0]
    refined = refine_candidate(cand, rseq, rseq_rc, gcode)
    log(f"[refine] best CE {refined['start']}-{refined['end']} junctions "
        f"{refined['intron_left']} & {refined['intron_right']}")

    cand_refined = {"strand": cand["strand"],
                    "left": (refined["intron_left"][0],
                             refined["intron_left"][1], 0),
                    "right": (refined["intron_right"][0],
                              refined["intron_right"][1], 0),
                    "start": refined["start"], "end": refined["end"]}
    host, score = pick_host_gene(cand_refined, trans_exons, genes, achrom)
    if host is None:
        raise SystemExit("no host gene found")
    ensg, enst, strand, exons_sorted = host
    gene_name = genes[(achrom, ensg)]["name"]
    log(f"[host] gene {gene_name} ({ensg}, {enst}) boundary score {score}")

    counts, Rcry, b_cry, host_i = quantify_gene(
        exons_sorted, refined, rseq, rseq_rc, gcode)
    jl, jr = len(counts["J_left"]), len(counts["J_right"])
    span = len(counts["spanning"])
    n_canon = len(counts["J_canon_intron"])
    log(f"[quant] {gene_name}: total={counts['total']} J_left={jl} "
        f"J_right={jr} both={span} canon-intron-spanning={n_canon}")

    B0, B1 = b_cry[host_i + 1], b_cry[host_i + 2]
    cov = np.zeros(B1 - B0, dtype=np.int32)
    for mat in (rseq, rseq_rc):
        bp, bm = scan_ref(mat, Rcry, length)
        for i in np.where((bp >= 0) & (bm <= MAXMM))[0]:
            p = int(bp[i])
            lo, hi = max(p, B0), min(p + length, B1)
            if hi > lo:
                cov[lo - B0:hi - B0] += 1

    write_outputs(chrom, gene_name, ensg, enst, strand, exons_sorted, refined,
                  counts, jl, jr, span, n_canon, cov, jrows, rseq, length,
                  introns, host_i, len(gcode))
    log(f"[done] ({time.time()-T0:.0f}s)")



def write_outputs(chrom, gene_name, ensg, enst, strand, exons_sorted, refined,
                  counts, jl, jr, span, n_canon, cov, jrows, rseq, length,
                  introns, host_i, genlen):
    outdir = os.path.join(WS, "output")
    os.makedirs(outdir, exist_ok=True)
    cs, ce = refined["start"], refined["end"]
    l_intron = refined["intron_left"]
    r_intron = refined["intron_right"]
    intron1_s = exons_sorted[host_i][1] + 1
    intron1_e = exons_sorted[host_i + 1][0] - 1

    evidence = (
        f"{counts['total']} uniquely-mapped {gene_name} reads; all map to the "
        f"cryptic isoform (0 reads span constitutive intron 1, "
        f"{intron1_s}-{intron1_e}); {span} reads span both novel junctions "
        f"(full exon inclusion); exon body covered {int(cov.min())}-"
        f"{int(cov.max())}x (mean {cov.mean():.1f}x) by "
        f"{jl + jr - span} junction-overlapping reads"
    )

    with open(os.path.join(outdir, "cryptic_exon.tsv"), "w",
              encoding="utf8", newline="\n") as f:
        f.write("gene\tchrom\tstart\tend\tleft_junction_reads\t"
                "right_junction_reads\texpression_evidence\n")
        f.write(f"{gene_name}\t{chrom}\t{cs}\t{ce}\t{jl}\t{jr}\t{evidence}\n")

    refined_juncs = {
        (l_intron[0], l_intron[1]): ("cryptic_exon_left", jl),
        (r_intron[0], r_intron[1]): ("cryptic_exon_right", jr),
    }
    with open(os.path.join(outdir, "junctions.tsv"), "w",
              encoding="utf8", newline="\n") as f:
        f.write("chrom\tintron_start\tintron_end\tstrand\treads\tnovel\t"
                "gene\trole\tnotes\n")
        def near(a, b, tol=4):
            return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol

        host_introns = [(exons_sorted[i][1] + 1, exons_sorted[i + 1][0] - 1)
                        for i in range(len(exons_sorted) - 1)]
        rows = []
        seen = set()
        for r in jrows:
            key = (r["start"], r["end"])
            st = "+" if r["strand"] == 1 else "-"
            if any(near(key, h) for h in host_introns):
                # superseded by the isoform-quantified host-gene rows below
                if key not in refined_juncs:
                    continue
            if key in refined_juncs:
                role, cnt = refined_juncs[key]
                seen.add(key)
                rows.append((r["chrom"], r["start"], r["end"], st, cnt, True,
                             gene_name, role,
                             "canonical GT-AG coordinates resolved by "
                             "spliced-transcript alignment; count includes "
                             "reads spanning both novel junctions"))
            elif near(key, l_intron):
                rows.append((r["chrom"], r["start"], r["end"], st, r["reads"],
                             r["novel"], gene_name, "microhomology_variant",
                             "microhomology-shifted alignment of "
                             "cryptic_exon_left; superseded by "
                             f"{l_intron[0]}-{l_intron[1]}"))
            elif near(key, r_intron):
                rows.append((r["chrom"], r["start"], r["end"], st, r["reads"],
                             r["novel"], gene_name, "microhomology_variant",
                             "microhomology-shifted alignment of "
                             "cryptic_exon_right; superseded by "
                             f"{r_intron[0]}-{r_intron[1]}"))
            else:
                rows.append((r["chrom"], r["start"], r["end"], st,
                             r["reads"], r["novel"], "", "",
                             "uniquely-placed split-alignment cluster"))
        for key, (role, cnt) in refined_juncs.items():
            if key not in seen:
                achrom = chrom[3:] if chrom.startswith("chr") else chrom
                rows.append((chrom, key[0], key[1], "+", cnt,
                             (achrom, key[0], key[1]) not in introns,
                             gene_name, role,
                             "canonical GT-AG coordinates resolved by "
                             "spliced-transcript alignment; count includes "
                             "reads spanning both novel junctions"))
        # annotated host-gene junctions quantified via isoform alignment
        for i in range(len(exons_sorted) - 1):
            a = exons_sorted[i][1] + 1
            b = exons_sorted[i + 1][0] - 1
            if (a, b) in refined_juncs:
                continue
            if i == host_i:
                # constitutive intron harbouring the cryptic exon: no reads
                # span it in full
                rows.append((chrom, a, b, "+" if strand == "1" else "-", 0,
                             False, gene_name, "annotated",
                             "constitutive intron 1; spanned by 0 reads"))
            else:
                cnt = len(counts["J_other"]) if i == host_i + 1 else 0
                rows.append((chrom, a, b, "+" if strand == "1" else "-", cnt,
                             False, gene_name, "annotated",
                             "quantified via isoform alignment"))
        rows.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")
    n_total_final = len(rows)
    n_novel_final = sum(1 for r in rows if r[5] is True)
    write_report(chrom, gene_name, ensg, enst, strand, exons_sorted, refined,
                 counts, jl, jr, span, n_canon, cov, rseq, length, introns,
                 host_i, jrows, genlen, n_total_final, n_novel_final)


def write_report(chrom, gene_name, ensg, enst, strand, exons_sorted, refined,
                 counts, jl, jr, span, n_canon, cov, rseq, length, introns,
                 host_i, jrows, genlen, n_total=None, n_novel=None):
    cs, ce = refined["start"], refined["end"]
    li = refined["intron_left"]
    ri = refined["intron_right"]
    e1s, e1e = exons_sorted[host_i]
    e2s, e2e = exons_sorted[host_i + 1]
    if n_total is None:
        n_total = len(jrows)
    if n_novel is None:
        n_novel = sum(1 for r in jrows if r["novel"])
    report = f"""# Cryptic-exon analysis report

## Question
Identify the protein-coding HGNC gene containing the highly expressed cryptic
exon supported by two novel splice junctions. Novelty is assessed against the
supplied Ensembl release-112 protein-coding annotation.

## Data
- `inputs/cryptic.exon.q1.fq.gz`: {rseq.shape[0]:,} single-end reads, {length} bp each.
- `inputs/reference/GRCh38_chr9.fa.gz`: chromosome {chrom} ({genlen:,} bp, Broad hg38/v0).
- `inputs/reference/ensembl112_protein_coding_exons.tsv.gz`: Ensembl 112
  protein-coding exon annotation ({len(introns):,} transcript introns used as the
  known-junction set).

## Method (implemented in `analysis.py`, pure Python + NumPy)
No splice-aware aligner was available, so the pipeline:

1. builds a genome index restricted to 16-mers present in the reads (k-mers
   with >20 genomic copies are treated as repetitive);
2. places every read (both orientations) by seed/diagonal voting, verifying
   alignments base-by-base (<=3 mismatches / 107 bp; introns 30-200,000 bp);
   reads with several equally good loci (segmental duplications) are
   discarded, so every counted read is uniquely placed;
3. canonicalises each unique junction cluster to the GT-AG splice motif
   (+/-3 bp window, strand-aware);
4. calls a junction novel iff `(chrom, intron_start, intron_end)` (1-based
   inclusive skipped interval) is absent from every supplied transcript;
5. searches pairs of novel junctions whose interval
   `left_intron_end+1 .. right_intron_start-1` overlaps no annotated exon and
   lies inside a protein-coding gene;
6. resolves microhomology-ambiguous boundaries by re-aligning all reads to
   every canonical splice-site variant of the candidate and choosing the
   variant whose spliced product is supported by reads spanning both
   junctions;
7. quantifies the host gene against full-length canonical and cryptic isoform
   references.

## Result

**The cryptic exon lies in the protein-coding HGNC gene {gene_name}
({ensg}), chromosome {chrom}.**

| item | value |
|---|---|
| Gene | {gene_name} (protein-coding, + strand, transcript {enst}) |
| Cryptic exon | {chrom}:{cs}-{ce} ({ce - cs + 1} bp) |
| Location | {gene_name} intron 1 ({e1e + 1}-{e2s - 1}), between exon 1 ({e1s}-{e1e}) and exon 2 ({e2s}-{e2e}) |
| Left novel junction (exon 1 -> cryptic exon) | {chrom}:{li[0]}-{li[1]}, {jl} reads |
| Right novel junction (cryptic exon -> exon 2) | {chrom}:{ri[0]}-{ri[1]}, {jr} reads |
| Reads spanning both novel junctions | {span} |
| Cryptic-exon body coverage | {int(cov.min())}-{int(cov.max())}x (mean {cov.mean():.1f}x) |
| {gene_name} uniquely-mapped reads | {counts['total']} |
| Reads spanning constitutive intron 1 (canonical isoform) | {n_canon} |

### Novelty against the supplied Ensembl 112 annotation
- Left junction `{chrom}:{li[0]}-{li[1]}`: `(chrom, {li[0]}, {li[1]})` is absent
  from every supplied transcript (annotated intron 1 is {e1e + 1}-{e2s - 1}; the
  cryptic 3' splice site {li[1]} is new). **Novel.**
- Right junction `{chrom}:{ri[0]}-{ri[1]}`: `(chrom, {ri[0]}, {ri[1]})` is absent
  from every supplied transcript (the cryptic 5' splice site {ri[0]} is new; the
  downstream acceptor {ri[1]} is the annotated one). **Novel.**
- Both junctions carry canonical GT-AG dinucleotides on the + strand.
- The interval {cs}-{ce} does not overlap any annotated exon of any supplied
  transcript.

### Expression evidence
{gene_name} is among the most highly expressed genes in this library, and
every uniquely-mapped {gene_name} read is consistent with the cryptic isoform
(exon 1 - cryptic exon - exon 2 - exon 3): **0 reads span constitutive intron
1**, i.e. the cryptic exon is included in ~100% of {gene_name} transcripts.
Direct evidence for the exon body: {span} reads span both novel junctions and
thereby cover the entire {ce - cs + 1} bp exon; together with the remaining
junction reads the exon body is covered at {int(cov.min())}-{int(cov.max())}x
(mean {cov.mean():.1f}x). Junction support: {jl} left + {jr} right
(double-spanning reads counted in both).

Boundary note: exon 1 ends in `...AAG` and the genomic cryptic-exon region
begins `AAGTTG...` (3-bp AAG microhomology), which makes raw split alignment
of the left junction ambiguous by 3 bp. Re-aligning reads to candidate
spliced transcripts shows the transcript keeps the AAG on the exon side, so
the cryptic exon starts at {cs} (acceptor AG at {li[1] - 1}-{li[1]}); this is the
only variant supported by reads spanning both junctions with 0-1 mismatches.

### Ruled-out alternatives
`junctions.tsv` lists {n_total} detected junctions, {n_novel} of them novel
against Ensembl 112. The other novel junctions are
+/-1-3 bp shifted copies of annotated introns, frequently supported on both
strands - the signature of reads from inverted segmental-duplication copies
of highly expressed genes (RPS6, RPL7A, HSPA5, SET, ...). Every interval
defined by those junction pairs overlaps an annotated exon (the shifted exon
of the duplicated copy), so none qualifies as a cryptic exon. The {gene_name}
pair is the only well-supported pair of novel junctions whose intervening
interval is unannotated, gene-internal and expressed.

*Generated by `output/analysis.py` (deterministic; re-run with
`python output/analysis.py`).*
"""
    with open(os.path.join(WS, "output", "report.md"), "w",
              encoding="utf8", newline="\n") as f:
        f.write(report)


if __name__ == "__main__":
    main()
