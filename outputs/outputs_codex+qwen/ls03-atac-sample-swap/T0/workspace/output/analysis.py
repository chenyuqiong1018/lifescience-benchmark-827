#!/usr/bin/env python3
"""
Detect organ-label swaps in the axolotl bulk ATAC-seq 10-kb bin count table.

Approach
--------
1. Parse the AmexG v6.0 GTF and index genes by ortholog symbol tokens found in
   gene_name ("SYM [hs]", "SYM [nr]" or plain "SYM").
2. For each of the 15 organs, define a curated set of marker genes (human/
   vertebrate ortholog symbols with well-established organ specificity).
3. Stream the bin x sample count table once. For every marker gene accumulate
   promoter counts (TSS +/- 2 kb) and local flank counts (TSS +/- 100 kb,
   excluding the promoter bins), plus per-sample library sizes.
4. Two independent accessibility metrics per gene/sample:
     a) log2 promoter CPM (per 10 kb bin), z-scored across the 15 samples
     b) log2 (promoter-per-bin / flank-per-bin), z-scored across samples
   Organ score(sample) = mean z across its marker genes; the consensus of (a)
   and (b) is the primary score.
5. Marker coherence: for each organ, which labeled sample best matches its
   marker set (and vice versa).
6. Pairwise swap evidence for ordered pair (A -> B):
     e(A->B) = score(A markers, column B) - score(A markers, column A)
   swap_score({A,B}) = e(A->B) + e(B->A). All 105 unordered pairs are ranked.
7. Genome-wide log2-CPM Pearson correlation between samples (two streaming
   passes) as an independent QC of gross sample anomalies.

Decision rule: swap_detected = true only if the top-ranked pair has both
directional terms positive, a clear margin over the runner-up pair, and the
two organs' marker sets mutually prefer each other's labeled columns.
"""
import gzip
import itertools
import json
import math
import os
import re
import collections

TSS_WIN = 2_000        # promoter window (bp) around TSS
FLANK = 100_000        # local background window (bp) around TSS
MAX_HITS_PER_SYMBOL = 6  # symbols matching more genes than this are ambiguous
EPS = 0.05             # CPM pseudocount

ORGANS = ['Bladder', 'Brain', 'Cloaca', 'GallBladder', 'Gill', 'Heart',
          'Intestine', 'Kidney', 'Limb', 'Liver', 'Lung', 'Pancreas',
          'Prostate', 'Spleen', 'Stomach']

MARKERS = {
    'Bladder': ['UPK1A', 'UPK1B', 'UPK2', 'UPK3A', 'UPK3B', 'KRT20', 'KRT8',
                'KRT18', 'PPARG', 'SLC26A6', 'ELF3', 'FGFR3', 'PLAC8', 'GRHL3'],
    'Brain': ['SNAP25', 'SYT1', 'SYN1', 'SYP', 'NCAM1', 'MAP2', 'TUBB3',
              'NEFL', 'NEFM', 'NEFH', 'ENO2', 'GAD1', 'GAD2', 'GRIN1',
              'GRIN2B', 'GFAP', 'MBP', 'PLP1', 'MOG', 'OLIG2', 'KCNA2',
              'DLG4', 'CAMK2A', 'GRIA2', 'SYT4'],
    'Cloaca': ['CDX2', 'HOXC10', 'HOXC11', 'HOXC12', 'HOXC13', 'WNT5A',
               'TP63', 'GATA3', 'MSX2'],
    'GallBladder': ['KRT7', 'KRT19', 'KRT8', 'KRT18', 'SOX9', 'EPCAM',
                    'ANXA4', 'SLC26A6', 'ONECUT1', 'HNF1B', 'MUC1', 'CFTR',
                    'SLC4A2'],
    'Gill': ['RHBG', 'RHCG', 'ATP1A1', 'ATP1B1', 'ATP6V1A', 'CA2', 'CA4',
             'CA12', 'SLC4A4', 'SLC4A1', 'SLC12A2', 'SLC9A3', 'SLC26A6',
             'SLC26A3', 'AQP3'],
    'Heart': ['MYH6', 'MYH7', 'TNNT2', 'TNNI3', 'ACTC1', 'MYL2', 'MYL3',
              'MYL4', 'MYL7', 'NPPA', 'NPPB', 'RYR2', 'ATP2A2', 'KCNJ2',
              'GJA1', 'GJA5', 'TBX5', 'GATA4', 'MEF2C', 'HAND2', 'ACTN2',
              'DES', 'CASQ2', 'SCN5A'],
    'Intestine': ['VIL1', 'CDX2', 'FABP1', 'FABP2', 'FABP6', 'SI', 'LCT',
                  'DPP4', 'MUC2', 'REG4', 'TREH', 'SLC5A1', 'GUCA2A',
                  'APOA4', 'CDX1', 'SLC15A1', 'KRT20'],
    'Kidney': ['AQP1', 'AQP2', 'SLC12A1', 'SLC12A3', 'UMOD', 'CALB1',
               'NPHS1', 'NPHS2', 'PODXL', 'LRP2', 'CUBN', 'SLC22A6',
               'SLC22A8', 'GGT1', 'PAX8', 'WT1', 'KCNJ1', 'CLCNKB',
               'HNF1B', 'SLC34A1'],
    'Limb': ['HOXA9', 'HOXA10', 'HOXA11', 'HOXA13', 'HOXD9', 'HOXD10',
             'HOXD11', 'HOXD12', 'HOXD13', 'MEIS1', 'MEIS2', 'TBX5', 'TBX4',
             'PITX1', 'MSX1', 'MSX2', 'DLX5', 'GREM1', 'FGF8', 'SHH',
             'PBX1', 'HAND2', 'HOXC11'],
    'Liver': ['ALB', 'AFP', 'APOA1', 'APOA2', 'APOB', 'APOC2', 'APOC3',
              'APOE', 'SERPINA1', 'AGT', 'TTR', 'FGA', 'FGB', 'FGG', 'HP',
              'HNF4A', 'HNF1A', 'ONECUT1', 'FGL1', 'CYP3A4', 'CYP2E1',
              'SERPINC1', 'PROC', 'AMBP', 'AHSG', 'ORM1', 'LECT2'],
    'Lung': ['SFTPA1', 'SFTPA2', 'SFTPB', 'SFTPC', 'SFTPD', 'ABCA3',
             'NKX2-1', 'SCGB1A1', 'SCGB3A1', 'SCGB3A2', 'AQP5', 'HOPX',
             'MUC1', 'SLC34A2', 'FXYD3', 'LAMP3', 'AGER'],
    'Pancreas': ['INS', 'GCG', 'SST', 'PPY', 'PDX1', 'MAFA', 'MAFB',
                 'NKX6-1', 'NKX2-2', 'ISL1', 'SLC30A8', 'IAPP', 'PCSK1',
                 'PCSK2', 'REG1A', 'REG1B', 'CELA3A', 'CELA2A', 'CPA1',
                 'CPA2', 'CPB1', 'CTRB1', 'PRSS1', 'PRSS2', 'PNLIP'],
    # TMPRSS3 is deliberately not used as a prostate marker: its expression
    # is not prostate-restricted (inner ear / multiple epithelia), which would
    # cross-contaminate the score in epithelium-rich samples.
    'Prostate': ['KLK3', 'KLK2', 'KLK4', 'NKX3-1', 'HOXB13', 'TMPRSS2',
                 'ACPP', 'MSMB', 'AR', 'STEAP4', 'SLC45A3', 'KLK1', 'KLK7',
                 'AZGP1'],
    'Spleen': ['HBA1', 'HBA2', 'HBB', 'HBG1', 'ALAS2', 'TFRC', 'GYPA',
               'PTPRC', 'CD3E', 'CD3D', 'CD79A', 'CD19', 'SPI1', 'LYZ',
               'MPO', 'ELANE', 'CSF1R', 'ITGAM', 'ITGAX', 'RAG1', 'RAG2',
               'GATA1', 'C1QA', 'CD74', 'FCER1G', 'TYROBP'],
    'Stomach': ['GHRL', 'PGA3', 'PGA4', 'PGA5', 'PGC', 'ATP4A', 'ATP4B',
                'MUC5AC', 'MUC6', 'GKN1', 'LIPF', 'GAST', 'HDC', 'CHIA',
                'MUC1', 'TFF1', 'TFF2', 'GIF', 'CLDN18', 'SHH', 'TFF3',
                'TFF3.2'],
}

HALLMARKS = {
    'Cloaca': ['CLDN18', 'SHH', 'TFF3.2', 'MUC2', 'PGA3', 'PGC'],
    'Stomach': ['CDX2', 'HOXC11', 'HOXC10', 'MSMB', 'HOXB13'],
}

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
INPUTS = os.path.join(ROOT, 'inputs')
if not os.path.isdir(INPUTS):
    # Fallback to cwd-relative paths. On some Windows consoles the Python
    # level cwd string can be mojibake-corrupted even though the process
    # working directory is correct, so relative paths still resolve.
    ROOT = '.'
    INPUTS = os.path.join('.', 'inputs')
    BASE = 'output'
GTF = os.path.join(INPUTS, 'AmexT_v47-AmexG_v6.0-DD.gtf.gz')
TSV = os.path.join(INPUTS, 'sample.swap.atac.q1.tsv.gz')
CHROMSIZES = os.path.join(INPUTS, 'sample.swap.atac.q1.chrom.sizes')


def log(msg):
    print(msg, flush=True)


def parse_chrom_pieces():
    pieces = []
    arm_len = collections.defaultdict(int)
    with open(CHROMSIZES) as fh:
        for line in fh:
            p, ln = line.rstrip('\n').split('\t')
            arm = re.sub(r'_\d+$', '', p)
            pieces.append((p, arm, arm_len[arm], int(ln)))
            arm_len[arm] += int(ln)
    return pieces, dict(arm_len)


def parse_genes():
    attr_re = re.compile(r'gene_id "([^"]+)"; gene_name "([^"]+)"')
    genes = []
    with gzip.open(GTF, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if cols[2] != 'gene':
                continue
            m = attr_re.search(cols[8])
            if not m:
                continue
            genes.append((cols[0], int(cols[3]), int(cols[4]), cols[6],
                          m.group(1), m.group(2)))
    return genes


def build_symbol_index(genes):
    idx = collections.defaultdict(list)
    for (c, s, e, st, gid, gname) in genes:
        for tok in gname.split('|'):
            tok = tok.strip()
            m = re.match(r'^(.*)\s+\[(\w+)\]$', tok)
            sym = m.group(1).strip() if m else tok
            idx[sym].append((c, s, e, st, gid))
    return idx


def resolve_markers(sym_index, arm_len):
    resolved, notes = {}, []
    for organ, symbols in MARKERS.items():
        picked = []
        for sym in symbols:
            hits = [h for h in sym_index.get(sym, []) if h[0] in arm_len]
            if not hits:
                continue
            if len(hits) > MAX_HITS_PER_SYMBOL:
                notes.append(f'{organ}: symbol {sym} skipped '
                             f'({len(hits)} hits, ambiguous)')
                continue
            for (c, s, e, st, gid) in hits:
                tss0 = (s - 1) if st == '+' else e
                picked.append(dict(symbol=sym, gid=gid, chrom=c, tss0=tss0,
                                   strand=st))
        resolved[organ] = picked
        n_sym = len({p['symbol'] for p in picked})
        notes.append(f'{organ}: {len(MARKERS[organ])} candidate symbols -> '
                     f'{n_sym} resolved, {len(picked)} genes')
    return resolved, notes


def build_targets(resolved, pieces):
    piece_by_arm = collections.defaultdict(list)
    for (p, arm, off, ln) in pieces:
        piece_by_arm[arm].append((p, off, ln))
    genes_all, organ_gene_idx = [], {}
    gi = 0
    for organ in ORGANS:
        idxs = []
        for g in resolved[organ]:
            genes_all.append(g)
            idxs.append(gi)
            gi += 1
        organ_gene_idx[organ] = idxs
    targets = collections.defaultdict(list)
    for gidx, g in enumerate(genes_all):
        arm, tss0 = g['chrom'], g['tss0']
        pws, pwe = tss0 - TSS_WIN, tss0 + TSS_WIN
        fws, fwe = tss0 - FLANK, tss0 + FLANK
        for (p, off, ln) in piece_by_arm[arm]:
            lo, hi = off, off + ln
            s_, e_ = max(pws, lo), min(pwe, hi)
            if s_ < e_:
                targets[p].append([s_ - off, e_ - off, gidx, 0])
            s_, e_ = max(fws, lo), min(pws, hi)
            if s_ < e_:
                targets[p].append([s_ - off, e_ - off, gidx, 1])
            s_, e_ = max(pwe, lo), min(fwe, hi)
            if s_ < e_:
                targets[p].append([s_ - off, e_ - off, gidx, 1])
    for p in targets:
        targets[p].sort(key=lambda t: t[0])
    return genes_all, organ_gene_idx, targets


def score_marker_genes(targets, n_genes, n_samples):
    prom = [[0] * n_samples for _ in range(n_genes)]
    flank = [[0] * n_samples for _ in range(n_genes)]
    prom_bins = [0] * n_genes
    flank_bins = [0] * n_genes
    colsum = [0] * n_samples
    cur, active, ptr, starts, tlist = None, [], 0, [], []
    with gzip.open(TSV, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        samples = header[3:]
        assert samples == ORGANS, samples
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            piece = parts[0]
            s = int(parts[1]); e = int(parts[2])
            vals = [int(x) for x in parts[3:]]
            for k in range(n_samples):
                colsum[k] += vals[k]
            if piece != cur:
                cur = piece
                tlist = targets.get(piece, [])
                starts = [t[0] for t in tlist]
                active, ptr = [], 0
            if not starts:
                continue
            while ptr < len(starts) and starts[ptr] < e:
                active.append(tlist[ptr]); ptr += 1
            if active:
                active = [t for t in active if t[1] > s]
            if not active:
                continue
            prom_g = {t[2] for t in active
                      if t[3] == 0 and t[0] < e and t[1] > s}
            for t in active:
                if t[0] >= e or t[1] <= s:
                    continue
                gidx, kind = t[2], t[3]
                if kind == 0:
                    for k in range(n_samples):
                        prom[gidx][k] += vals[k]
                    prom_bins[gidx] += 1
                elif gidx not in prom_g:
                    for k in range(n_samples):
                        flank[gidx][k] += vals[k]
                    flank_bins[gidx] += 1
    return prom, flank, prom_bins, flank_bins, colsum, samples


def zscore_rows(rows):
    out = []
    for r in rows:
        mu = sum(r) / len(r)
        var = sum((x - mu) ** 2 for x in r) / len(r)
        sd = math.sqrt(var) if var > 0 else 1.0
        out.append([(x - mu) / sd for x in r])
    return out


def compute_scores(prom, flank, prom_bins, flank_bins, colsum, genes_all,
                   organ_gene_idx):
    n = len(colsum)
    G = len(genes_all)
    log_cpm, log_fold = [], []
    valid = []
    for g in range(G):
        if prom_bins[g] < 1:
            valid.append(False)
            log_cpm.append(None); log_fold.append(None)
            continue
        valid.append(True)
        cpm_row, fold_row = [], []
        for k in range(n):
            p = prom[g][k] / colsum[k] * 1e6 / prom_bins[g]
            f = ((flank[g][k] / colsum[k] * 1e6 / flank_bins[g])
                 if flank_bins[g] else 0.0)
            cpm_row.append(math.log2(p + EPS))
            fold_row.append(math.log2(p / (f + EPS) + 1e-3))
        log_cpm.append(cpm_row); log_fold.append(fold_row)
    z_cpm_rows = zscore_rows([r for r in log_cpm if r is not None])
    z_fold_rows = zscore_rows([r for r in log_fold if r is not None])
    z_cpm = [None] * G; z_fold = [None] * G
    vi = 0
    for g in range(G):
        if valid[g]:
            z_cpm[g] = z_cpm_rows[vi]; z_fold[g] = z_fold_rows[vi]; vi += 1
    score_cpm = {}; score_fold = {}; score = {}
    for organ in ORGANS:
        idx = [g for g in organ_gene_idx[organ] if valid[g]]
        rows_c = [z_cpm[g] for g in idx]
        rows_f = [z_fold[g] for g in idx]
        sc = [sum(r[k] for r in rows_c) / len(rows_c) for k in range(n)]
        sf = [sum(r[k] for r in rows_f) / len(rows_f) for k in range(n)]
        score_cpm[organ] = sc
        score_fold[organ] = sf
        score[organ] = [(a + b) / 2 for a, b in zip(sc, sf)]
    return score, score_cpm, score_fold, valid, z_cpm, z_fold


def pairwise_swap_scores(score):
    E = {}
    for a in ORGANS:
        for b in ORGANS:
            if a != b:
                E[(a, b)] = (score[a][ORGANS.index(b)] -
                             score[a][ORGANS.index(a)])
    pairs = []
    for a, b in itertools.combinations(ORGANS, 2):
        s = E[(a, b)] + E[(b, a)]
        pairs.append(dict(organ_a=a, organ_b=b, swap_score=s,
                          e_ab=E[(a, b)], e_ba=E[(b, a)]))
    pairs.sort(key=lambda d: d['swap_score'], reverse=True)
    for rank, d in enumerate(pairs, 1):
        d['rank'] = rank
    return pairs, E


def genome_wide_correlation():
    colsum = [0.0] * 15
    with gzip.open(TSV, 'rt') as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            for k in range(15):
                colsum[k] += int(parts[3 + k])
    S = [0.0] * 15; Sq = [0.0] * 15
    C = [[0.0] * 15 for _ in range(15)]
    M = 0
    with gzip.open(TSV, 'rt') as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            vals = []
            for k in range(15):
                cpm = int(parts[3 + k]) / colsum[k] * 1e6
                vals.append(math.log2(cpm + 0.1))
            for k in range(15):
                S[k] += vals[k]; Sq[k] += vals[k] * vals[k]
            for i in range(15):
                vi = vals[i]
                row = C[i]
                for j in range(i, 15):
                    row[j] += vi * vals[j]
            M += 1
    corr = [[0.0] * 15 for _ in range(15)]
    mu = [S[k] / M for k in range(15)]
    var = [max(Sq[k] / M - mu[k] ** 2, 1e-12) for k in range(15)]
    for i in range(15):
        for j in range(i, 15):
            cov = C[i][j] / M - mu[i] * mu[j]
            r = cov / math.sqrt(var[i] * var[j])
            corr[i][j] = corr[j][i] = r
    return corr


def main():
    pieces, arm_len = parse_chrom_pieces()
    log(f'chromosome arms: {len(arm_len)}, pieces: {len(pieces)}')
    genes = parse_genes()
    log(f'GTF genes: {len(genes)}')
    sym_index = build_symbol_index(genes)
    resolved, notes = resolve_markers(sym_index, arm_len)
    for nt in notes:
        log('  ' + nt)
    genes_all, organ_gene_idx, targets = build_targets(resolved, pieces)
    log(f'marker genes scored: {len(genes_all)}')
    prom, flank, prom_bins, flank_bins, colsum, samples = score_marker_genes(
        targets, len(genes_all), len(ORGANS))
    log('library sizes: ' + ', '.join(f'{s}={c/1e6:.0f}M'
                                       for s, c in zip(samples, colsum)))
    score, score_cpm, score_fold, valid, z_cpm, z_fold = compute_scores(
        prom, flank, prom_bins, flank_bins, colsum, genes_all, organ_gene_idx)

    log('\norgan marker score matrix (rows=marker set, cols=labeled sample):')
    log('set'.ljust(13) + ''.join(s[:9].rjust(10) for s in ORGANS))
    for organ in ORGANS:
        log(organ.ljust(13) + ''.join(f'{v:10.2f}' for v in score[organ]))

    log('\nmarker-coherence summary:')
    best_sample_of = {}
    for organ in ORGANS:
        vals = score[organ]
        bi = max(range(len(ORGANS)), key=lambda k: vals[k])
        self_rank = sum(1 for v in vals
                        if v > vals[ORGANS.index(organ)]) + 1
        best_sample_of[organ] = ORGANS[bi]
        log(f'  {organ:12s} self_z={vals[ORGANS.index(organ)]:+.2f} '
            f'best_sample={ORGANS[bi]:12s} (z={vals[bi]:+.2f}) '
            f'self_rank={self_rank}/15')
    for k, s in enumerate(ORGANS):
        vals = [score[o][k] for o in ORGANS]
        bi = max(range(len(ORGANS)), key=lambda i: vals[i])
        log(f'  sample {s:12s} best_matches={ORGANS[bi]:12s} '
            f'(z={vals[bi]:+.2f}, own z={vals[k]:+.2f})')

    pairs, E = pairwise_swap_scores(score)
    log('\ntop 10 pairs by swap_score:')
    for d in pairs[:10]:
        log(f"  {d['swap_score']:+.3f}  {d['organ_a']} <-> {d['organ_b']}  "
            f"(e({d['organ_a']}->{d['organ_b']})={d['e_ab']:+.2f}, "
            f"e({d['organ_b']}->{d['organ_a']})={d['e_ba']:+.2f})")

    hallmark_table = {}
    sym2g = collections.defaultdict(list)
    for g in range(len(genes_all)):
        if valid[g]:
            sym2g[genes_all[g]['symbol']].append(g)
    log('\nhallmark gene consensus z-profiles:')
    for organ, syms in HALLMARKS.items():
        for sym in syms:
            for g in sym2g.get(sym, []):
                zc = z_cpm[g]; zf = z_fold[g]
                cons = [(a + b) / 2 for a, b in zip(zc, zf)]
                am = ORGANS[max(range(15), key=lambda k: cons[k])]
                hallmark_table.setdefault(organ, []).append(
                    (genes_all[g]['symbol'], genes_all[g]['gid'], cons, am))
                log('  ' + (sym + '/' + genes_all[g]['gid'][-4:]).ljust(16) +
                    ''.join(f'{v:8.2f}' for v in cons) + f'  -> {am}')

    corr = genome_wide_correlation()
    log('\ngenome-wide sample correlation (Pearson, log2 CPM):')
    log(''.ljust(12) + ''.join(s[:8].rjust(9) for s in ORGANS))
    for i, s in enumerate(ORGANS):
        log(s.ljust(12) + ''.join(f'{corr[i][j]:9.3f}' for j in range(15)))

    top = pairs[0]
    second = pairs[1]
    a, b = top['organ_a'], top['organ_b']
    both_positive = top['e_ab'] > 0 and top['e_ba'] > 0
    mutual_best = (best_sample_of[a] == b and best_sample_of[b] == a)
    displaced = (best_sample_of[a] != a and best_sample_of[b] != b)
    margin = top['swap_score'] - second['swap_score']
    score_ok = top['swap_score'] >= 1.0
    detected = bool(both_positive and score_ok and margin >= 0.2 and
                    mutual_best and displaced)
    if detected:
        confidence = min(0.97,
                         0.60 +
                         0.20 * min(margin / 0.5, 1.0) +
                         0.10 * min(top['swap_score'] / 2.0, 1.0) +
                         0.07 * min(abs(top['e_ab']) + abs(top['e_ba']),
                                    2.0) / 2.0)
        confidence = round(confidence, 2)
        ev_parts = [
            f"column labeled {b} shows the strongest accessibility of {a} "
            f"marker genes (z {score[a][ORGANS.index(b)]:+.2f} vs "
            f"{score[a][ORGANS.index(a)]:+.2f} in column {a})",
            f"column labeled {a} shows the strongest accessibility of {b} "
            f"marker genes (z {score[b][ORGANS.index(a)]:+.2f} vs "
            f"{score[b][ORGANS.index(b)]:+.2f} in column {b})",
        ]
        for organ in (a, b):
            hits = hallmark_table.get(organ, [])
            if hits:
                top3 = sorted(hits, key=lambda h: -max(h[2]))[:3]
                ev_parts.append(
                    'hallmark genes for ' + organ + ': ' +
                    '; '.join(f"{h[0]} peaks in column {h[3]} "
                              f"(z={max(h[2]):+.2f})" for h in top3))
        ev_parts.append('the remaining 13 organs best-match their own '
                        'labeled columns')
        evidence = '; '.join(ev_parts)
    else:
        confidence = 0.2
        evidence = (f"top pair {a}/{b} (score {top['swap_score']:+.2f}) vs "
                    f"next {second['organ_a']}/{second['organ_b']} "
                    f"({second['swap_score']:+.2f}); criteria: "
                    f"both_directions_positive={both_positive}, "
                    f"mutual_best={mutual_best}, displaced={displaced}, "
                    f"margin={margin:.2f}, score_ok={score_ok}")

    os.makedirs(BASE, exist_ok=True)
    call = {
        'swap_detected': detected,
        'organ_a': a,
        'organ_b': b,
        'confidence': confidence,
        'evidence': evidence,
    }
    with open(os.path.join(BASE, 'swap_call.json'), 'w') as fh:
        json.dump(call, fh, indent=2)
    log('\nswap_call.json: ' + json.dumps(call, indent=2))

    with open(os.path.join(BASE, 'sample_similarity.csv'), 'w') as fh:
        fh.write('organ_a,organ_b,swap_score,rank,evidence_type\n')
        for d in pairs:
            ev = 'promoter_accessibility+marker_coherence'
            if d['rank'] == 1 and detected:
                ev += '+hallmark_genes'
            fh.write(f"{d['organ_a']},{d['organ_b']},"
                     f"{d['swap_score']:.6f},{d['rank']},{ev}\n")
    log(f'sample_similarity.csv written ({len(pairs)} pairs)')


if __name__ == '__main__':
    main()
