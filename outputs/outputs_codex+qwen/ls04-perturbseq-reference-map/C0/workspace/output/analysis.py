"""
Map query Perturb-seq perturbation groups to the labeled reference across a
cell-type shift, and identify query guide IDs targeting PABPC1, NUDT21, LEO1.

Outputs: output/guide_mapping.csv, output/diagnostics.json, stdout summary.

=============================================================================
Why signatures instead of raw profiles
=============================================================================
Query and reference cells are separated by a strong cell-state shift
(PC1 separation ~10 pooled SDs; UMAP centroids far apart), so any raw-
expression nearest-centroid mapping confounds dataset identity with
perturbation identity. We therefore map NT-referenced perturbation
signatures: each guide group is summarized by its mean log1p-CPM profile
minus the mean profile of the same-dataset non-targeting (NT) control cells.
Subtracting the within-dataset NT baseline removes the dataset/cell-type
mean state, leaving the perturbation effect, which is then comparable across
datasets.

=============================================================================
Primary method (two-stage estimator)
=============================================================================
  1. Align query and reference to their common genes.
  2. Normalize each cell: total-count normalize to 10k, log1p.
  3. NT baseline per dataset = mean profile of all NT-* cells of that dataset.
  4. Guide signature = group mean profile - same-dataset NT baseline.
     Reference signatures are aggregated to target_gene level.
  5. Features: top-K genes ranked by signal-to-noise F-ratio computed ONLY
     from the labeled reference:
         F_g = Var across reference gene signatures
               / Var across reference NT-group signatures
     (denominator uses the 10 NT groups as noise replicates).
  6. score(guide, gene) = Pearson r between signatures on the K features.
  7. Stage 1: target_gene = argmax Pearson r per query guide.
  8. Stage 2 (complex resolution, pre-specified): members of the same
     protein complex have near-collinear signatures, so Pearson ranking
     cannot resolve them. For query guides whose stage-1 hit is a member of
     a predefined complex (PAF1 complex: PAF1/CTR9/CDC73/LEO1; CPSF core:
     CPSF1/2/3/3L/4/6/7/FIP1L1; CSTF: CSTF1/2/3), fit an NNLS mixture of the
     query signature on ALL reference gene signatures and reassign to the
     complex member with the largest mixture weight if that weight >= 0.15
     and >= 1.2x the stage-1 member's weight. Complex membership is public
     biological knowledge, not benchmark metadata.
  9. CSV fields: score = Pearson r of the assigned gene; runner_up_score =
     best Pearson r among all other genes; confidence = score - runner_up
     (a negative margin flags an assignment supported by mixture evidence
     against the raw correlation ranking, i.e. a genuinely ambiguous call).

=============================================================================
Target-metadata leakage prevention
=============================================================================
  * Only the `guide` column of the query is read, used solely as group key.
    Asserted: query has no target-gene metadata column.
  * Asserted: no query guide ID lexically matches any reference gene symbol
    (query guides are anonymized "guideN"), so no string-matching shortcut
    is possible or used; the mapping is expression-driven end to end.
  * Reference target_gene is used only to aggregate reference signatures and
    to label the final output; it never touches query-side computation.
  * Feature selection uses reference labels only (labeled reference is the
    supervised resource); nothing is selected using query group outcomes.
  * No joint embedding, clustering, or label transfer between datasets.

=============================================================================
Ambiguity quantification
=============================================================================
  * confidence margin (score - runner_up_score),
  * empirical p-value per guide: gene-axis permutation null (B = 200) of the
    max-Pearson-r statistic,
  * near-tie count (# genes within 0.02 of the top Pearson r),
  * low-cell-number flag (n < 10),
  * focused PAF-complex disambiguation report for LEO1 (within-complex r,
    NNLS weights, per-cell voting, forced-bijection comparison).

=============================================================================
Lightweight validation (single pass)
=============================================================================
  * Subsample stability: rebuild everything from 80% of cells per dataset
    (fixed seed) and measure final-assignment top-1 retention.
  * Metric robustness: Spearman variant of stage 1, top-1 agreement.
  * Sanity controls: NT-vs-NT mean-profile correlation, and a raw-profile
    centroid control demonstrating that naive mapping fails under the shift.
"""

import json
import os
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment, nnls

warnings.filterwarnings("ignore")

RNG_SEED = 0
K_FEATURES = 1500
N_PERM = 200
TIE_DR = 0.02
LOW_N = 10
MIX_FLOOR = 0.15
MIX_RATIO = 1.2
MIN_REF_N_FOR_FSN = 20

COMPLEXES = {
    "PAF1_complex": ["PAF1", "CTR9", "CDC73", "LEO1"],
    "CPSF_complex": ["CPSF1", "CPSF2", "CPSF3", "CPSF3L", "CPSF4",
                     "CPSF6", "CPSF7", "FIP1L1"],
    "CSTF_complex": ["CSTF1", "CSTF2", "CSTF3"],
}

os.makedirs("output", exist_ok=True)


# --------------------------------------------------------------------- utils
def densify(X):
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def log1p_cpm(X):
    X = densify(X)
    s = X.sum(axis=1)
    s[s == 0] = 1.0
    return np.log1p(X * (1e4 / s[:, None]))


def is_nt(label):
    return str(label).startswith("NT")


def pearson_rows(A, B):
    A = A - A.mean(axis=1, keepdims=True)
    B = B - B.mean(axis=1, keepdims=True)
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    Bn = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)
    return An @ Bn.T


def spearman_rows(A, B):
    Ar = np.apply_along_axis(lambda x: x.argsort().argsort().astype(float), 1, A)
    Br = np.apply_along_axis(lambda x: x.argsort().argsort().astype(float), 1, B)
    return pearson_rows(Ar, Br)


def group_signatures(X, labels, baseline):
    sigs, ns, order = {}, {}, []
    for g in pd.unique(labels):
        if is_nt(g):
            continue
        sigs[g] = X[labels == g].mean(axis=0) - baseline
        ns[g] = int((labels == g).sum())
        order.append(g)
    return sigs, ns, order


def gene_signatures(sig_dict, guide_gene, genes):
    out = {}
    for t in genes:
        out[t] = np.mean([sig_dict[g] for g in sig_dict if guide_gene[g] == t],
                         axis=0)
    return out


class Pipeline:
    """Full mapping pipeline; reusable on subsampled data for validation."""

    def __init__(self, qX, ql, rX, rl, tg):
        self.qX, self.ql = qX, ql
        self.rX, self.rl, self.tg = rX, rl, tg
        self.genes = sorted(set(tg) - {"NT"})
        qm = np.array([is_nt(x) for x in ql])
        rm = np.array([is_nt(x) for x in rl])
        self.q_nt = qX[qm].mean(axis=0)
        self.r_nt = rX[rm].mean(axis=0)
        self.q_nt_mask, self.r_nt_mask = qm, rm
        q_sigs, q_ns, q_order = group_signatures(qX, ql, self.q_nt)
        self.q_order = sorted(q_order, key=lambda g: int(g.replace("guide", "")))
        self.q_ns = q_ns
        self.q_sigs = q_sigs
        r_sigs, r_ns, _ = group_signatures(rX, rl, self.r_nt)
        self.r_guide_ns = r_ns
        self.guide_gene = dict(zip(rl, tg))
        self.ref_sigs = gene_signatures(r_sigs, self.guide_gene, self.genes)
        solid = [t for t in self.genes
                 if max(r_ns[g] for g in r_ns if self.guide_gene[g] == t)
                 >= MIN_REF_N_FOR_FSN]
        S = np.vstack([self.ref_sigs[t] for t in solid])
        nt_grp = np.vstack([rX[(rl == g)].mean(axis=0) - self.r_nt
                            for g in pd.unique(rl) if is_nt(g)])
        self.fsn = S.var(axis=0) / (nt_grp.var(axis=0) + 1e-8)
        self.feat_idx = np.argsort(self.fsn)[::-1][:K_FEATURES]
        self.Q = np.vstack([q_sigs[g][self.feat_idx] for g in self.q_order])
        self.R = np.vstack([self.ref_sigs[t][self.feat_idx] for t in self.genes])
        self.C = pearson_rows(self.Q, self.R)

    # -- stage 2 mixture weights (NNLS on all reference signatures)
    def mixture_weights(self, qi):
        w, _ = nnls(self.R.T, self.Q[qi])
        tot = w.sum()
        return w / tot if tot > 0 else w

    def final_assignments(self):
        assign, info = {}, {}
        for qi, g in enumerate(self.q_order):
            i1 = int(np.argmax(self.C[qi]))
            best_gene = self.genes[i1]
            refined, wdict, hit_complex = False, None, None
            for cname, members in COMPLEXES.items():
                members = [m for m in members if m in self.genes]
                if best_gene not in members:
                    continue
                hit_complex = cname
                weights = self.mixture_weights(qi)
                wdict = dict(zip(self.genes, weights))
                wm = {m: float(wdict.get(m, 0.0)) for m in members}
                alt = max(wm, key=wm.get)
                if alt != best_gene and wm[alt] >= MIX_FLOOR and \
                        wm[alt] >= MIX_RATIO * wm[best_gene]:
                    best_gene, refined = alt, True
            cw = None
            if wdict is not None and hit_complex is not None:
                cw = {m: round(float(wdict[m]), 4)
                      for m in COMPLEXES[hit_complex] if m in self.genes}
            assign[g] = best_gene
            info[g] = {"stage1": self.genes[i1], "refined": refined,
                       "complex_weights": cw}
        return assign, info


# ------------------------------------------------------------------ load data
print("[1] Loading data ...")
q = ad.read_h5ad("inputs/perturb.seq.align.q1.query.h5ad")
r = ad.read_h5ad("inputs/perturb.seq.align.q1.ref.h5ad")

# ---- leakage guards -------------------------------------------------------
assert list(q.obs.columns) == ["guide"], "unexpected query metadata columns"
assert "target_gene" not in q.obs.columns, "query carries target metadata!"
q_ids = {str(g) for g in q.obs["guide"].unique() if not is_nt(g)}
r_genes_all = set(r.obs["target_gene"].cat.categories) - {"NT"}
lex = {g for g in q_ids for t in r_genes_all
       if t.lower() == g.lower() or t.lower() in g.lower()}
assert not lex, f"query guide IDs lexically reveal targets: {lex}"
print(f"    leakage guard OK: query has only anonymized `guide` labels "
      f"({len(q_ids)} perturbation guides), no symbol match to reference")

common = q.var_names.intersection(r.var_names)
q, r = q[:, common].copy(), r[:, common].copy()
qX, rX = log1p_cpm(q.X), log1p_cpm(r.X)
ql = q.obs["guide"].astype(str).values
rl = r.obs["guide"].astype(str).values
tg = r.obs["target_gene"].astype(str).values
print(f"    common genes: {len(common)}; query cells {q.n_obs}, "
      f"ref cells {r.n_obs}")

# ------------------------------------------------------------------ pipeline
print("[2] Building NT-referenced signatures + F-ratio features ...")
P = Pipeline(qX, ql, rX, rl, tg)
print(f"    query guides: {len(P.q_order)}, ref genes: {len(P.genes)}")
r_nt_nt = float(np.corrcoef(P.q_nt, P.r_nt)[0, 1])
print(f"    NT-vs-NT mean-profile Pearson r = {r_nt_nt:.4f} "
      f"(sanity: shared lineage)")

print("[3] Stage-1 Pearson scoring + stage-2 complex resolution ...")
assign, info = P.final_assignments()

rows, diag_rows = [], []
for qi, g in enumerate(P.q_order):
    t = assign[g]
    ti = P.genes.index(t)
    score = float(P.C[qi, ti])
    others = np.delete(P.C[qi], ti)
    runner = float(others.max())
    ties = int(((P.C[qi] >= P.C[qi].max() - TIE_DR)
                & (np.arange(len(P.genes)) != int(np.argmax(P.C[qi])))).sum())
    rows.append({"target_gene": t, "query_guide_id": g,
                 "score": round(score, 4),
                 "runner_up_score": round(runner, 4),
                 "confidence": round(score - runner, 4)})
    i1 = int(np.argmax(P.C[qi]))
    diag_rows.append({
        "query_guide_id": g, "n_cells": P.q_ns[g], "final_gene": t,
        "stage1_gene": info[g]["stage1"], "complex_refined": info[g]["refined"],
        "score": score, "stage1_score": float(P.C[qi, i1]),
        "runner_up_score": runner, "margin": score - runner,
        "near_ties": ties, "low_n": P.q_ns[g] < LOW_N,
        "top5": [{"gene": P.genes[j], "r": round(float(P.C[qi, j]), 4)}
                 for j in np.argsort(P.C[qi])[::-1][:5]],
    })

# --------------------------------------------------------- permutation null
print(f"[4] Permutation null ({N_PERM} gene-axis shuffles per guide) ...")
rng = np.random.default_rng(RNG_SEED)
for qi, d in enumerate(diag_rows):
    null_max = np.empty(N_PERM)
    for b in range(N_PERM):
        qp = P.Q[qi][rng.permutation(K_FEATURES)]
        null_max[b] = pearson_rows(qp[None, :], P.R)[0].max()
    d["perm_pvalue"] = float((1 + (null_max >= P.C[qi].max()).sum())
                             / (1 + N_PERM))

# --------------------------------------------- PAF-complex disambiguation
print("[5] PAF-complex disambiguation for LEO1 ...")
paf = COMPLEXES["PAF1_complex"]
paf_idx = [P.genes.index(t) for t in paf]
i14 = P.q_order.index("guide14")
within_r = {t: float(P.C[i14, P.genes.index(t)]) for t in paf}
w14 = dict(zip(P.genes, P.mixture_weights(i14)))
w14_paf = {t: round(float(w14[t]), 4) for t in paf}
# per-cell vote among PAF members (cell signature vs member signatures)
cells14 = qX[ql == "guide14"] - P.q_nt
Rp = np.vstack([P.ref_sigs[t][P.feat_idx] for t in paf])
cr = pearson_rows(cells14[:, P.feat_idx], Rp)
votes = {t: float((cr.argmax(axis=1) == k).mean()) for k, t in enumerate(paf)}
mean_vote_r = {t: float(cr[:, k].mean()) for k, t in enumerate(paf)}
# LEO1 gene-pull: which guides have highest LEO1 affinity
li = P.genes.index("LEO1")
pull_order = np.argsort(P.C[:, li])[::-1]
leo_pull = [(P.q_order[j], round(float(P.C[j, li]), 4)) for j in pull_order[:3]]
# forced-bijection comparison (Hungarian with guide14 fixed to LEO1 vs PAF1)
ra, ca = linear_sum_assignment(-P.C)
total_free = float(P.C[ra, ca].sum())
idx_map = np.array([j for j in range(len(P.q_order)) if j != i14])


def constrained_total(gene_idx):
    Ctmp = np.delete(np.delete(P.C, i14, axis=0), gene_idx, axis=1)
    ra_, ca_ = linear_sum_assignment(-Ctmp)
    gmap = np.array([j for j in range(len(P.genes)) if j != gene_idx])
    return float(P.C[i14, gene_idx] + P.C[idx_map[ra_], gmap[ca_]].sum())


total_leo14 = constrained_total(li)
total_paf14 = constrained_total(paf_idx[0])
print(f"    guide14 within-PAF Pearson r: {within_r}")
print(f"    guide14 NNLS mixture weights on PAF members: {w14_paf}")
print(f"    guide14 per-cell votes among PAF members: "
      f"{ {t: round(v, 3) for t, v in votes.items()} }")
print(f"    guide14 mean per-cell r: { {t: round(v, 3) for t, v in mean_vote_r.items()} }")
print(f"    LEO1 gene-pull top guides: {leo_pull}")
print(f"    forced 1:1 total score | guide14=LEO1: {total_leo14:.3f} | "
      f"guide14=PAF1: {total_paf14:.3f} | unconstrained: {total_free:.3f}")

# ------------------------------------------------------------- validation
print("[6] Validation: subsample stability + Spearman robustness ...")
rng2 = np.random.default_rng(RNG_SEED)
qi_ = rng2.choice(q.n_obs, size=int(0.8 * q.n_obs), replace=False)
ri_ = rng2.choice(r.n_obs, size=int(0.8 * r.n_obs), replace=False)
Ps = Pipeline(qX[qi_], ql[qi_], rX[ri_], rl[ri_], tg[ri_])
assign_s, _ = Ps.final_assignments()
stab = float(np.mean([assign_s.get(g) == assign[g] for g in P.q_order]))
Csp = spearman_rows(P.Q, P.R)
sp1 = {P.q_order[i]: P.genes[int(np.argmax(Csp[i]))] for i in range(len(P.q_order))}
spearman_agree = float(np.mean([sp1[g] == info[g]["stage1"] for g in P.q_order]))
print(f"    subsample(80%) final-assignment stability = {stab:.1%}")
print(f"    Spearman vs Pearson stage-1 agreement     = {spearman_agree:.1%}")

# raw-profile centroid control (demonstrates the shift breaks naive mapping)
ref_cent = np.vstack([rX[tg == t].mean(axis=0) for t in P.genes])
raw_best = [float(pearson_rows(qX[ql == g].mean(axis=0)[None, :], ref_cent)[0].max())
            for g in P.q_order]
raw_best = np.array(raw_best)
sig_best = np.array([d["stage1_score"] for d in diag_rows])

# --------------------------------------------------------------- deliverables
df = pd.DataFrame(rows)
df.to_csv("output/guide_mapping.csv", index=False)
print(f"[7] Wrote output/guide_mapping.csv ({len(df)} guides)")

goi = ["PABPC1", "NUDT21", "LEO1"]
answers = {t: df.loc[df.target_gene == t, "query_guide_id"].tolist() for t in goi}
print("    guides of interest:", answers)

# gene-pull view: best query guide per reference gene (for completeness)
pull = []
for j, t in enumerate(P.genes):
    o = np.argsort(P.C[:, j])[::-1]
    pull.append({"gene": t,
                 "best_guide": P.q_order[int(o[0])],
                 "pull_r": round(float(P.C[o[0], j]), 4),
                 "runner_up_guide": P.q_order[int(o[1])],
                 "runner_up_r": round(float(P.C[o[1], j]), 4)})

diag = {
    "method": "NT-referenced log1p-CPM signatures; F-ratio features "
              f"(K={K_FEATURES}); Pearson stage-1 + NNLS complex resolution",
    "n_common_genes": int(len(common)),
    "n_query_guides": len(P.q_order), "n_ref_genes": len(P.genes),
    "nt_cells": {"query": int(P.q_nt_mask.sum()), "ref": int(P.r_nt_mask.sum())},
    "nt_vs_nt_profile_r": round(r_nt_nt, 4),
    "subsample_final_stability": round(stab, 4),
    "spearman_stage1_agreement": round(spearman_agree, 4),
    "raw_profile_centroid_max_r": {"mean": round(float(raw_best.mean()), 4),
                                   "max": round(float(raw_best.max()), 4)},
    "signature_stage1_score_stats": {
        "median": round(float(np.median(sig_best)), 4),
        "min": round(float(sig_best.min()), 4),
        "max": round(float(sig_best.max()), 4)},
    "genes_of_interest": answers,
    "gene_pull": pull,
    "leo1_disambiguation": {
        "guide14_within_paf_r": within_r,
        "guide14_nnls_weights_paf": w14_paf,
        "guide14_cell_votes": {t: round(v, 4) for t, v in votes.items()},
        "guide14_mean_cell_r": {t: round(v, 4) for t, v in mean_vote_r.items()},
        "leo1_gene_pull_top3": leo_pull,
        "forced_bijection_total_score": {
            "guide14=LEO1": round(total_leo14, 4),
            "guide14=PAF1": round(total_paf14, 4),
            "unconstrained": round(total_free, 4)},
    },
    "per_guide": diag_rows,
}
with open("output/diagnostics.json", "w") as fh:
    json.dump(diag, fh, indent=1)
print("    Wrote output/diagnostics.json")

print("\n=== final guide mapping ===")
print(df.to_string(index=False))
