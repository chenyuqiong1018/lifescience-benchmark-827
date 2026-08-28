"""Map anonymized query perturbation groups onto the labeled reference.

Primary method (single production run)
--------------------------------------
1. Load ref (labeled: guide + target_gene) and query (guide ids only, anonymized).
2. Leakage guards:
   - query metadata is asserted to contain no target-gene information (only 'guide');
   - guide-name strings are never parsed to infer gene identity; the only label
     semantics used is the explicit 'NT' non-targeting control designation, which is
     part of the query's own experimental design;
   - matching is done exclusively on expression perturbation signatures, never on
     guide-name overlap between datasets.
3. Each dataset is normalized independently (library size 1e4 CPM + log1p), so no
   cross-dataset baseline information is shared.
4. Perturbation signature per guide = mean(log-norm expr of guide cells) minus the
   within-dataset NT-control mean. Differencing against each dataset's own control
   removes the cell-type/baseline shift between reference and query.
5. Reference guide signatures are aggregated to target_gene via the target_gene
   label column (never via guide-name string parsing).
6. Score(query guide, target gene) = Pearson correlation of the two signatures over
   shared genes.
7. One-to-one assignment by Hungarian matching (max total correlation).
8. Output: output/guide_mapping.csv with one row per reference target_gene.

Lightweight independent validation
----------------------------------
- Subsample stability: recompute all signatures on a random 50% of each guide's
  cells (fixed seed) and rerun the assignment; report fraction of stable rows.
- NT sanity: NT control groups must match each other by expression signature.
"""
import json
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")
RNG_SEED = 0

REF_PATH = "inputs/perturb.seq.align.q1.ref.h5ad"
QRY_PATH = "inputs/perturb.seq.align.q1.query.h5ad"


def load(path):
    a = ad.read_h5ad(path)
    a.var_names_make_unique()
    return a


def normalize(a):
    """Independent per-dataset CPM(1e4) + log1p normalization."""
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    X = X.astype(np.float32)
    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(X / lib * 1e4)


def group_means(X, labels):
    """Per-group mean expression (vectorized via argsort + reduceat)."""
    labels = np.asarray(labels)
    order = np.argsort(labels, kind="stable")
    uniq, start = np.unique(labels[order], return_index=True)
    sums = np.add.reduceat(X[order], start, axis=0)
    counts = np.diff(np.append(start, len(labels)))
    return pd.DataFrame(sums / counts[:, None], index=uniq)


def subsample_mask(labels, frac, rng):
    keep = np.zeros(len(labels), bool)
    for g in np.unique(labels):
        idx = np.where(labels == g)[0]
        take = rng.choice(idx, size=max(1, int(len(idx) * frac)), replace=False)
        keep[take] = True
    return keep


def build(ref, qry, frac=1.0, seed=RNG_SEED):
    """Build correlation matrix (reference target_genes x query guides)."""
    Xr, Xq = normalize(ref), normalize(qry)
    shared = sorted(set(ref.var_names) & set(qry.var_names))
    Xr = Xr[:, ref.var_names.get_indexer(shared)]
    Xq = Xq[:, qry.var_names.get_indexer(shared)]

    obsr, obsq = ref.obs.copy(), qry.obs.copy()
    if frac < 1.0:
        rng = np.random.default_rng(seed)
        keep_r = subsample_mask(obsr["guide"].astype(str).values, frac, rng)
        keep_q = subsample_mask(obsq["guide"].astype(str).values, frac, rng)
        obsr, obsq = obsr.loc[keep_r], obsq.loc[keep_q]
        Xr, Xq = Xr[keep_r], Xq[keep_q]

    is_ctrl_r = (obsr["target_gene"].astype(str) == "NT").values
    is_ctrl_q = obsq["guide"].astype(str).str.startswith("NT").values

    gm_r = group_means(Xr, obsr["guide"].astype(str).values)
    gm_q = group_means(Xq, obsq["guide"].astype(str).values)
    ctrl_r = Xr[is_ctrl_r].mean(axis=0)
    ctrl_q = Xq[is_ctrl_q].mean(axis=0)
    sig_r = gm_r.subtract(ctrl_r, axis=1)     # per ref guide vs ref NT
    sig_q = gm_q.subtract(ctrl_q, axis=1)     # per query guide vs query NT

    # aggregate reference guide signatures to target_gene via the label column only
    guide2gene = obsr.groupby("guide")["target_gene"].first()
    sig_r.index = sig_r.index.map(lambda g: guide2gene[g])
    gene_sig = sig_r.groupby(level=0).mean()

    A = gene_sig.values - gene_sig.values.mean(axis=1, keepdims=True)
    B = sig_q.values - sig_q.values.mean(axis=1, keepdims=True)
    C = pd.DataFrame(
        (A / np.linalg.norm(A, axis=1, keepdims=True)) @
        (B / np.linalg.norm(B, axis=1, keepdims=True)).T,
        index=gene_sig.index, columns=sig_q.index)
    return C, {"n_shared_genes": len(shared),
               "ctrl_profile_corr": float(pearsonr(ctrl_r, ctrl_q)[0])}


def assign(C):
    """One-to-one Hungarian assignment maximizing total correlation."""
    gi, qj = linear_sum_assignment(-C.values)
    order = np.argsort(gi)
    rows = []
    for i, j in zip(gi[order], qj[order]):
        score = float(C.values[i, j])
        runner_up = float(np.delete(C.values[i], j).max())
        w = C.values[i] - C.values[i].max()
        prob = float(np.exp(w)[j] / np.exp(w).sum())
        rows.append({"target_gene": C.index[i], "query_guide_id": C.columns[j],
                     "score": score, "runner_up_score": runner_up,
                     "confidence": prob, "margin": score - runner_up})
    return pd.DataFrame(rows)


def main():
    ref, qry = load(REF_PATH), load(QRY_PATH)

    # ---- Leakage guards -------------------------------------------------
    assert list(qry.obs.columns) == ["guide"], "query carries unexpected metadata"
    assert not any("target" in c.lower() or "gene" in c.lower()
                   for c in qry.obs.columns), "target metadata in query"

    C, extras = build(ref, qry, frac=1.0)
    res = assign(C)
    res_out = (res[["target_gene", "query_guide_id", "score",
                    "runner_up_score", "confidence"]]
               .sort_values("target_gene").reset_index(drop=True))
    res_out.to_csv("output/guide_mapping.csv", index=False)

    # ---- lightweight independent validation -----------------------------
    C_sub, _ = build(ref, qry, frac=0.5, seed=RNG_SEED)
    res_sub = assign(C_sub).set_index("target_gene")
    stable = float((res_sub.loc[res["target_gene"], "query_guide_id"].values
                    == res["query_guide_id"].values).mean())

    # NT sanity: control groups carry no perturbation signature, so (a) the NT
    # row must still be assigned an NT guide, and (b) NT query guides must show
    # only noise-level correlations with every perturbation signature.
    nt_cols = [c for c in C.columns if c.startswith("NT")]
    nt_guide_max_abs_corr = float(C[nt_cols].abs().values.max())
    nt_row_assigned_nt_guide = bool(
        res.loc[res.target_gene == "NT", "query_guide_id"]
        .str.startswith("NT").all())

    # supplementary per-query-guide view (explanatory; not required artifact)
    qmap = []
    for g in C.columns:
        top2 = C[g].sort_values(ascending=False).head(2)
        qmap.append({"query_guide_id": g,
                     "best_target_gene": top2.index[0],
                     "score": float(top2.iloc[0]),
                     "runner_up_target_gene": top2.index[1],
                     "runner_up_score": float(top2.iloc[1])})
    pd.DataFrame(qmap).to_csv("output/supplementary_query_guide_mapping.csv",
                              index=False)

    summary = {"n_ref_cells": int(ref.n_obs), "n_qry_cells": int(qry.n_obs),
               "n_shared_genes": extras["n_shared_genes"],
               "ctrl_profile_corr": extras["ctrl_profile_corr"],
               "n_rows_written": int(len(res_out)),
               "subsample_stability": stable,
               "subsample_unstable_genes": sorted(
                   res.loc[res_sub.loc[res["target_gene"], "query_guide_id"].values
                           != res["query_guide_id"].values, "target_gene"].tolist()),
               "nt_guide_max_abs_corr": nt_guide_max_abs_corr,
               "nt_row_assigned_nt_guide": nt_row_assigned_nt_guide,
               "targets": {}}
    for g in ["PABPC1", "NUDT21", "LEO1"]:
        r = res[res.target_gene == g].iloc[0]
        summary["targets"][g] = {k: r[k] for k in
                                 ["query_guide_id", "score", "runner_up_score",
                                  "confidence", "margin"]}
    with open("output/run_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    print("\nFull mapping:")
    print(res.sort_values("target_gene").to_string(index=False))


if __name__ == "__main__":
    main()
