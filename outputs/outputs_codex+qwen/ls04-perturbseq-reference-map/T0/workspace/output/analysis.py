"""
Map query Perturb-seq guide groups to a labeled reference across a cell-type shift.

Primary method (single production run):
  1. Align query and reference on shared genes (query var intentionally drops most
     reference target genes -> mapping cannot use target-gene expression).
  2. Per-dataset CP10k + log1p, joint HVG selection, joint PCA for description.
  3. Cell-state clusters via KMeans on each dataset's own PCA (cell-type proxy).
  4. Cell-type-adjusted perturbation effect per guide: OLS
     expr_g ~ intercept + cluster dummies + perturbation_indicator
     fitted on that guide's cells + same-dataset NT control cells, per HVG gene.
  5. Query->reference mapping: Pearson correlation of effect vectors on shared HVGs,
     excluding the candidate target gene itself (no cis/target-metadata leakage).
  6. Ambiguity: score, runner-up score, score gap, and nonparametric bootstrap
     agreement (resample cells, refit effects, re-argmax).

Leak prevention:
  - Query guide IDs are treated as opaque labels (never parsed as gene symbols).
  - Similarity never uses the candidate target gene's own expression.
  - Effects are computed separately per dataset; no joint normalization of effects.

Validation (lightweight, in the same run):
  - Reference split-half self-mapping (hold out half the cells of each ref guide,
    re-rank all 33 targets; top-1 accuracy is the oracle-recoverable ceiling).
  - NT negative controls: NT-vs-NT pseudo-perturbations must have low max similarity.

Outputs: output/guide_mapping.csv, output/results.json
"""
import json, warnings, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")
RNG_SEED = 0
rng = np.random.default_rng(RNG_SEED)
N_BOOT = 100
K_CLUSTERS = 10
N_TOP_GENES = 2000

def load_aligned():
    ref = ad.read_h5ad("inputs/perturb.seq.align.q1.ref.h5ad")
    qry = ad.read_h5ad("inputs/perturb.seq.align.q1.query.h5ad")
    shared = sorted(set(ref.var_names) & set(qry.var_names))
    ref, qry = ref[:, shared].copy(), qry[:, shared].copy()
    ref.obs["batch"], qry.obs["batch"] = "ref", "qry"
    for a in (ref, qry):
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    targets = sorted(set(ref.obs["target_gene"]) - {"NT"})
    dropped = sorted(set(ref.obs["target_gene"].dropna()) - set(shared) - {"NT"})
    return ref, qry, shared, targets, dropped

def hvg_pca(ref, qry, shared):
    a = ad.concat([ref, qry], merge="same")
    sc.pp.highly_variable_genes(a, n_top_genes=N_TOP_GENES, flavor="seurat")
    hv = list(a.var_names[a.var["highly_variable"]])
    X = a[:, hv].X.toarray() if hasattr(a[:, hv].X, "toarray") else np.asarray(a[:, hv].X)
    Xs = (X - X.mean(0)) / (X.std(0) + 1e-8)
    n_ref = ref.n_obs
    pca = PCA(n_components=30, random_state=RNG_SEED).fit(Xs)
    emb = pca.transform(Xs)
    return hv, emb[:n_ref], emb[n_ref:], a.obs["batch"].values

def cluster_kmeans(emb, k=K_CLUSTERS):
    return KMeans(n_clusters=k, random_state=RNG_SEED, n_init=10).fit_predict(emb)

def ols_effects(X, labels, cluster, group, nt_mask, hv_index):
    """Cell-type-adjusted effect of `group` vs NT. X: cells x genes (log1p CP10k)."""
    sel = (labels == group) | nt_mask
    Xg, cg, ind = X[sel], cluster[sel], (labels[sel] == group).astype(float)
    C = np.zeros((len(cg), cluster.max() + 1)); C[np.arange(len(cg)), cg] = 1.0
    D = np.column_stack([np.ones(len(cg)), C[:, 1:], ind])
    beta, *_ = np.linalg.lstsq(D, Xg, rcond=None)
    return beta[-1]

def pearson_excl(a, b, excl_idx):
    keep = np.ones(a.shape[0], bool)
    if excl_idx is not None:
        keep[excl_idx] = False
    a, b = a[keep], b[keep]
    a = a - a.mean(); b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else 0.0

def main():
    ref, qry, shared, targets, dropped = load_aligned()
    hv, emb_ref, emb_qry, batch = hvg_pca(ref, qry, shared)
    clu_ref = cluster_kmeans(emb_ref)
    clu_qry = cluster_kmeans(emb_qry)

    Xr = ref[:, hv].X; Xr = Xr.toarray() if hasattr(Xr, "toarray") else np.asarray(Xr)
    Xq = qry[:, hv].X; Xq = Xq.toarray() if hasattr(Xq, "toarray") else np.asarray(Xq)
    hv_pos = {g: i for i, g in enumerate(hv)}

    ref_targets = [t for t in targets if True]  # all 33 labeled targets
    nt_ref = (ref.obs["target_gene"] == "NT").values
    nt_qry = qry.obs["guide"].str.startswith("NT").values
    qguides = sorted([g for g in qry.obs["guide"].unique() if not g.startswith("NT")],
                     key=lambda s: int(s.replace("guide", "")))

    # ---- reference effects (one guide per target) ----
    eff_ref = {}
    for t in ref_targets:
        gname = ref.obs.loc[ref.obs["target_gene"] == t, "guide"].unique()
        assert len(gname) == 1
        eff_ref[t] = ols_effects(Xr, ref.obs["guide"].values, clu_ref, gname[0], nt_ref, hv_pos)

    # ---- query effects ----
    qlabels = qry.obs["guide"].values
    eff_q = {g: ols_effects(Xq, qlabels, clu_qry, g, nt_qry, hv_pos) for g in qguides}

    # ---- similarity with target-gene exclusion ----
    def sim_matrix(eff_q_g):
        scores = {}
        for t in ref_targets:
            excl = hv_pos.get(t)
            scores[t] = pearson_excl(eff_q_g, eff_ref[t], excl)
        return scores

    rows = []
    boot_detail = {}
    for g in qguides:
        sc_scores = sim_matrix(eff_q[g])
        ranked = sorted(sc_scores.items(), key=lambda kv: -kv[1])
        best_t, best_s = ranked[0]
        run_t, run_s = ranked[1]
        # bootstrap agreement
        cells_g = np.where(qlabels == g)[0]
        nt_idx = np.where(nt_qry)[0]
        votes = 0
        for _ in range(N_BOOT):
            sb = rng.choice(cells_g, size=len(cells_g), replace=True)
            nb = rng.choice(nt_idx, size=len(nt_idx), replace=True)
            keep = np.concatenate([sb, nb])
            lab_b = np.concatenate([np.full(len(sb), g), np.full(len(nb), "NT")])
            clu_b = clu_qry[keep]
            C = np.zeros((len(clu_b), clu_qry.max() + 1)); C[np.arange(len(clu_b)), clu_b] = 1
            D = np.column_stack([np.ones(len(clu_b)), C[:, 1:], (lab_b == g).astype(float)])
            beta, *_ = np.linalg.lstsq(D, Xq[keep], rcond=None)
            eb = beta[-1]
            b_scores = {t: pearson_excl(eb, eff_ref[t], hv_pos.get(t)) for t in ref_targets}
            votes += max(b_scores, key=b_scores.get) == best_t
        conf = votes / N_BOOT
        rows.append(dict(query_guide_id=g, target_gene=best_t, score=round(best_s, 4),
                         runner_up_gene=run_t, runner_up_score=round(run_s, 4),
                         confidence=round(conf, 3), n_cells=int(len(cells_g)),
                         gap=round(best_s - run_s, 4)))
        boot_detail[g] = dict(top3=[(t, round(s, 4)) for t, s in ranked[:3]])

    mapping = pd.DataFrame(rows)[["target_gene", "query_guide_id", "score",
                                  "runner_up_score", "confidence", "runner_up_gene", "gap", "n_cells"]]
    mapping.to_csv("output/guide_mapping.csv", index=False)

    # ---- validation 1: reference split-half self-mapping ----
    val = {}
    correct = 0; tot = 0
    for t in ref_targets:
        gname = ref.obs.loc[ref.obs["target_gene"] == t, "guide"].unique()[0]
        idx = np.where(ref.obs["guide"] == gname)[0]
        half = rng.permutation(len(idx)) < len(idx) // 2
        A, B = idx[half], idx[~half]
        def eff_from(sel_idx):
            keep = np.concatenate([sel_idx, np.where(nt_ref)[0]])
            lab = np.concatenate([np.full(len(sel_idx), gname),
                                  np.full(int(nt_ref.sum()), "NT")])
            C = np.zeros((len(keep), clu_ref.max() + 1)); C[np.arange(len(keep)), clu_ref[keep]] = 1
            D = np.column_stack([np.ones(len(keep)), C[:, 1:], (lab == gname).astype(float)])
            beta, *_ = np.linalg.lstsq(D, Xr[keep], rcond=None)
            return beta[-1]
        eA, eB = eff_from(A), eff_from(B)
        sims = {}
        for t2 in ref_targets:
            ref_e = eB if t2 == t else eff_ref[t2]
            sims[t2] = pearson_excl(eA, ref_e, hv_pos.get(t2))
        pred = max(sims, key=sims.get)
        correct += pred == t; tot += 1
        val[t] = dict(pred=pred, ok=bool(pred == t),
                      score=round(sims[t], 4),
                      top2=sorted(sims.items(), key=lambda kv: -kv[1])[:2])
    val_acc = correct / tot

    # ---- validation 2: NT negative controls ----
    nt_guides_q = sorted(qry.obs.loc[nt_qry, "guide"].unique())
    nt_null = []
    for g in nt_guides_q:
        e = ols_effects(Xq, qlabels, clu_qry, g, nt_qry & (qlabels != g), hv_pos)
        s = sim_matrix(e)
        nt_null.append(max(s.values()))
    nt_guides_r = sorted(ref.obs.loc[nt_ref, "guide"].unique())
    rlabs = ref.obs["guide"].values
    nt_null_r = []
    for g in nt_guides_r:
        e = ols_effects(Xr, rlabs, clu_ref, g, nt_ref & (rlabs != g), hv_pos)
        s = {t: pearson_excl(e, eff_ref[t], hv_pos.get(t)) for t in ref_targets}
        nt_null_r.append(max(s.values()))

    # ---- descriptive: cell-type shift ----
    comp = pd.crosstab(pd.Series(cluster_kmeans(np.vstack([emb_ref, emb_qry]))),
                       pd.Series(batch), normalize="columns")

    res = dict(
        shared_genes=len(shared), dropped_target_genes_from_query=dropped,
        targets_in_query_var=sorted(set(targets) & set(shared)),
        n_ref=int(ref.n_obs), n_qry=int(qry.n_obs),
        ref_self_map_top1_accuracy=round(val_acc, 4),
        ref_self_map_detail=val,
        nt_null_max_sim_query=[round(x, 4) for x in nt_null],
        nt_null_max_sim_ref=[round(x, 4) for x in nt_null_r],
        query_nt_guides=nt_guides_q,
        query_guide_mapping=rows,
        cell_type_shift_joint_cluster_composition=comp.round(3).to_dict(),
        targets_asked={t: next((r["query_guide_id"] for r in rows if r["target_gene"] == t), None)
                       for t in ["PABPC1", "NUDT21", "LEO1"]},
    )
    with open("output/results.json", "w") as f:
        json.dump(res, f, indent=1, default=str)

    print("=== query guide -> reference target ===")
    print(mapping.to_string(index=False))
    print("\nPABPC1/NUDT21/LEO1 query guides:", res["targets_asked"])
    print("ref split-half top-1 accuracy:", val_acc)
    print("NT null max sims (query):", np.round(nt_null, 3))
    print("NT null max sims (ref)  :", np.round(nt_null_r, 3))

if __name__ == "__main__":
    main()

