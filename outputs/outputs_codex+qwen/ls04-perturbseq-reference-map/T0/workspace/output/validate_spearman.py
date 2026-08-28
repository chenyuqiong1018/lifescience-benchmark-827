"""Lightweight independent validation: rank-based (Spearman) similarity on the same
cell-type-adjusted effect vectors; compare argmax mapping to the primary (Pearson) result."""
import numpy as np, pandas as pd
from scipy.stats import rankdata
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import analysis as A

ref, qry, shared, targets, dropped = A.load_aligned()
hv, emb_r, emb_q, batch = A.hvg_pca(ref, qry, shared)
clu_r, clu_q = A.cluster_kmeans(emb_r), A.cluster_kmeans(emb_q)
hv_pos = {g: i for i, g in enumerate(hv)}
Xr = ref[:, hv].X; Xr = Xr.toarray() if hasattr(Xr, "toarray") else np.asarray(Xr)
Xq = qry[:, hv].X; Xq = Xq.toarray() if hasattr(Xq, "toarray") else np.asarray(Xq)
nt_ref = (ref.obs["target_gene"] == "NT").values
nt_qry = qry.obs["guide"].str.startswith("NT").values
qlabels = qry.obs["guide"].values
qguides = sorted([g for g in np.unique(qlabels) if not g.startswith("NT")],
                 key=lambda s: int(s.replace("guide", "")))
eff_ref = {}
for t in targets:
    gname = ref.obs.loc[ref.obs["target_gene"] == t, "guide"].unique()
    eff_ref[t] = A.ols_effects(Xr, ref.obs["guide"].values, clu_r, gname[0], nt_ref, hv_pos)
eff_q = {g: A.ols_effects(Xq, qlabels, clu_q, g, nt_qry, hv_pos) for g in qguides}

def spearman_excl(a, b, excl):
    keep = np.ones(len(a), bool)
    if excl is not None: keep[excl] = False
    ra, rb = rankdata(a)[keep], rankdata(b)[keep]
    ra -= ra.mean(); rb -= rb.mean()
    d = np.linalg.norm(ra) * np.linalg.norm(rb)
    return float(ra @ rb / d) if d > 0 else 0.0

res = {}
for g in qguides:
    s = {t: spearman_excl(eff_q[g], eff_ref[t], hv_pos.get(t)) for t in targets}
    rk = sorted(s.items(), key=lambda kv: -kv[1])
    res[g] = rk

prim = pd.read_csv("output/guide_mapping.csv")
prim_map = dict(zip(prim.query_guide_id, prim.target_gene))
agree = sum(res[g][0][0] == prim_map[g] for g in qguides)
print(f"Spearman vs Pearson top-1 agreement: {agree}/{len(qguides)}")
print("disagreements (guide: pearson->spearman):")
for g in qguides:
    if res[g][0][0] != prim_map[g]:
        print(f"  {g}: {prim_map[g]} -> {res[g][0][0]} (spearman top3: {res[g][:3]})")
for g in ["guide13", "guide14", "guide15", "guide18"]:
    print(f"{g}: pearson={prim_map[g]:6s} spearman_top3={[(t, round(s,3)) for t,s in res[g][:3]]}")
