"""Spatial deconvolution of Visium spot Spot_710-1 using a matched scRNA-seq reference.

Primary method: non-negative least squares (NNLS) fit of the spot's CPM-normalized
expression profile against cell-type centroid profiles built from the single-cell
reference (average CPM per annotated cell type). Weights are renormalized to sum to 1.

Validation (independent of the primary fit):
  1) Exhaustive best-subset search over 1-3 cell-type combinations (NNLS per subset).
  2) Pearson correlation of the spot profile against each cell-type centroid.
  3) k-nearest reference-cell voting (top-20 Pearson-correlated single cells).
  4) Spatial context: NNLS composition of the 8 neighboring spots.
"""
import gzip
import json
import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.optimize import nnls
from scipy.stats import pearsonr

WS = os.path.dirname(os.path.abspath(__file__))           # .../workspace/output
ROOT = os.path.dirname(WS)                                # .../workspace
IN = os.path.join(ROOT, "inputs")
SPOT = "Spot_710-1"

# ---------------------------------------------------------------- load sc reference
meta = pd.read_csv(os.path.join(IN, "spatial_q_sc_metadata.csv"), index_col=0)
sc = pd.read_csv(os.path.join(IN, "spatial_q_sc_counts.csv"), index_col=0)
assert sc.index.equals(meta.index), "sc counts/metadata cell mismatch"
genes = sc.columns.values

sc_cpm = sc.div(sc.sum(axis=1), axis=0) * 1e6             # CPM per cell
cell_types = list(meta["cell_type"].value_counts().index)
centroids = np.vstack([sc_cpm[meta["cell_type"] == ct].mean(axis=0).values
                       for ct in cell_types])              # types x genes

# ---------------------------------------------------------------- load Visium spot
mat = mmread(os.path.join(IN, "matrix.mtx.gz")).tocsc()    # genes x spots
with gzip.open(os.path.join(IN, "barcodes.tsv.gz"), "rt") as fh:
    barcodes = [l.strip() for l in fh]
feat = pd.read_csv(os.path.join(IN, "features.tsv.gz"), sep="\t", header=None)
assert list(feat[1].values) == list(genes), "gene order mismatch sc vs Visium"

idx = barcodes.index(SPOT)
x = np.asarray(mat[:, idx].todense()).ravel().astype(float)
x_cpm = x / x.sum() * 1e6

def r2(fit):
    return float(1.0 - ((x_cpm - fit) ** 2).sum() / ((x_cpm - x_cpm.mean()) ** 2).sum())

# ---------------------------------------------------------------- primary: NNLS
A = centroids.T                                            # genes x types
w_raw, _ = nnls(A, x_cpm)
w = w_raw / w_raw.sum()
r2_full = r2(A @ w_raw)

# ---------------------------------------------------------------- validation
# (1) best-subset search over 1-3 cell types
res = []
for n in [1, 2, 3]:
    for combo in combinations(range(len(cell_types)), n):
        Asub = centroids[list(combo)].T
        ww, _ = nnls(Asub, x_cpm)
        if ww.sum() <= 0:
            continue
        res.append((r2(Asub @ ww), [cell_types[i] for i in combo], ww / ww.sum()))
res.sort(key=lambda t: -t[0])
best_r2, best_combo, best_w = res[0]

# (2) centroid Pearson correlations
corr = {ct: float(pearsonr(x_cpm, centroids[i])[0]) for i, ct in enumerate(cell_types)}

# (3) top-k nearest reference cells by Pearson correlation (vectorized)
M = sc_cpm.values
mc = M - M.mean(axis=1, keepdims=True)
xc = x_cpm - x_cpm.mean()
per_cell = (mc @ xc) / (np.linalg.norm(mc, axis=1) * np.linalg.norm(xc) + 1e-12)
k = 20
votes = meta.iloc[np.argsort(per_cell)[::-1][:k]]["cell_type"].value_counts()

# (4) spatial context: NNLS composition of the 8 neighboring spots
pos = pd.read_csv(os.path.join(IN, "tissue_positions.csv"))
r = pos[pos["barcode"] == SPOT].iloc[0]
nb = pos[(pos["array_row"].between(r["array_row"] - 1, r["array_row"] + 1))
         & (pos["array_col"].between(r["array_col"] - 1, r["array_col"] + 1))
         & (pos["barcode"] != SPOT) & (pos["in_tissue"] == 1)]
nb_comp = {}
for bc in nb["barcode"]:
    v = np.asarray(mat[:, barcodes.index(bc)].todense()).ravel().astype(float)
    v = v / v.sum() * 1e6
    wj, _ = nnls(A, v)
    nb_comp[bc] = {ct: float(wj[i] / wj.sum()) for i, ct in enumerate(cell_types)}
nb_mean = {ct: float(np.mean([nb_comp[b][ct] for b in nb_comp])) for ct in cell_types}

# ---------------------------------------------------------------- write CSV
# Report supported components (NNLS weight >= 0.05); renormalize to sum to 1.
keep = [i for i in range(len(cell_types)) if w[i] >= 0.05]
wrep = w_raw[keep] / w_raw[keep].sum()
rows = []
for new_i, i in enumerate(keep):
    ct = cell_types[i]
    in_combo = "yes" if ct in best_combo else "no"
    ev = (f"NNLS weight {wrep[new_i]:.3f} (full-model fit R2={r2_full:.3f}); "
          f"member of best-fit cell-type subset (R2={best_r2:.3f}): {in_combo}; "
          f"Pearson r vs {ct} centroid {corr[ct]:.3f}; "
          f"{int(votes.get(ct, 0))}/{k} top-correlated reference cells are {ct}; "
          f"mean NNLS weight in 8 neighboring spots {nb_mean[ct]:.3f}")
    rows.append({"cell_type": ct, "weight": round(float(wrep[new_i]), 3), "evidence": ev})
out = pd.DataFrame(rows, columns=["cell_type", "weight", "evidence"])
out.to_csv(os.path.join(WS, "spot_710_composition.csv"), index=False)

# ---------------------------------------------------------------- summary + stats
print("Spot:", SPOT, "total UMIs:", int(x.sum()),
      "grid (array_row,array_col):", int(r["array_row"]), int(r["array_col"]))
print("Full NNLS weights:", {ct: round(float(w[i]), 4) for i, ct in enumerate(cell_types)})
print("Reported:", dict(zip(out.cell_type, out.weight)), "sum:", round(out.weight.sum(), 3))
print("Full-model R2:", round(r2_full, 4))
print("Best subset:", best_combo, "R2:", round(best_r2, 4), "w:", np.round(best_w, 3))
print("Runner-up:", res[1][1], "R2:", round(res[1][0], 4))
print("Centroid Pearson r:", {ct: round(v, 3) for ct, v in corr.items()})
print("Top-%d cell votes:" % k, votes.to_dict())
print("Neighbor mean composition:", {ct: round(v, 3) for ct, v in nb_mean.items()})

with open(os.path.join(WS, "deconv_stats.json"), "w") as fh:
    json.dump({
        "spot": SPOT, "umi": int(x.sum()),
        "grid": [int(r["array_row"]), int(r["array_col"])],
        "r2_full": r2_full,
        "weights_full": {ct: float(w[i]) for i, ct in enumerate(cell_types)},
        "best_combo": best_combo, "best_r2": best_r2,
        "best_w": {ct: float(vv) for ct, vv in zip(best_combo, best_w)},
        "runner_up": {"combo": res[1][1], "r2": res[1][0]},
        "corr": corr, "votes": votes.to_dict(), "k": k,
        "nb_mean": nb_mean, "reported": out.to_dict(orient="records"),
    }, fh, indent=1)
print("Wrote", os.path.join(WS, "spot_710_composition.csv"))
