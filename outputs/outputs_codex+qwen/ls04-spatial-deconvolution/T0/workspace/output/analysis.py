"""Deconvolution of Visium spot Spot_710-1 using a single-cell reference.

Primary method
    Non-negative least squares (NNLS) of the spot's count profile against
    per-cell-type mean signatures (CPM-normalized) built from the scRNA-seq
    reference.

Independent validation
    1. Marker-gene mass ratio: for each cell type, sum spot counts over its
       top-20 fold-change marker genes and divide by the per-cell reference
       marker mass -> cell-count equivalents per type.
    2. Marker-block inspection: expression of each cell type's differential
       gene block in the spot.
    3. Cosine nearest-neighbour voting against single cells (reported as a
        diagnostic; dominated by the shared expression program in this data).

Outputs
    output/spot_710_composition.csv with columns cell_type, weight, evidence.
    Weights are nonnegative and sum to 1 within 0.01.
"""
import os
import tarfile

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.optimize import nnls

HERE = os.path.dirname(os.path.abspath(__file__))          # .../workspace/output
WS = os.path.dirname(HERE)                                  # .../workspace
TAR = os.path.join(WS, "inputs", "spatial.sim.tar.gz")
EXT = os.path.join(WS, "extracted")
SPOT = "Spot_710-1"


def extract_inputs():
    needed = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz",
              "spatial_q_sc_counts.csv", "spatial_q_sc_metadata.csv",
              "tissue_positions.csv"]
    if not all(os.path.exists(os.path.join(EXT, f)) for f in needed):
        os.makedirs(EXT, exist_ok=True)
        with tarfile.open(TAR, "r:gz") as t:
            t.extractall(EXT)


def load_data():
    sc = pd.read_csv(os.path.join(EXT, "spatial_q_sc_counts.csv"), index_col=0)
    meta = pd.read_csv(os.path.join(EXT, "spatial_q_sc_metadata.csv"), index_col=0)
    feat = pd.read_csv(os.path.join(EXT, "features.tsv.gz"), sep="\t", header=None)
    barc = pd.read_csv(os.path.join(EXT, "barcodes.tsv.gz"), header=None)[0]
    mat = pd.DataFrame(
        sio.mmread(os.path.join(EXT, "matrix.mtx.gz")).toarray(),
        index=feat[1].values, columns=barc.values)
    genes = [g for g in mat.index if g in sc.columns]
    return sc[genes], meta, mat.loc[genes]


def main():
    extract_inputs()
    sc, meta, mat = load_data()
    spot = mat[SPOT].astype(float)
    print(f"Spot {SPOT}: total UMI = {spot.sum():.0f}")

    # ---- primary method: NNLS on CPM-normalized cell-type signatures ----
    sig = {}
    for ct, idx in meta.groupby("cell_type").groups.items():
        prof = sc.loc[idx].sum(axis=0).astype(float)
        sig[ct] = prof / prof.sum() * 1e4
    S = pd.DataFrame(sig)                                   # genes x cell types
    s_norm = spot / spot.sum() * 1e4
    w, _ = nnls(S.values, s_norm.values)
    w = pd.Series(w, index=S.columns)
    w_norm = w / w.sum()
    pred = S.values @ w.values
    ss = np.sum((s_norm.values - s_norm.values.mean()) ** 2)
    r2 = 1.0 - np.sum((s_norm.values - pred) ** 2) / ss
    print("NNLS weights (full profile):")
    print(w_norm.sort_values(ascending=False).round(4))
    print(f"NNLS fit R2 = {r2:.3f}")

    # ---- validation 1: marker-gene mass -> cell-count equivalents ----
    means = sc.groupby(meta.cell_type).mean()               # cell types x genes
    fc = means.div(means.mean(axis=0) + 1e-3, axis=1)
    top20 = {ct: fc.loc[ct].sort_values(ascending=False).head(20).index
             for ct in means.index}
    mass = pd.Series({
        ct: float(spot.reindex(g).sum()) / float(means.loc[ct, g].sum())
        for ct, g in top20.items()})
    print("Marker-mass cell-count equivalents:")
    print(mass.sort_values(ascending=False).round(3))

    # ---- validation 2: expression of each cell type's marker block ----
    block_stat = {}
    for ct, g in top20.items():
        v = spot.reindex(g)
        block_stat[ct] = (int((v > 0).sum()), round(float(v.mean()), 2))
    print("Marker block (top-20 FC genes): (#expressed, mean spot count)")
    for ct, v in sorted(block_stat.items(), key=lambda kv: -kv[1][1]):
        print(f"  {ct}: {v}")

    # ---- validation 3 (diagnostic): cosine nearest-neighbour voting ----
    def unit(v):
        n = np.linalg.norm(v)
        return v / n if n > 0 else v
    sv = unit(spot.values)
    cos = sc.apply(lambda r: float(np.dot(unit(r.values.astype(float)), sv)), axis=1)
    top30 = cos.sort_values(ascending=False).head(30).index
    votes = meta.loc[top30, "cell_type"].value_counts(normalize=True)
    print("Cosine top-30 neighbour votes (shared-program dominated):")
    print(votes.round(3))

    # ---- decide composition ----
    # Keep components with NNLS weight >= 0.05 (background level in this
    # simulation is ~0.25 cell-count-equivalents / weight < 0.05).
    keep = w_norm[w_norm >= 0.05]
    keep = keep / keep.sum()
    print("Final composition:")
    print(keep.sort_values(ascending=False).round(4))

    # ---- evidence strings ----
    def evidence(ct):
        n_exp, mean_c = block_stat[ct]
        genes = ", ".join(list(top20[ct][:3]))
        return (
            f"NNLS_fit_R2={r2:.2f}; nnls_weight={w_norm[ct]:.3f}; "
            f"marker_mass_fraction={(mass / mass.sum())[ct]:.2f}; "
            f"top20_marker_genes_expressed={n_exp}/20 (mean count {mean_c}); "
            f"example markers {genes}"
        )

    out = pd.DataFrame({
        "cell_type": keep.sort_values(ascending=False).index,
        "weight": keep.sort_values(ascending=False).round(3).values,
        "evidence": [evidence(ct) for ct in keep.sort_values(ascending=False).index],
    })
    os.makedirs(os.path.join(WS, "output"), exist_ok=True)
    out_path = os.path.join(WS, "output", "spot_710_composition.csv")
    out.to_csv(out_path, index=False)
    print(f"wrote {out_path}")

    # ---- schema assertion ----
    chk = pd.read_csv(out_path)
    assert list(chk.columns) == ["cell_type", "weight", "evidence"], chk.columns
    assert (chk["weight"] >= 0).all()
    assert abs(chk["weight"].sum() - 1.0) <= 0.01, chk["weight"].sum()
    assert set(chk["cell_type"]) <= set(meta.cell_type.unique())
    print("schema assertion passed: columns ok, weights nonnegative, "
          f"sum={chk['weight'].sum():.3f}, rows={len(chk)}")


if __name__ == "__main__":
    main()
