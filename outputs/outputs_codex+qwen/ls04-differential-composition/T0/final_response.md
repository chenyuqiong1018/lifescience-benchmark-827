## Answer: **horizontal cell** is severely depleted in sample 2

| | sample 1 | sample 2 |
|---|---|---|
| horizontal cells | 237 / 6,295 (3.76%) | 49 / 5,004 (0.98%) |
| s2/s1 fraction ratio | — | **0.260** |

Per the frozen rule in `ANNOTATION_RULE.md`, the call is the type with the smallest sample-2/sample-1 fraction ratio among types ≥1% in sample 1. Horizontal cells win (0.2601), with astrocytes a close second (0.2737, also strongly reduced); rods/cones/bipolar/RGC/Müller glia ratios are all ≥0.85.

**Method** — raw MTX counts → per-cell CP10k + log1p → each cell scored by the arithmetic mean of each type's marker genes (`X_norm @ M` marker matrix, `argmax` with panel-order tie rule) → fractions over total cells. No cells/genes removed, per the frozen rule. All 82 marker genes mapped uniquely.

**QC (documented, no filtering)** — sample 1: 6,295 cells, median library 1,999; sample 2: 5,004 cells, median library 3,197; no empty droplets; mito fraction ≤1.5%. Sample 2 is sequenced ~1.6× deeper (depth effect mitigated by CP10k, noted as a limitation).

**Uncertainty** — the horizontal-vs-astrocyte ordering sits ~5% apart in ratio space (astrocytes also carry the highest tie fraction, 42%/24%, due to shared markers like SLC1A3); horizontal-cell assignments themselves are confident (mean margin 0.73/0.37). Binomial noise on the 49-cell sample-2 count cannot overturn the depletion conclusion.

**Verification** — an independent re-implementation (`verify.py`, COO entry-wise normalization + `np.add.at` scoring) reproduced all 32 fractions within rounding and the exact depleted call.

**Deliverables** (all in `output/`): `composition.csv` (sample,cell_type,n_cells,fraction; sums to 1.0 per sample), `depleted_call.json`, `analysis.py`, `report.md`, plus supporting `qc.json` and `annotation_evidence.json`.