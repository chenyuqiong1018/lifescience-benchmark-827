All deliverables are complete and validated.

## Summary

**Spot_710-1 is a three-way mixture: Endothelial 0.359, B_Cell 0.321, Macrophage 0.320** (weights sum to 1.000).

**Method (primary):** Cell-type centroids (mean CPM profile per annotated type from the 1,200-cell sc reference) fit to the spot's CPM profile via non-negative least squares (`scipy.optimize.nnls`), weights renormalized to sum to 1.

**Validation (independent checks, all concordant):**
- Exhaustive best-subset search: {B_Cell, Macrophage, Endothelial} is the best 1–3-type combination (R² = 0.610, weights 0.321/0.321/0.358 — matching the full model); runner-up subset only R² = 0.351 and every single-type-only fit has negative R² → a pure cell type is ruled out, and the mixture is well supported.
- Centroid Pearson r positive only for Endothelial (0.465), B_Cell (0.409), Macrophage (0.395).
- kNN voting (top-20 cells all Endothelial) is consistent with Endothelial as the largest component; the low single-cell correlation ceiling (~0.42) confirms no single type explains the spot.
- Spatial context: the 8 neighboring spots are Fibroblast/T-cell dominated, placing this spot at the edge of a distinct immune/vascular domain (recorded as context, not used in the fit).

**Deliverables (workspace `output/`):**
- `spot_710_composition.csv` — columns `cell_type,weight,evidence`; 3 rows; weights nonnegative, sum = 1.000 (schema asserted).
- `analysis.py` — full reproducible pipeline (also emits `deconv_stats.json`).
- `report.md` — methods, validation, result table, limitations.