# Deconvolution of Visium spot Spot_710-1

## Question
Identify the cell type or cell-type mixture captured at Visium spot `Spot_710-1`
(array_row 20, array_col 23, in_tissue=1, 3,693 UMIs), using the matched
single-cell reference (1,200 cells; 6 annotated types; 1,000 shared genes).

## Data
- `inputs/spatial_q_sc_counts.csv` + `spatial_q_sc_metadata.csv`: 1,200 reference
  cells (200 each of Tumor_Core, Fibroblast_Stroma, T_Cell, B_Cell, Macrophage,
  Endothelial), 1,000 genes.
- `inputs/matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` +
  `tissue_positions.csv`: 900 Visium spots x 1,000 genes (gene order identical to
  the reference).

## Method (primary)
Cell-type centroids were built as the mean CPM-normalized profile of all
reference cells of each type. The CPM-normalized spot profile was then fit as a
non-negative linear combination of the six centroids using scipy `nnls`; weights
were renormalized to sum to 1. Components with weight < 0.05 were dropped and the
remainder renormalized.

## Validation (independent checks)
1. **Best-subset search**: exhaustive NNLS over all 1-3 cell-type subsets.
   Best subset = {B_Cell, Macrophage, Endothelial}, R2 = 0.610 with weights
   0.321 / 0.321 / 0.358 - identical (within rounding) to the full-model
   weights. The runner-up subset (T_Cell, B_Cell, Endothelial) reaches only
   R2 = 0.351, and every single cell-type-only fit has negative R2, ruling out a
   pure cell type.
2. **Centroid correlations**: Pearson r of the spot vs centroids is positive only
   for Endothelial (0.465), B_Cell (0.409) and Macrophage (0.395); Tumor_Core,
   Fibroblast_Stroma and T_Cell are ~0 or negative.
3. **Nearest-cell voting**: all 20 most correlated reference cells are
   Endothelial (max single-cell r ~ 0.42). This is consistent with Endothelial
   being the largest single component, but the low correlation ceiling shows no
   single cell type explains the spot, matching the mixture conclusion.
4. **Spatial context**: the 8 neighboring spots are dominated by
   Fibroblast_Stroma (~0.50) and T_Cell (~0.31) with little B_Cell/Macrophage,
   so Spot_710-1 sits at the edge of a distinct immune/vascular domain. Neighbor
   composition was recorded as context evidence only and was not used to fit the
   target spot.

## Result
Spot_710-1 is a three-way mixture:

| cell_type  | weight |
|------------|--------|
| Endothelial| 0.359  |
| B_Cell     | 0.321  |
| Macrophage | 0.320  |

Weights sum to 1.000. Tumor_Core, Fibroblast_Stroma and T_Cell receive
negligible support (full-model NNLS weights 0.000, 0.000 and 0.012).

## Limitations
- Centroid-based NNLS assumes the spot is a linear mixture of mean cell-type
  profiles; cell-state heterogeneity within types is not modeled.
- Full-model fit R2 = 0.61, so part of the spot variance is unexplained
  (dropout/noise and within-type heterogeneity).
- kNN voting favors Endothelial because single-cell profiles are noisy; fit
  quality (R2) and subset search were treated as the decisive evidence.
- Neighboring spots differ strongly in composition, so spatial smoothing was
  intentionally not applied to the target spot.

Reproducible via `python output/analysis.py` (writes
`output/spot_710_composition.csv` and `output/deconv_stats.json`).
