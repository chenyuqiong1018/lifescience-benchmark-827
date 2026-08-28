# Spot_710-1 cell-type composition (Visium deconvolution)

## Conclusion

Spot_710-1 is a **three-way mixture** of **Endothelial (~0.35)**, **B_Cell (~0.32)**, **Macrophage (~0.32)**. No single cell type explains the spot; the mixture is strongly supported by residual-sum-of-squares model selection and reproduced by an independent cell-level deconvolution.

## Inputs

- Visium matrix: 1000 genes x 900 spots (target Spot_710-1: 3693 total UMI counts).
- scRNA-seq reference: 1200 cells x 1000 genes; 6 annotated cell types, 200 cells each.
- Shared genes used: 1000.

## Primary method

1. CPM (library-size) normalization of reference cells and the spot (a spot's count profile is a linear mixture of cell CPM profiles).
2. Reference profile per cell type = mean CPM across its 200 cells.
3. Nonnegative least squares `spot_cpm ~ R @ w`, weights renormalized to sum to 1.

## Estimated composition

| cell_type | weight | cosine vs profile | pearson vs profile |
|---|---|---|---|
| Endothelial | 0.3544 | 0.723 | 0.465 |
| B_Cell | 0.3174 | 0.693 | 0.409 |
| Macrophage | 0.3166 | 0.687 | 0.395 |
| T_Cell | 0.0116 | 0.487 | -0.013 |
| Fibroblast_Stroma | 0.0000 | 0.457 | -0.062 |
| Tumor_Core | 0.0000 | 0.466 | -0.058 |

Mixture reconstruction R2 = 0.610 (residual floor is gene-level Poisson sampling noise).

## Validation

**A. Cell-level NNLS** (spot against all 1200 individual cells, weights aggregated by type) reproduces the same mixture: Endothelial 0.35, B_Cell 0.31, Macrophage 0.30, T_Cell 0.02.

**B. Correlation ranking**: all 30 most correlated reference cells are Endothelial; among the top 100, Endothelial=98, B_Cell=2 — consistent with Endothelial as a major component of the mixture.

**C. Model selection (RSS, lower is better)**:

| model | RSS |
|---|---|
| best pure type (Endothelial) | 8.174e+08 |
| best pair (B_Cell + Endothelial) | 4.778e+08 |
| best triplet (B_Cell + Endothelial + Macrophage) | 2.783e+08 |
| full 6-type NNLS | 2.781e+08 |
| total SS (null) | 7.132e+08 |

The triplet improves RSS by ~39% over the best pair and ~66% over the best pure type, so a mixture (not a single type) is supported. Reference profiles are mutually orthogonal (pairwise |pearson| < 0.1), so the NNLS weights are not a collinearity artifact.

## Deliverables

- `output/spot_710_composition.csv` — cell_type, weight, evidence (weights nonnegative, sum to 1 within 0.01).
- `output/analysis.py` — this analysis.
- `output/report.md` — this report.
