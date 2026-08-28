Done. **Spot_710-1 is a three-way mixture of Endothelial (~0.35), B_Cell (~0.32), and Macrophage (~0.32)** — not a single cell type.

**Method** (one primary + one lightweight validation, as required):
- **Primary:** CPM normalization of the spot and the 1200-cell scRNA-seq reference → mean-CPM profiles per cell type → nonnegative least squares deconvolution with weights renormalized to sum to 1 (reconstruction R²=0.61; residual floor is per-gene sampling noise).
- **Validation:** independent cell-level NNLS reproduced the same mixture (Endo 0.35 / B 0.31 / Macro 0.30); correlation ranking of the spot against all 1200 cells puts Endothelial dominant among top hits; RSS model selection decisively favors the triplet (pure best 8.2e8 → best pair 4.8e8 → triplet 2.78e8), and reference profiles are mutually orthogonal (|r|<0.1), ruling out collinearity artifacts.

**Deliverables** (all schema-checked: exact columns `cell_type,weight,evidence`, 6 rows, weights ≥ 0, sum = 1.0000):

| cell_type | weight |
|---|---|
| Endothelial | 0.3544 |
| B_Cell | 0.3174 |
| Macrophage | 0.3166 |
| T_Cell | 0.0116 (noise) |
| Fibroblast_Stroma | 0.0 |
| Tumor_Core | 0.0 |

- `output/spot_710_composition.csv` — composition with per-type evidence strings
- `output/analysis.py` — full reproducible analysis (primary NNLS + all three validation checks + report generation)
- `output/report.md` — conclusion, inputs, method, estimated composition, and validation details