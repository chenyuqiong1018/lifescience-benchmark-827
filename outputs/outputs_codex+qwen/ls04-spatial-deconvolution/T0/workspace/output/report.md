# Deconvolution report: Spot_710-1

## Question
Identify the cell type or mixture captured at Visium spot `Spot_710-1`
(array_row 20, array_col 23, in_tissue, total 3,693 UMIs) using the supplied
single-cell reference.

## Data
- Visium spot-by-gene matrix: 900 spots x 1,000 genes
  (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`, `tissue_positions.csv`).
- Single-cell reference: 1,200 cells x 1,000 genes
  (`spatial_q_sc_counts.csv`) with 6 cell types, 200 cells each:
  Tumor_Core, Fibroblast_Stroma, T_Cell, B_Cell, Macrophage, Endothelial.

## Method (primary)
Per-cell-type signatures were built as mean CPM-normalized profiles (sum of
counts over the 200 reference cells per type, scaled to 10,000). The spot
count vector was CPM-normalized and decomposed by non-negative least squares
(NNLS, `scipy.optimize.nnls`): spot ~ signature x weights. Weights were
normalized to sum to 1.

## Result: a three-way mixture
| cell_type   | weight |
|-------------|--------|
| Endothelial | 0.359  |
| B_Cell      | 0.321  |
| Macrophage  | 0.320  |

NNLS fit R^2 = 0.61 on the full profile. T_Cell (0.012), Fibroblast_Stroma
(0.000) and Tumor_Core (0.000) received negligible/zero mass and were treated
as background, so the three retained components were renormalized to sum to 1.

## Evidence
1. **Marker-gene blocks.** Each cell type has a block of high fold-change
   marker genes. In Spot_710-1 all three blocks are strongly expressed:
   - Endothelial (e.g. Gene_278=12, Gene_280=14, Gene_283=10, Gene_297=11):
     20/20 top markers expressed, mean count 7.8.
   - B_Cell (e.g. Gene_156=10, Gene_167=7, Gene_190=7, Gene_194=7):
     20/20 top markers expressed, mean count 6.1.
   - Macrophage (e.g. Gene_202=8, Gene_222=8, Gene_224=8, Gene_233=8):
     20/20 top markers expressed, mean count 6.0.
   The other three types' blocks sit at background (14-17/20 expressed,
   mean counts 1.3-1.6).
2. **Marker-mass cell-count equivalents** (spot marker mass / per-cell
   reference marker mass over top-20 markers): Endothelial 1.27, Macrophage
   1.20, B_Cell 1.04 vs ~0.25-0.28 for the three absent types. Renormalized
   over the three components this gives 0.36 / 0.34 / 0.30 - consistent with
   the NNLS weights within ~0.05.
3. **Diagnostic - cosine nearest-neighbour voting.** All 30 nearest single
   cells are Endothelial (cos ~ 0.65). This is misleading on its own:
   ~99% of spot UMIs come from a shared expression program common to all
   cells in this simulation, so full-transcript cosine similarity is
   dominated by the shared component, while the sparse marker blocks carry
   the cell-type signal. The marker-based analyses (1-2) resolve the
   mixture, which is why the spot is not called pure Endothelial.

## Conclusion
Spot_710-1 is a **mixture of Endothelial (~36%), B_Cell (~32%) and
Macrophage (~32%)**, approximately three equal cell populations. A single
cell-type call is not supported: three distinct cell-type marker programs are
co-expressed at comparable strength in the spot.

## Caveats
- Weights are transcriptome-fraction estimates from NNLS on mean signatures;
  the two independent estimates agree within ~0.05 per component.
- The shared-expression structure of this dataset limits the resolution of
  similarity-based (nearest-neighbour) assignment; marker-gene-based
  deconvolution was used as the arbiter.

## Reproducibility
Run `python output/analysis.py` from the workspace root. It extracts
`inputs/spatial.sim.tar.gz` if needed, performs the NNLS deconvolution and
validation, writes `output/spot_710_composition.csv` (columns:
cell_type, weight, evidence; weights nonnegative, sum = 1.000), and asserts
the schema.
