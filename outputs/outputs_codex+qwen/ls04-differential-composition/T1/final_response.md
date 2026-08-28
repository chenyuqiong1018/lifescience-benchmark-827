Done. The cell population severely depleted in sample 2 is the **horizontal cell**.

## Key result
| | sample 1 | sample 2 | ratio |
|---|---|---|---|
| Horizontal cell | 237 / 6,295 = **3.76%** | 49 / 5,004 = **0.98%** | **0.260** |

Per the frozen rule in `ANNOTATION_RULE.md`, this is the smallest sample-2/sample-1 fraction ratio among types with ≥1% in sample 1 (runner-up: astrocyte, 0.274). Sample 2 would be expected to hold ≈188 horizontal cells if composition were preserved — only 49 were observed (~3.8-fold deficit).

## Deliverables (all in `output/`)
- **composition.csv** — `sample,cell_type,n_cells,fraction`, 32 rows (2 samples × 16 panel types); fractions sum to 1 per sample using matrix column count as denominator
- **depleted_call.json** — depleted type, fractions, counts, ratio, criterion, candidate list, and full ratio table for reference
- **analysis.py** — reproducible implementation of the frozen rule (transpose → CP10K → log1p → marker-mean scoring → argmax with panel-row tie-break; no cells/genes removed)
- **report.md** — documents QC, normalization, annotation evidence, uncertainty
- *analysis_summary.json* — full QC/scoring diagnostics backing the report

## QC & validation highlights
- **Integrity:** both matrices are non-negative integer counts (9.25M / 9.95M nnz), no zero-library cells; all 84 panel markers present, no duplicate gene symbols
- **QC (documentation-only, per frozen rule):** median library size 1,999 (s1) vs 3,197 (s2); mitochondrial % median ≤0.04%; ≤0.7% low-content cells per sample, retained as required
- **Annotation evidence:** >93% of cells have assignment margin ≥0.1 log units; neuronal types show strong consistent marker scores
- **Independent verification (different code path, float64):** reproduced identical per-type counts, totals matching MTX headers, same depleted call; horizontal-cell marker expression across all cells drops 0.088 → 0.063, corroborating true depletion rather than a scoring artifact
- **Uncertainty documented:** astrocyte is a near-tie secondary loss (0.274), shared markers RLBP1/SLC1A3, marker-only (not cluster-based) annotation, no filtering, small-count types excluded by design