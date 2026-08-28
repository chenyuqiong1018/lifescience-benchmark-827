# Differential composition analysis — retinal single-cell samples 1 vs 2

**Date:** 2026-08-28 · **Task:** compare two retinal scRNA-seq count matrices and identify the cell population severely depleted in sample 2.

## 1. Result summary

| | sample1 | sample2 | ratio (s2/s1) |
|---|---|---|---|
| **Horizontal cell fraction** | **3.76%** (237 / 6,295) | **0.98%** (49 / 5,004) | **0.260** |

**Depleted population: `horizontal cell`.** Its sample-2/sample-1 fraction ratio (0.260) is the smallest among all cell types holding ≥1% fraction in sample 1 (runner-up: astrocyte, 0.274). Under composition preservation, sample 2 would be expected to contain ≈188 horizontal cells (5,004 × 3.76%); only 49 were observed (~3.8-fold fewer). Mean normalized expression of the horizontal-cell marker set across all cells also drops (0.088 → 0.063), corroborating the loss rather than a scoring artifact.

Deliverables: `output/composition.csv`, `output/depleted_call.json`, `output/analysis.py`, this report; `output/analysis_summary.json` holds full QC/scoring diagnostics. Reproduce with `python output/analysis.py` from the workspace root.

## 2. Inputs

| File | Content |
|---|---|
| `differential.composition.q1.1.mtx.gz` | sample 1, genes × cells MatrixMarket integer counts, 36,601 × **6,295**, 9,246,637 non-zeros |
| `differential.composition.q1.2.mtx.gz` | sample 2, genes × cells MatrixMarket integer counts, 36,601 × **5,004**, 9,953,348 non-zeros |
| `differential.composition.q1.genes.txt.gz` | 36,601 rows, `gene_ids`/`gene_symbols`, no duplicated symbols |
| `MARKER_PANEL.tsv` | 16 retinal cell types, 84 marker entries (82 unique) |
| `ANNOTATION_RULE.md` | frozen analysis rule (below) |

## 3. Methods (frozen rule, per `ANNOTATION_RULE.md`)

1. **Orientation:** each MTX matrix is transposed to cells × genes.
2. **Normalization:** per-cell library-size normalization × 10,000 (CP10K), then `log1p`. No gene filtering.
3. **Scoring:** per cell type, score = arithmetic mean of normalized expression over its marker genes.
4. **Assignment:** each cell → type with the largest score; ties follow marker-panel row order (`np.argmax` returns the first maximum). **No cells or genes are removed.**
5. **Fractions:** `n_cells / matrix column count` (all cells are assigned exactly one type, so fractions sum to 1).
6. **Depleted call:** type with the smallest sample-2/sample-1 fraction ratio among types with fraction ≥ 1% in sample 1.

All marker genes in the panel were found in the gene list; no imputation needed. `RLBP1` (RPE cell / muller glia cell) and `SLC1A3` (astrocyte / muller glia cell) are shared between panels by design.

## 4. QC (documentation only — frozen rule forbids filtering)

Integrity: both matrices contain non-negative integer counts; no empty cells (zero-library cells: 0 in both samples).

| Metric | sample1 | sample2 |
|---|---|---|
| Cells (matrix columns) | 6,295 | 5,004 |
| Library size: min / p5 / median / p95 / max | 205 / 830 / 1,999 / 9,250 / 45,431 | 204 / 1,837 / 3,197 / 9,954 / 48,129 |
| Genes detected: min / median / max | 181 / 1,166 / 8,802 | 174 / 1,672 / 9,052 |
| Mitochondrial % (MT- genes): median / p95 | 0.00 / 0.13 | 0.04 / 0.13 |
| Cells with library < 500 | 43 (0.7%) | 29 (0.6%) |
| Cells with genes detected < 200 | 2 | 5 |
| Cells with mito > 20% | 0 | 0 |

Interpretation: both samples are high-quality with low ambient mitochondrial signal; the small tails of low-content cells (≤0.7%) were **retained** as required by the frozen rule and are negligible relative to the observed effect. Sample 2 has a higher median library size (3,197 vs 1,999); CP10K normalization removes library-size scale before scoring, so this does not confound the composition comparison.

## 5. Composition results (`output/composition.csv`)

| Cell type | s1 n | s1 fraction | s2 n | s2 fraction | ratio s2/s1 |
|---|---:|---:|---:|---:|---:|
| rod cell | 3,976 | 63.16% | 3,558 | 71.10% | 1.126 |
| bipolar cell | 918 | 14.58% | 620 | 12.39% | 0.850 |
| muller glia cell | 512 | 8.13% | 433 | 8.65% | 1.064 |
| astrocyte | 285 | 4.53% | 62 | 1.24% | 0.274 |
| **horizontal cell** | **237** | **3.76%** | **49** | **0.98%** | **0.260** ← depleted |
| cone cell | 225 | 3.57% | 167 | 3.34% | 0.934 |
| RGC | 116 | 1.84% | 106 | 2.12% | 1.150 |
| microglial cell | 11 | 0.17% | 3 | 0.06% | 0.343* |
| T cells | 4 | 0.06% | 0 | 0.00% | —* |
| macrophage | 3 | 0.05% | 0 | 0.00% | —* |
| RPE cell | 2 | 0.03% | 3 | 0.06% | —* |
| Schwann cell | 2 | 0.03% | 2 | 0.04% | —* |
| Pericyte | 2 | 0.03% | 1 | 0.02% | —* |
| endothelial cell | 2 | 0.03% | 0 | 0.00% | —* |
| B cell | 0 | 0.00% | 0 | 0.00% | n/a |
| fibroblast | 0 | 0.00% | 0 | 0.00% | n/a |

\* types below the 1%-in-sample-1 threshold are excluded from the depletion criterion by the frozen rule (their ratios are unstable at ≤11 cells and are reported in `depleted_call.json` for reference). The seven candidate types are listed in `depleted_call.json`.

Retinal neurons dominate both samples (rods, bipolar, cones, horizontal, RGC, Müller glia), consistent with retinal tissue; rare vascular/immune/stromal types are near-noise-level in both samples.

## 6. Annotation evidence

- Assignment is deterministic under the frozen marker-mean rule; every cell receives exactly one label.
- Mean winning marker score per assigned type (log1p-CP10K units), sample1 → sample2: rod 1.66 → 2.88; muller glia 1.48 → 2.11; bipolar 1.23 → 1.45; horizontal 1.03 → 0.73; cone 0.75 → 1.27; RGC 0.98 → 1.27; astrocyte 0.37 → 0.45. Neuronal types show strong, consistent marker signal in both samples.
- Assignment confidence margin (winning score − runner-up): median 0.99 (s1) / 2.18 (s2) log units; only 6.9% (s1) and 4.3% (s2) of cells have margin < 0.1, i.e., >93% of cells are unambiguously assigned under this rule.
- The horizontal-cell call is supported at the population level too: mean expression of its 5 markers (LHX1, ONECUT1, ONECUT2, PROX1, CALB1) over *all* cells falls from 0.088 to 0.063 (independent dense-matrix recomputation).

## 7. Uncertainty and limitations

- **Near-tie at the top of the depletion ranking:** astrocyte (ratio 0.274) is close to horizontal cell (0.260). The frozen rule makes the call deterministic, but both populations drop ~3.5–4-fold in relative abundance; astrocyte loss is a real secondary signal worth reporting alongside the primary call.
- **Marker-only annotation:** scoring uses 5–6 markers per type and no transcriptome-wide clustering/differential expression; it is the benchmark's frozen reproducibility rule, not a full biological annotation workflow (per `ANNOTATION_RULE.md`).
- **Shared markers:** RLBP1 and SLC1A3 each contribute to two panels; a cell high in a shared marker gains score in both types (ties resolved by panel row order).
- **No cell filtering:** low-quality tails (≤43 cells/sample with library < 500) were retained per the frozen rule; they are far too few to alter any call.
- **Small-count types:** T cells, macrophage, endothelial (≤4 cells in s1) give unstable ratios and were excluded from the criterion by design.
- **No statistical test:** the depletion call is a point estimate of fraction ratios; no confidence intervals were requested by the frozen rule. The ~3.8-fold deficit of horizontal cells (49 observed vs ≈188 expected) is far beyond sampling noise for n ≈ 5,000 cells.
- **Numerical precision:** the production pipeline uses float32; an independent float64/matrix-multiplication recomputation reproduced identical per-type counts and the same depleted call (ratio 0.2601).

## 8. Verification performed

Independent re-implementation (different library-size computation via CSR `indptr` reduce-at, marker scoring via sparse matrix multiplication against a gene×type weight matrix, expected cell totals parsed from MTX headers) confirmed: per-type assignment counts for both samples, fractions summing to 1, totals matching the matrix headers (6,295 / 5,004), and the depleted call (`horizontal cell`, ratio ≈ 0.260).
