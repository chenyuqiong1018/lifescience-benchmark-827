# Differential composition analysis — retinal single-cell samples 1 vs 2

**Date:** 2026-08-28
**Question:** Which cell population is severely depleted in sample 2 relative to sample 1?

## Result

**`horizontal cell` is the population severely depleted in sample 2.**

| metric | sample 1 | sample 2 |
|---|---|---|
| horizontal cells (n) | 237 / 6,295 | 49 / 5,004 |
| fraction | 3.765% | 0.979% |
| sample-2 / sample-1 fraction ratio | — | **0.2601** |

Per the frozen rule, the depleted call is the listed cell type with the smallest
sample-2/sample-1 fraction ratio among types with >= 1% fraction in sample 1.
Horizontal cells have the smallest ratio (0.2601). Astrocytes are a close second
(ratio 0.2737) and are also strongly reduced; see Uncertainty below.

Full ratios for the 7 eligible types (>= 1% in sample 1):

| cell type | frac s1 | frac s2 | ratio s2/s1 |
|---|---|---|---|
| **horizontal cell** | **0.0376** | **0.0098** | **0.2601** |
| astrocyte | 0.0453 | 0.0124 | 0.2737 |
| bipolar cell | 0.1458 | 0.1239 | 0.8496 |
| cone cell | 0.0357 | 0.0334 | 0.9337 |
| muller glia cell | 0.0813 | 0.0865 | 1.0639 |
| rod cell | 0.6316 | 0.7110 | 1.1257 |
| RGC | 0.0184 | 0.0212 | 1.1495 |

Types below the 1% eligibility threshold in sample 1 (T cells, macrophage,
endothelial cell also drop to 0 cells in sample 2, but are excluded by the rule):
microglial cell, RPE cell, T cells, macrophage, Schwann cell, Pericyte, B cell,
fibroblast, endothelial cell.

## Inputs

| file | content |
|---|---|
| `inputs/differential.composition.q1.1.mtx.gz` | sample 1 raw UMI counts, 36,601 genes x 6,295 cells, 9,246,637 nonzeros |
| `inputs/differential.composition.q1.2.mtx.gz` | sample 2 raw UMI counts, 36,601 genes x 5,004 cells, 9,953,348 nonzeros |
| `inputs/differential.composition.q1.genes.txt.gz` | shared gene table (`gene_ids`, `gene_symbols`; no duplicate symbols) |
| `inputs/MARKER_PANEL.tsv` | 16 retinal cell types with marker gene lists (82 unique marker genes) |
| `inputs/ANNOTATION_RULE.md` | frozen annotation + normalization + depletion-call rule |

All 82 marker gene symbols were found in the gene table with a unique row each, so
no marker was dropped and no ambiguity arose in gene-to-row mapping.

## QC (descriptive only — no filtering)

The frozen rule explicitly forbids removing cells or genes, so QC is documented but
not used to filter.

| metric | sample 1 | sample 2 |
|---|---|---|
| cells | 6,295 | 5,004 |
| genes | 36,601 | 36,601 |
| nonzero entries | 9,246,637 | 9,953,348 |
| sparsity | 0.960 | 0.946 |
| library size per cell (min / median / mean / max) | 205 / 1,999 / 3,058 / 45,431 | 204 / 3,197 / 4,236 / 48,129 |
| genes detected per cell (min / median / max) | 181 / 1,166 / 8,802 | 174 / 1,672 / 9,052 |
| mitochondrial fraction per cell (median / max, 13 MT- genes) | 0.0000 / 0.0114 | 0.0004 / 0.0147 |
| empty droplets (zero counts) | 0 | 0 |

Interpretation:
- No empty droplets; every cell has >= ~200 UMIs, so the divide-by-library-size step
  is well defined (a guard would substitute library size 1 for any zero, not needed).
- Mitochondrial content is very low in both samples (max ~1.5%), indicating little
  apoptotic/stressed-cell signal by that criterion.
- Sample 2 cells were sequenced somewhat deeper (median library 3,197 vs 1,999).
  CP10k normalization mitigates this, but residual depth-related composition effects
  cannot be fully excluded (see Uncertainty).

## Normalization

Frozen rule, applied independently per matrix:
1. Transpose to cells x genes.
2. Divide each cell by its library size (total UMIs of the cell).
3. Multiply by 10,000 (CP10k).
4. `log1p`.

Implemented on the sparse matrix via a right diagonal scaling (per-cell divisors),
then `log1p` on the sparse matrix (`log1p(0)=0`, so sparsity is preserved).
No genes or cells removed.

## Annotation

Frozen rule: for each of the 16 panel cell types, score every cell by the arithmetic
mean of its log1p-CP10k expression over the type's marker genes; assign the cell to
the highest-scoring type; ties resolved by marker-panel row order. Implemented as a
single vectorized product `scores = X_norm @ M` with `M[g, t] = 1/(#markers of t)` if
gene g is a marker of type t, then `np.argmax` (first maximum => panel-order tie rule).
Denominator for all fractions is the matrix column count (total cells).

### Annotation evidence

Per-type evidence is saved in `output/annotation_evidence.json`
(mean marker score of assigned cells, mean margin over the runner-up type, and the
fraction of assigned cells with a zero winner-runner-up margin).

Highlights:
- Major neuronal populations are confidently assigned: rod cells (mean marker score
  1.66 in s1 / 2.88 in s2, margin 1.22 / 2.25), bipolar cells, cone cells, Muller
  glia and RGCs all show large margins (> 0.36) and near-zero tie fractions.
- Horizontal cells are assigned with good confidence: mean margin 0.73 (s1) and
  0.37 (s2); only 6.8% (s1) and 2.0% (s2) of assigned horizontal cells are exact ties.
- Globally, 2.73% (s1) and 0.52% (s2) of all cells are exact score ties (mostly cells
  expressing none of the relevant markers, defaulting to the earliest panel row,
  astrocyte). Astrocyte assignments carry the highest tie fraction (42.5% s1 /
  24.2% s2), consistent with astrocyte/Muller-glia marker overlap (e.g. SLC1A3) and
  low marker expression; the astrocyte fraction should therefore be read with caution.
- Rare types (<= 11 cells) have unstable assignments by construction and were
  excluded from the depletion call by the >= 1% eligibility rule.

## Uncertainty and limitations

1. **The depletion call is close.** Horizontal cells (ratio 0.2601) beat astrocytes
   (ratio 0.2737) by ~5% relative. Both populations are severely reduced in sample 2
   (roughly 4-fold and 3.6-fold respectively). The frozen rule selects a single
   answer deterministically; biologically, sample 2 shows concurrent loss of
   horizontal cells and astrocytes.
2. **Marker-score annotation is coarse.** The rule is frozen for reproducibility and
   is not a full biological annotation workflow: no background/negative-marker
   penalization, no clustering, no reference mapping. Overlapping markers
   (SLC1A3 in astrocyte and Muller glia; RLBP1 in RPE and Muller glia) can split or
   merge related types. Cells with no detectable marker expression default to the
   earliest panel row (astrocyte) via the tie rule.
3. **Depth difference between samples** (median library 3,197 vs 1,999) is handled by
   CP10k but can still shift marker-mean assignments for lowly expressed markers.
4. **Small counts for the called type in sample 2** (49 cells) make the sample-2
   fraction estimate noisy (binomial SE ~ sqrt(0.0098*0.9902/5004) ~ 0.0014, i.e.
   ~14% relative); even at the upper end the fraction stays far below sample 1's
   3.76%, so the depletion conclusion is robust. The astrocyte-vs-horizontal-cell
   ordering, however, sits inside the combined estimation noise of the two ratios.
5. No batch correction or differential-abundance significance testing was performed;
   the deliverable is a descriptive composition comparison under the frozen rule.

## Verification

- Production run: `output/analysis.py` (one run; deterministic).
- Independent verification run (`verify.py`, workspace root) re-derived normalization
  and assignments through a different code path (COO entry-wise normalization with
  `np.bincount` library sizes and `np.add.at` score accumulation instead of sparse
  diagonal algebra). All 32 fractions matched (within CSV rounding, 1e-11), and the
  depleted call reproduced exactly: `horizontal cell`, ratio 0.2601.
- Fractions sum to 1.0 within each sample; every cell is assigned exactly once.

## Artifacts

| path | description |
|---|---|
| `output/composition.csv` | sample, cell_type, n_cells, fraction (all 16 types x 2 samples) |
| `output/depleted_call.json` | depleted call + rule, fractions, ratios, eligible types |
| `output/analysis.py` | reproducible analysis script (frozen rule implementation) |
| `output/report.md` | this report |
| `output/qc.json` | per-sample QC metrics (descriptive) |
| `output/annotation_evidence.json` | per-type marker-score evidence and tie/uncertainty statistics |
