# Differential composition analysis — retinal single-cell samples 1 vs 2

## Question
Compare two retinal single-cell RNA-count matrices and identify the cell population
severely depleted in sample 2.

## Inputs
| File | Content |
|---|---|
| `inputs/differential.composition.q1.1.mtx.gz` | sample 1, Matrix Market, 36,601 genes x 6,295 cells, 9,246,637 non-zeros |
| `inputs/differential.composition.q1.2.mtx.gz` | sample 2, Matrix Market, 36,601 genes x 5,004 cells, 9,953,348 non-zeros |
| `inputs/differential.composition.q1.genes.txt.gz` | shared feature table (`gene_ids`, `gene_symbols`), 36,601 rows, all symbols unique |
| `inputs/MARKER_PANEL.tsv` | 16 cell types with marker gene lists (81 marker genes) |
| `inputs/ANNOTATION_RULE.md` | frozen composition rule (see Normalization/Annotation) |

## QC
- Matrix orientation verified: MTX dimensions match the 36,601-row gene table; row i of each
  matrix corresponds to row i of `genes.txt.gz`.
- **All 81 panel marker genes are present** in the shared gene symbol list; no duplicates, so the
  symbol-to-row mapping is unambiguous. No cells or genes were removed (per the frozen rule).
- Per-sample QC metrics (no filtering applied):

| Metric | sample1 | sample2 |
|---|---|---|
| Cells (matrix columns) | 6,295 | 5,004 |
| Genes | 36,601 | 36,601 |
| Library size — mean / median | 3,058 / 1,999 | 4,236 / 3,197 |
| Library size — min / max | 205 / 45,431 | 204 / 48,129 |
| Genes detected per cell — median (mean) | 1,166 (1,469) | 1,672 (1,989) |

- Sample 2 has ~1.6x higher median library size than sample 1. The per-cell library-size
  normalization below absorbs this global depth difference, so composition comparisons are not
  driven by sequencing depth; residual gene-length/GC biases common to both samples largely cancel
  in a within-sample fraction.
- Cells with extreme library sizes (up to ~45–48k UMIs) are kept by the frozen rule; such cells can
  be doublets, which would slightly blur rare-type fractions. This is a known limitation, not
  corrected here by design.

## Normalization
Frozen rule applied exactly: transpose to cells-by-genes, divide each cell by its library size,
multiply by 10,000 (CP10K), then `log1p`. Implemented with sparse row scaling
(`diag(1e4/lib) @ X_cells_x_genes`) followed by in-place `log1p` on non-zero entries
(`log1p(0) = 0`, so sparsity is preserved).

## Cell-type annotation
Per the frozen rule: for each of the 16 panel cell types, each cell receives the arithmetic mean of
the type's marker genes; the cell is assigned to the type with the largest mean score, ties broken
by panel row order (`numpy.argmax` returns the first maximum, matching this convention). Fractions
use the matrix column count as denominator.

Annotation evidence:
- Assigned cells show clearly separated mean marker scores for the major retinal types, e.g.
  sample 1 median assigned score: rod 1.66, Muller glia 1.48, bipolar 1.23, horizontal 1.03,
  RGC 0.98, cone 0.75; immune/stromal types that receive only stray cells score < 0.45.
- Horizontal cell depletion is corroborated at the raw-marker level: cells with >=2 detected
  horizontal-marker UMIs drop from 282 (sample 1) to 119 (sample 2); LHX1-expressing cells drop
  from 122 to 4, ONECUT1 from 275 to 103, ONECUT2 from 272 to 73. (PROX1/CALB1 are less specific —
  CALB1 is also in photoreceptors/other neurons — but the specific markers LHX1/ONECUT1/ONECUT2
  all collapse in sample 2.)

## Composition results (`output/composition.csv`)

| cell_type | sample1 n | sample1 frac | sample2 n | sample2 frac | ratio s2/s1 |
|---|---|---|---|---|---|
| rod cell | 3,976 | 0.631612 | 3,558 | 0.711031 | 1.126 |
| bipolar cell | 918 | 0.145830 | 620 | 0.123901 | 0.850 |
| muller glia cell | 512 | 0.081334 | 433 | 0.086531 | 1.064 |
| astrocyte | 285 | 0.045274 | 62 | 0.012390 | 0.274 |
| **horizontal cell** | **237** | **0.037649** | **49** | **0.009792** | **0.260** |
| cone cell | 225 | 0.035743 | 167 | 0.033373 | 0.934 |
| RGC | 116 | 0.018427 | 106 | 0.021183 | 1.150 |
| microglial cell | 11 | 0.001747 | 3 | 0.000600 | n/a (<1% in s1) |
| T cells | 4 | 0.000635 | 0 | 0.0 | n/a |
| macrophage | 3 | 0.000477 | 0 | 0.0 | n/a |
| RPE cell | 2 | 0.000318 | 3 | 0.000600 | n/a |
| Schwann cell | 2 | 0.000318 | 2 | 0.000400 | n/a |
| Pericyte | 2 | 0.000318 | 1 | 0.000200 | n/a |
| endothelial cell | 2 | 0.000318 | 0 | 0.0 | n/a |
| B cell | 0 | 0.0 | 0 | 0.0 | n/a |
| fibroblast | 0 | 0.0 | 0 | 0.0 | n/a |

## Depleted call (`output/depleted_call.json`)
Per the frozen rule — among types with >=1% fraction in sample 1, the smallest
sample2/sample1 fraction ratio — the severely depleted population is:

> **horizontal cell**: 3.765% in sample 1 (237/6,295) -> 0.979% in sample 2 (49/5,004),
> ratio = 0.260 (~4-fold depletion, dropping below the 1% level).

Full eligible ranking by ratio: horizontal cell 0.260, astrocyte 0.274, bipolar cell 0.850,
cone cell 0.934, Muller glia 1.064, rod cell 1.126, RGC 1.150.

## Uncertainty and caveats
- **Near-tie at the top**: astrocyte (ratio 0.274) is nearly as depleted as horizontal cell
  (0.260); the separation is ~5% of the ratio value and rests on 237 vs 62 and 49 vs 62 cell
  counts, so sampling noise is non-trivial. The frozen rule selects horizontal cell; both
  populations are strongly reduced in sample 2.
- **Assignment ambiguity**: cells whose top-two type scores differ by <0.1 are 434/6,295 (6.9%)
  in sample 1 and 213/5,004 (4.3%) in sample 2 (score gap <0.01: 190 and 48 cells). Reassigning
  these borderline cells could shift small-type counts by a few dozen cells at most; the ~4-fold
  horizontal-cell depletion is larger than this noise band.
- **Marker overlaps**: RLBP1 is shared by RPE and Muller glia, SLC1A3 by astrocyte and Muller
  glia, PROX1/CALB1 are broad neuronal markers; argmax with panel-order tie-break is the frozen
  convention for these collisions.
- **No cell QC filtering** (frozen rule): doublets/ambient-RNA-contaminated cells remain and can
  inflate rare-type assignments; immune/stromal assignments at <0.2% fractions should be read as
  "present at trace levels" rather than precise counts.
- **Method scope**: marker-score argmax without clustering or differential-expression validation
  is a reproducible benchmark rule, not a full biological annotation workflow (as stated in
  `ANNOTATION_RULE.md`).
- The two samples differ in total cell counts (6,295 vs 5,004) and sequencing depth; using
  fractions (per the rule) and CP10K normalization makes the comparison compositionally valid,
  but absolute cell-number differences cannot be interpreted biologically from this data alone.

## Reproducibility
- `output/analysis.py` — full pipeline (QC, normalization, annotation, composition, depleted
  call); writes `composition.csv`, `depleted_call.json`, and `qc.json` (all QC metrics above).
  Deterministic; no random sampling.
- Verification: an independent dense (`numpy`) reimplementation of the frozen rule reproduced all
  32 per-sample cell counts exactly (check passed).
