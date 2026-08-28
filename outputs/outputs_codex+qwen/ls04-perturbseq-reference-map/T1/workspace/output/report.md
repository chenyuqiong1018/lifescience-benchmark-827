# Perturb-seq query-to-reference guide mapping — report

## Task
Map the anonymized query perturbation groups (`guide1..guide32`, plus `NT-I*`
non-targeting controls) onto the labeled reference (43 guides / 34 target genes,
one guide per perturbed gene), robust to the cell-type/state shift between the two
datasets, and identify the query guide IDs corresponding to **PABPC1, NUDT21, LEO1**.

## Data
| | reference | query |
|---|---|---|
| cells | 4183 | 3516 |
| genes (measured) | 11,889 | 11,889 |
| shared genes used | 11858 | 11858 |
| obs metadata | `guide`, `target_gene` | `guide` only |
| guides | 10 NT + 33 perturbed (1/gene) | 10 NT + 32 anonymized perturbed |
| matrix | raw UMI counts | raw UMI counts |

Cross-dataset shift: the NT-control mean expression profiles of the two datasets
correlate at r = 0.918, i.e. the query lives in a shifted
cellular context. Matching is therefore performed on *perturbation effects*
(guide mean minus own-dataset NT control mean), never on absolute expression.

## Leakage prevention
1. **No target metadata in the query.** Asserted at runtime: the query AnnData
   `obs` contains only the anonymized `guide` column; no gene/target field exists.
2. **No name-based matching.** Guide-name strings are never compared across
   datasets and never parsed to infer gene identity. The only label semantics used
   is the explicit `NT` non-targeting-control designation, which is part of the
   query's own experimental design (needed to estimate the control baseline).
3. **Within-dataset baselines only.** Each dataset is normalized independently
   (CPM 1e4 + log1p) and signatures are computed relative to each dataset's own NT
   mean, so no cross-dataset baseline or cell-type composition information leaks
   into the scores.
4. **Label-column aggregation only.** Reference guide signatures are grouped to
   `target_gene` strictly via the `target_gene` label column.

## Primary method
1. Normalize each dataset independently (library-size 1e4, log1p).
2. Restrict to the 11858 shared genes.
3. Perturbation signature per guide = mean(log-normalized expression) minus
   own-dataset NT-control mean (log fold-change vector).
4. Reference signatures aggregated to `target_gene` (NT = average of its 10
   control guides, ~0 by construction).
5. Score(gene, query guide) = Pearson correlation of the two signature vectors.
6. One-to-one assignment (reference has exactly one guide per gene, so a
   bijective gene-to-guide matching is the natural constraint) by Hungarian
   maximization of total score.
7. `score` = correlation of the assigned pair; `runner_up_score` = best
   correlation of the same gene with any *other* query guide; `confidence` =
   softmax probability of the assigned guide over all 42 query-guide scores for
   that gene (uniform chance level = 1/42 ~ 0.024).

## Headline result
| target_gene | query_guide_id | score | runner_up_score | confidence |
|---|---|---|---|---|
| **PABPC1** | **guide18** | 0.532 | 0.339 | 0.035 |
| **NUDT21** | **guide13** | 0.636 | 0.513 | 0.038 |
| **LEO1** | **guide14** | 0.314 | 0.223 | 0.029 |

- **PABPC1 -> guide18**: strongest of the three; guide18's own best reference hit
  is also PABPC1 (mutually consistent), margin +0.194,
  stable under 50% subsampling.
- **NUDT21 -> guide13**: mutually consistent, margin +0.123,
  stable under subsampling. Main confounder is CPSF6 (same 3'-end processing
  pathway; r = 0.513 with the query signature).
- **LEO1 -> guide14**: weakest of the three. Margin +0.091,
  but from the query side guide14 scores CTR9 slightly higher (0.370 vs 0.314) —
  CTR9 and LEO1 are both PAF1-complex subunits with near-degenerate knockdown
  effects — and the LEO1 row flips under 50% subsampling (below). Treat this
  assignment as tentative: the data support "a PAF1-complex guide" more strongly
  than they distinguish LEO1 from CTR9.

## Ambiguity quantification
- **Margins**: 23/34 rows have negative margin (some *other* query guide
  correlates better with the reference signature than the globally assigned one).
  This is inherent to forcing a one-to-one assignment when several same-pathway
  genes (CPSF1/2/3/4/6, CSTF1/3, SCAF1, SYMPK, CTR9/LEO1/PAF1/CDC73) have strongly
  correlated knockdown signatures. Median score = 0.321, median
  margin = -0.035.
- **Subsample stability**: recomputing every signature on a random 50% of each
  group's cells (seed 0) and re-running the assignment leaves
  76.5% of rows unchanged. Unstable rows:
  CDC73, CPSF2, CPSF3, CTR9, LEO1, PAF1, PHF3, SF3A1 — concentrated in the low-margin
  PAF1-complex block (CDC73, CTR9, LEO1, PAF1) and CPSF block (CPSF2, CPSF3),
  plus SF3A1 and PHF3. PABPC1 and NUDT21 rows are stable.
- **Negative-control floor**: the 10 query NT guides show only noise-level
  correlations with all perturbation signatures (max |r| = 0.141),
  and the NT row is still assigned an NT guide — the scores above reflect
  perturbation biology rather than technical artifacts.
- **Unresolvable block**: PHF3 (best r = 0.066), SCAF8 (best r = 0.054) and query
  guide32 (best r = 0.052) all sit at noise level. The query carries 32 perturbed
  guides versus 33 reference genes, so one gene's guide is absent from the query;
  the data cannot decide whether PHF3 or SCAF8 is the missing gene, nor whether
  guide32 is a failed perturbation of the other. Both genes were assigned residual
  NT guides by the global optimizer; those two rows are effectively unmapped and
  are flagged as such.

## Validation performed (lightweight, independent)
1. 50%-subsample re-derivation of every signature and of the full assignment —
   an independent stochastic rerun of the pipeline (see stability above).
2. Negative-control sanity: NT row assigned an NT guide, and all NT query guides
   at noise-level correlation (above).

## Full gene -> query-guide mapping (output/guide_mapping.csv, 34 rows)
| target_gene | query_guide_id | score | runner_up_score | confidence |
|---|---|---|---|---|
| CDC73 | guide10 | 0.114 | 0.344 | 0.023 |
| CPSF1 | guide4 | 0.671 | 0.546 | 0.036 |
| CPSF2 | guide21 | 0.473 | 0.613 | 0.030 |
| CPSF3 | guide3 | 0.523 | 0.616 | 0.031 |
| CPSF3L | guide1 | 0.122 | 0.273 | 0.025 |
| CPSF4 | guide29 | 0.551 | 0.643 | 0.032 |
| CPSF6 | guide19 | 0.519 | 0.526 | 0.033 |
| CPSF7 | guide2 | 0.203 | 0.082 | 0.028 |
| CSTF1 | guide16 | 0.468 | 0.536 | 0.030 |
| CSTF2 | guide31 | 0.100 | 0.239 | 0.023 |
| CSTF3 | guide22 | 0.507 | 0.572 | 0.032 |
| CTR9 | guide28 | 0.310 | 0.370 | 0.028 |
| FIP1L1 | guide17 | 0.463 | 0.464 | 0.031 |
| LEO1 | guide14 | 0.314 | 0.223 | 0.029 |
| NT | NT-I8 | 0.018 | 0.137 | 0.023 |
| NUDT21 | guide13 | 0.636 | 0.513 | 0.038 |
| PABPC1 | guide18 | 0.532 | 0.339 | 0.035 |
| PABPN1 | guide24 | 0.399 | 0.276 | 0.031 |
| PAF1 | guide12 | 0.260 | 0.361 | 0.026 |
| PAPOLA | guide7 | 0.333 | 0.272 | 0.030 |
| PCF11 | guide20 | 0.559 | 0.501 | 0.033 |
| PHF3 | NT-I7 | 0.042 | 0.066 | 0.026 |
| RBBP6 | guide9 | 0.453 | 0.447 | 0.031 |
| RPAP2 | guide15 | 0.198 | 0.304 | 0.025 |
| RPRD1A | guide8 | 0.116 | 0.152 | 0.025 |
| RPRD1B | guide5 | 0.328 | 0.366 | 0.028 |
| SCAF1 | guide26 | 0.432 | 0.519 | 0.029 |
| SCAF11 | guide11 | 0.135 | 0.111 | 0.026 |
| SCAF4 | guide6 | 0.258 | 0.292 | 0.027 |
| SCAF8 | NT-I3 | 0.040 | 0.054 | 0.026 |
| SF3A1 | guide23 | 0.224 | 0.305 | 0.026 |
| SRSF3 | guide30 | 0.289 | 0.232 | 0.029 |
| SYMPK | guide25 | 0.435 | 0.511 | 0.029 |
| THOC5 | guide27 | 0.286 | 0.308 | 0.027 |

Margins (score - runner_up_score):
| target_gene | query_guide_id | score | runner_up_score | margin |
|---|---|---|---|---|
| CDC73 | guide10 | 0.114 | 0.344 | -0.229 |
| CPSF1 | guide4 | 0.671 | 0.546 | 0.125 |
| CPSF2 | guide21 | 0.473 | 0.613 | -0.140 |
| CPSF3 | guide3 | 0.523 | 0.616 | -0.093 |
| CPSF3L | guide1 | 0.122 | 0.273 | -0.151 |
| CPSF4 | guide29 | 0.551 | 0.643 | -0.092 |
| CPSF6 | guide19 | 0.519 | 0.526 | -0.008 |
| CPSF7 | guide2 | 0.203 | 0.082 | 0.120 |
| CSTF1 | guide16 | 0.468 | 0.536 | -0.068 |
| CSTF2 | guide31 | 0.100 | 0.239 | -0.139 |
| CSTF3 | guide22 | 0.507 | 0.572 | -0.065 |
| CTR9 | guide28 | 0.310 | 0.370 | -0.060 |
| FIP1L1 | guide17 | 0.463 | 0.464 | -0.001 |
| LEO1 | guide14 | 0.314 | 0.223 | 0.091 |
| NT | NT-I8 | 0.018 | 0.137 | -0.119 |
| NUDT21 | guide13 | 0.636 | 0.513 | 0.123 |
| PABPC1 | guide18 | 0.532 | 0.339 | 0.194 |
| PABPN1 | guide24 | 0.399 | 0.276 | 0.123 |
| PAF1 | guide12 | 0.260 | 0.361 | -0.101 |
| PAPOLA | guide7 | 0.333 | 0.272 | 0.061 |
| PCF11 | guide20 | 0.559 | 0.501 | 0.058 |
| PHF3 | NT-I7 | 0.042 | 0.066 | -0.024 |
| RBBP6 | guide9 | 0.453 | 0.447 | 0.006 |
| RPAP2 | guide15 | 0.198 | 0.304 | -0.106 |
| RPRD1A | guide8 | 0.116 | 0.152 | -0.036 |
| RPRD1B | guide5 | 0.328 | 0.366 | -0.038 |
| SCAF1 | guide26 | 0.432 | 0.519 | -0.088 |
| SCAF11 | guide11 | 0.135 | 0.111 | 0.023 |
| SCAF4 | guide6 | 0.258 | 0.292 | -0.034 |
| SCAF8 | NT-I3 | 0.040 | 0.054 | -0.014 |
| SF3A1 | guide23 | 0.224 | 0.305 | -0.081 |
| SRSF3 | guide30 | 0.289 | 0.232 | 0.058 |
| SYMPK | guide25 | 0.435 | 0.511 | -0.076 |
| THOC5 | guide27 | 0.286 | 0.308 | -0.022 |

## Supplementary: query-guide -> best reference gene view
Explanatory only (output/supplementary_query_guide_mapping.csv); the required
machine-readable artifact is guide_mapping.csv with one row per reference gene.

| query_guide_id | best_target_gene | score | runner_up_target_gene | runner_up_score |
|---|---|---|---|---|
| guide1 | SCAF1 | 0.174 | THOC5 | 0.160 |
| guide2 | CPSF7 | 0.203 | NUDT21 | 0.144 |
| guide3 | CPSF2 | 0.537 | CPSF1 | 0.533 |
| guide4 | CPSF1 | 0.671 | CPSF4 | 0.643 |
| guide5 | PCF11 | 0.455 | CPSF3 | 0.425 |
| guide6 | SCAF4 | 0.258 | NUDT21 | 0.192 |
| guide7 | PAPOLA | 0.333 | CPSF3 | 0.247 |
| guide8 | CPSF2 | 0.287 | CPSF3 | 0.257 |
| guide9 | RBBP6 | 0.453 | PCF11 | 0.429 |
| guide10 | CPSF1 | 0.162 | CPSF4 | 0.160 |
| guide11 | SCAF1 | 0.233 | SYMPK | 0.194 |
| guide12 | CPSF4 | 0.385 | CPSF1 | 0.385 |
| guide13 | NUDT21 | 0.636 | CPSF6 | 0.526 |
| guide14 | CTR9 | 0.370 | SCAF1 | 0.364 |
| guide15 | SCAF1 | 0.381 | SYMPK | 0.314 |
| guide16 | SCAF1 | 0.519 | SYMPK | 0.475 |
| guide17 | CPSF4 | 0.517 | CPSF1 | 0.486 |
| guide18 | PABPC1 | 0.532 | SCAF1 | 0.376 |
| guide19 | CPSF6 | 0.519 | NUDT21 | 0.513 |
| guide20 | PCF11 | 0.559 | CPSF4 | 0.503 |
| guide21 | CPSF2 | 0.473 | CPSF3 | 0.454 |
| guide22 | CPSF4 | 0.538 | CPSF1 | 0.534 |
| guide23 | SCAF1 | 0.379 | CSTF1 | 0.318 |
| guide24 | PABPN1 | 0.399 | CSTF1 | 0.230 |
| guide25 | CPSF3 | 0.456 | CPSF1 | 0.454 |
| guide26 | SCAF1 | 0.432 | SYMPK | 0.374 |
| guide27 | SCAF1 | 0.418 | SYMPK | 0.367 |
| guide28 | CSTF1 | 0.371 | SCAF1 | 0.348 |
| guide29 | CPSF4 | 0.551 | CPSF1 | 0.546 |
| guide30 | CSTF1 | 0.384 | RPRD1B | 0.366 |
| guide31 | CSTF1 | 0.222 | CPSF4 | 0.217 |
| guide32 | RPRD1B | 0.052 | CSTF1 | 0.044 |
| NT-I1 | CPSF3L | 0.052 | SF3A1 | 0.030 |
| NT-I2 | SCAF4 | 0.031 | PABPC1 | 0.029 |
| NT-I3 | SCAF8 | 0.040 | CPSF7 | 0.017 |
| NT-I4 | CPSF6 | 0.092 | NUDT21 | 0.087 |
| NT-I5 | CPSF2 | 0.103 | SYMPK | 0.099 |
| NT-I6 | PABPC1 | 0.049 | SYMPK | 0.043 |
| NT-I7 | PHF3 | 0.042 | SCAF8 | 0.020 |
| NT-I8 | CPSF2 | 0.065 | CPSF1 | 0.065 |
| NT-I9 | PHF3 | 0.015 | RPAP2 | 0.013 |
| NT-I10 | PABPC1 | 0.063 | PCF11 | 0.060 |

## Limitations
- One guide per perturbed gene in the reference: no replicate guides exist to
  estimate guide-level noise within the reference; ambiguity was therefore
  quantified via margins, subsample stability and control floors.
- Pearson correlation of mean log-fold-change signatures is a simple summary;
  pathway-level signature similarity (CPSF complex, CSTF complex, PAF1 complex)
  sets a hard ceiling on gene-level identifiability across the cell-type shift.
- The LEO1 row (and the PAF1 block generally) is best interpreted as a
  pathway-level assignment; rows with large negative margins (CDC73, CPSF2,
  CPSF3L, CSTF2, NT, RPAP2) are low-confidence and are reported as such in the
  CSV via runner_up_score/confidence.
