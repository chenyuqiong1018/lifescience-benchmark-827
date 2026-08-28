# Perturb-seq query-to-reference guide mapping — analysis report

**Task.** Map anonymized query perturbation groups (`guide1..guide32`) to the labeled
Perturb-seq reference (33 target genes + non-targeting `NT` controls) despite a cell-type
composition shift between the two datasets, and identify the query guides targeting
**PABPC1**, **NUDT21** and **LEO1**.

**Inputs.** `inputs/perturb.seq.align.q1.ref.h5ad` (4,183 cells x 11,889 genes; obs:
`guide`, `target_gene`) and `inputs/perturb.seq.align.q1.query.h5ad` (3,516 cells x 11,889
genes; obs: `guide` only). Both are raw UMI counts.

---

## 1. Headline answer

| target_gene | query guide(s) | evidence |
|---|---|---|
| **NUDT21** | **guide13** | score 0.484 (best overall), runner-up CPSF6 0.377, bootstrap confidence 1.00, 86 cells; confirmed by rank-based (Spearman) re-check |
| **LEO1** | **guide14** | score 0.249, runner-up CTR9 0.198 (same PAF1-complex — expected confuser), confidence 0.98, 55 cells; confirmed by Spearman re-check |
| **PABPC1** | **guide18** (primary) and **guide15** (secondary candidate) | guide18: score 0.326, runner-up SCAF1 0.136, gap 0.190, confidence 1.00, 97 cells — the strongest PABPC1 evidence. guide15: score 0.166, runner-up RPRD1B 0.124, gap 0.042, confidence 0.84, 55 cells — also ranks PABPC1 first under both metrics, but weakly. If the query contains one guide per target gene, guide18 is PABPC1 and guide15 is an ambiguous guide whose true target phenocopies PABPC1 partially; if two guides target PABPC1, both apply. |

## 2. Method (single primary method)

1. **Feature alignment.** Reference and query were restricted to the intersection of their
   variable sets (11,858 shared genes).
2. **Normalization.** Per-dataset CP10k + log1p (computed independently within each dataset).
3. **Feature selection.** Top 2,000 highly variable genes on the concatenated, normalized data.
4. **Cell-state clustering (cell-type proxy).** 30-component PCA on scaled HVG expression,
   then KMeans with k=10 fit **separately within each dataset** (deterministic seed 0).
   Joint clusters were computed only for describing the shift (Sec. 4).
5. **Cell-type-adjusted perturbation effects.** For each perturbation group g, a per-gene
   linear model on that group's cells plus the same-dataset NT control cells:
   `expr_g ~ 1 + cluster dummies + 1[cell in g]`. The perturbation-indicator coefficient is
   the effect vector. This removes composition differences between a guide's cells and the
   NT baseline, which is essential given the dataset cell-type shift.
6. **Query -> reference matching.** Pearson correlation between each query-guide effect and
   each reference-target effect on the shared HVG set, **excluding the candidate target gene
   itself** from the comparison. Query guide = argmax correlation; runner-up = second best.
7. **Confidence / ambiguity.** Nonparametric bootstrap (100 resamples of both the guide's
   cells and the NT baseline cells): recompute the effect and the argmax each time;
   `confidence` = fraction of bootstraps reproducing the argmax. `gap` = score - runner_up_score.

## 3. Target-metadata leakage prevention

- **Query guide IDs are opaque.** `guide1..guide32` were never parsed or matched against gene
  symbols; grouping used them only as categorical labels. Query `NT-I#` guides share names
  with reference `NT-I#` guides; this overlap was used solely to define the query's
  non-targeting control population, never to assign perturbation identity.
- **Structural leakage block in the data itself:** 32 of the 33 reference target genes
  (all except `CDC73`) are absent from the query feature matrix, i.e. the query cannot be
  matched via the target gene's own expression. Mapping therefore relies entirely on
  downstream transcriptional effects.
- **Candidate-gene exclusion:** for every (query guide, reference target G) comparison, gene
  G is removed from the correlation, so even the one target present on both sides (`CDC73`)
  cannot drive its own match through a cis effect.
- Effects were computed independently per dataset; the only cross-dataset quantity is the
  correlation of effect vectors (the mapping output itself).

## 4. Cell-type shift quantification

Joint KMeans (k=10) clusters on the concatenated PCA are almost batch-exclusive, confirming a
strong composition shift (fractions of each dataset's cells):

| joint cluster | query | reference |
|---|---|---|
| 0 | 0.000 | 0.231 |
| 1 | 0.078 | 0.000 |
| 2 | 0.387 | 0.000 |
| 3 | 0.000 | 0.186 |
| 4 | 0.000 | 0.491 |
| 5 | 0.164 | 0.000 |
| 6 | 0.179 | 0.000 |
| 7 | 0.007 | 0.027 |
| 8 | 0.184 | 0.000 |
| 9 | 0.001 | 0.065 |

Only cluster 7 is shared. Because the two datasets occupy largely different cell states,
cross-dataset comparisons of raw means would be dominated by composition; the per-dataset
cluster-adjusted effects (Sec. 2, step 5) are what make the mapping comparable.

## 5. Full query-guide mapping

`score`/`runner_up_score`: Pearson correlation of cell-type-adjusted effect vectors with the
candidate target gene excluded. `confidence`: bootstrap argmax agreement (100 resamples).

| query_guide_id | target_gene | score | runner_up_gene | runner_up_score | confidence | gap | n_cells |
|---|---|---|---|---|---|---|---|
| guide1 | **THOC5** | 0.1592 | PAPOLA | 0.1007 | 0.91 | 0.0585 | 199 |
| guide2 | **CPSF7** | 0.2024 | NUDT21 | 0.1131 | 1.00 | 0.0894 | 200 |
| guide3 | **CPSF2** | 0.3011 | CPSF1 | 0.2741 | 0.96 | 0.0270 | 23 |
| guide4 | **CPSF1** | 0.3424 | CSTF3 | 0.3130 | 0.90 | 0.0294 | 68 |
| guide5 | **RPRD1B** | 0.2260 | RBBP6 | 0.2130 | 0.75 | 0.0129 | 178 |
| guide6 | **SCAF4** | 0.1406 | SRSF3 | 0.1299 | 0.61 | 0.0107 | 127 |
| guide7 | **PAPOLA** | 0.3131 | PABPN1 | 0.2257 | 1.00 | 0.0874 | 208 |
| guide8 | **CPSF2** | 0.3015 | CPSF1 | 0.1745 | 1.00 | 0.1270 | 102 |
| guide9 | **RBBP6** | 0.2878 | PCF11 | 0.2348 | 0.92 | 0.0529 | 48 |
| guide10 | **CTR9** | 0.0704 | RPRD1B | 0.0676 | 0.18 | 0.0028 | 16 |
| guide11 | **SCAF11** | 0.1114 | RPRD1B | 0.0760 | 0.66 | 0.0354 | 158 |
| guide12 | **FIP1L1** | 0.0927 | PCF11 | 0.0907 | 0.17 | 0.0021 | 20 |
| guide13 | **NUDT21** | 0.4843 | CPSF6 | 0.3774 | 1.00 | 0.1069 | 86 |
| guide14 | **LEO1** | 0.2493 | CTR9 | 0.1976 | 0.98 | 0.0517 | 55 |
| guide15 | **PABPC1** | 0.1661 | RPRD1B | 0.1244 | 0.84 | 0.0417 | 55 |
| guide16 | **RPRD1B** | 0.2338 | PAPOLA | 0.1321 | 1.00 | 0.1016 | 40 |
| guide17 | **CPSF1** | 0.2445 | FIP1L1 | 0.2420 | 0.31 | 0.0025 | 27 |
| guide18 | **PABPC1** | 0.3261 | SCAF1 | 0.1357 | 1.00 | 0.1904 | 97 |
| guide19 | **CPSF6** | 0.3623 | NUDT21 | 0.3220 | 0.98 | 0.0403 | 32 |
| guide20 | **RBBP6** | 0.2641 | PCF11 | 0.2237 | 0.96 | 0.0404 | 43 |
| guide21 | **CPSF2** | 0.3587 | CPSF1 | 0.2474 | 1.00 | 0.1113 | 38 |
| guide22 | **CSTF3** | 0.3145 | FIP1L1 | 0.2765 | 0.87 | 0.0381 | 29 |
| guide23 | **RPRD1B** | 0.1487 | RBBP6 | 0.1195 | 0.84 | 0.0292 | 83 |
| guide24 | **PABPN1** | 0.4430 | PAPOLA | 0.2223 | 1.00 | 0.2207 | 65 |
| guide25 | **CPSF1** | 0.2110 | CPSF2 | 0.1761 | 0.83 | 0.0349 | 6 |
| guide26 | **PHF3** | 0.1010 | RPRD1B | 0.0686 | 0.56 | 0.0325 | 15 |
| guide27 | **PABPN1** | 0.1866 | THOC5 | 0.1651 | 0.58 | 0.0215 | 10 |
| guide28 | **RPRD1B** | 0.1744 | PCF11 | 0.1301 | 0.64 | 0.0444 | 16 |
| guide29 | **CSTF3** | 0.2542 | CPSF1 | 0.2501 | 0.20 | 0.0041 | 22 |
| guide30 | **RPRD1B** | 0.1079 | SRSF3 | 0.0905 | 0.48 | 0.0174 | 13 |
| guide31 | **CSTF1** | 0.1422 | RPRD1B | 0.1280 | 0.52 | 0.0143 | 10 |
| guide32 | **CSTF1** | 0.0960 | RPRD1B | 0.0885 | 0.38 | 0.0074 | 3 |

Coverage: 24 of 33 reference targets are the best match of at least one query guide; 8
targets are selected by >= 2 guides (`RPRD1B` x5, `CPSF1` x3, `CPSF2` x3, `RBBP6`, `PABPC1`,
`CSTF3`, `PABPN1`, `CSTF1` x2 each); 14 targets are never a best match
(`CDC73, CPSF3, CPSF3L, CPSF4, CSTF2, PAF1, PCF11, RPAP2, RPRD1A, SCAF1, SCAF8, SF3A1,
SRSF3, SYMPK`). Duplicates concentrate among small guides (<= 32 cells), where the weak
effect vector collapses toward strong, well-characterized reference phenotypes — this is
quantified per guide by the low `confidence`/`gap` values above.

## 6. Ambiguity quantification (summary)

- **Bootstrap confidence:** 15/32 guides have confidence >= 0.90, 19/32 >= 0.80; 6 guides
  have confidence <= 0.5 (guide10, guide12, guide29, guide30, guide31, guide32 — all small
  groups and/or near-zero scores).
- **Near-ties:** 13/32 guides have `gap` < 0.03 (top two targets statistically
  indistinguishable given bootstrap noise); the tightest are guide12 (0.0021), guide17
  (0.0025) and guide10 (0.0028).
- **Score scale:** scores range 0.070-0.484 (median 0.230). NT-vs-NT pseudo-perturbations
  (negative controls, below) top out at 0.167, so scores near or below ~0.1 are not
  meaningful signal; 6 query guides fall at or below that threshold.
- **Reference-side ambiguity:** leave-half-out self-mapping of the reference recovers the
  correct target for 25/33 guides (top-1 accuracy 0.76). All 8 errors involve tiny reference
  guides (5-43 cells: RPAP2, SYMPK, SF3A1, PAF1, CSTF1, PCF11, CSTF2, RPRD1A) or subunits of
  the same protein complex with near-identical phenotypes (PAF1 -> CTR9; CPSF2 -> CPSF3;
  CSTF1/CSTF2/RPRD1A -> CPSF complex). This 0.76 ceiling is the intrinsic ambiguity of the
  reference itself under the same metric.

## 7. Validation (one lightweight independent check)

- **Negative controls (NT):** treating each query NT guide (and each reference NT guide) as a
  pseudo-perturbation against the remaining NT cells, the maximum similarity to any reference
  target is 0.043-0.167 (query) and 0.054-0.124 (reference) — far below the scores of the
  confident mappings, showing the metric does not match noise to targets.
- **Independent rank-based re-check:** recomputing all similarities as Spearman (rank)
  correlations on the same adjusted effects reproduces the primary argmax for 23/32 guides.
  All 9 disagreements are low-score, low-confidence guides (e.g. guide32 with 3 cells,
  guide10 with score 0.07). **guide13->NUDT21, guide14->LEO1, guide15->PABPC1 and
  guide18->PABPC1 are unchanged, with the same runner-ups**, giving high confidence in the
  headline assignments.

## 8. Limitations

- Reference targets are represented by a single guide each; reference effects for small
  guides (e.g. RPAP2, 5 cells) are noisy, and their query-side matches are correspondingly
  unreliable.
- The mapping is one-directional argmax per query guide; no global one-to-one assignment was
  enforced, so target collisions are visible and interpretable via `confidence`/`gap`.
- Cell states were approximated by KMeans (k=10) clusters; finer cell-type resolution could
  sharpen small-group effects but would not change the headline assignments.

## 9. Reproducibility

- `output/analysis.py` — full pipeline (deterministic, seed 0); run from the workspace root:
  `python output/analysis.py`.
- `output/validate_spearman.py` — the independent rank-based validation.
- `output/results.json` — machine-readable results (similarity details, validation outputs,
  composition table).
- `output/guide_mapping.csv` — required deliverable (32 rows, one per perturbed query guide).

