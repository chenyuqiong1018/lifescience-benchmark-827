Done. All deliverables are in `output/` and validated.

## Headline answer
| Target | Query guide | Evidence |
|---|---|---|
| **NUDT21** | **guide13** | score 0.484 (highest of all), runner-up CPSF6 0.377, bootstrap conf 1.00 |
| **LEO1** | **guide14** | score 0.249, runner-up CTR9 0.198 (same PAF1-complex — expected confuser), conf 0.98 |
| **PABPC1** | **guide18** (primary), guide15 (secondary) | guide18: score 0.326, gap 0.190, conf 1.00; guide15 also ranks PABPC1 first but weakly (score 0.166, gap 0.042, conf 0.84) — flagged as ambiguous |

## Method (single primary method)
- Aligned both datasets to 11,858 shared genes; per-dataset CP10k+log1p; 2,000 HVGs; PCA + KMeans (k=10) cell-state clusters fit **within each dataset**.
- Cell-type-adjusted effect per guide: per-gene OLS `expr ~ cluster dummies + perturbation indicator` against same-dataset NT cells — necessary because the joint-cluster table shows a severe cell-type shift (9/10 clusters are batch-exclusive).
- Mapping = Pearson correlation of effect vectors with the **candidate target gene excluded**; confidence = bootstrap argmax agreement (100 resamples of guide + NT cells).

## Leakage prevention
- Query guide IDs treated as opaque (never parsed); NT-name overlap used only to define controls.
- 32/33 reference targets (all but CDC73) are structurally absent from the query matrix, so matching uses downstream effects only; additionally every candidate gene is excluded from its own similarity computation.

## Validation & ambiguity
- NT pseudo-perturbation negative controls top out at sim ≤ 0.17 vs. real scores up to 0.48.
- Reference split-half self-mapping ceiling: 25/33 correct (0.76); all errors are tiny ref guides or same-complex subunits (PAF1↔CTR9, CPSF2↔CPSF3).
- Independent Spearman re-check: 23/32 argmax agreement; all 9 flips are low-confidence small guides; **all three headline targets unchanged**.
- Ambiguity quantified per guide via `runner_up_score`, `gap`, bootstrap `confidence`: 15/32 guides ≥0.90 conf, 6 guides ≤0.50 conf (all small groups), 13/32 near-ties (gap <0.03).

**Files:** `output/guide_mapping.csv` (32 rows, required columns first), `output/analysis.py`, `output/report.md`, plus `output/results.json` and `output/validate_spearman.py`.