# Model Ranking Report

## Scope and inputs

Per `inputs/SCORING_RULE.md`, this is a benchmark-informed local extension: the supplied
values are already-computed comparison metrics against one frozen reference. Only the two
input files were used:

- `inputs/model_metrics.csv` — `tm_score`, `lddt`, `rmsd_a`, `aligned_residues`,
  `chain_mapping_complete` for models A, B, C.
- `inputs/residue_errors.csv` — per-interval `mean_error_a` rows (Å) for residues 1–180
  and 181–240 per model.

No coordinate files are present, so **no coordinate-level, interface, or experimental
properties are claimed**. `global_score` is `tm_score` directly (no invented composite),
and `interface_score` is defined by the rule purely as the chain-mapping completeness
indicator (1/0), not as a measured interface quality.

## Ranking rule applied

Stable tuple, in priority order: (1) complete chain mapping before incomplete;
(2) higher `tm_score`; (3) higher `lddt`; (4) lower `rmsd_a`; (5) more
`aligned_residues`; (6) lower critical-region mean error; (7) lexical `model_id`.

`critical_residue_risk` = length-weighted mean of `mean_error_a` over rows overlapping
residues 181–240. In this fixture each model has exactly one row covering the whole
181–240 interval, so each reported risk equals that row's value.

## Result

| rank | model_id | global_score (TM-score) | interface_score | critical_residue_risk (Å) | decision |
|---|---|---|---|---|---|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

Supporting metrics (secondary diagnostics, not mixed into the score):

| model | lddt | rmsd_a (Å) | aligned_residues | chain_mapping_complete |
|---|---|---|---|---|
| model_A | 0.76 | 2.4 | 238 | true |
| model_B | 0.79 | 2.0 | 211 | false |
| model_C | 0.65 | 3.8 | 240 | true |

## How chain-mapping completeness affects the ranking

Completeness is the **first, dominant key** of the tuple, so it overrides every
metric-level comparison. This is decisive here:

- `model_B` has the best `rmsd_a` (2.0 Å), the highest `lddt` (0.79), and a higher
  `tm_score` (0.74) than model_C (0.62). On metrics alone it would rank second, but its
  `chain_mapping_complete = false` demotes it below both complete models and assigns
  `decision = reject_incomplete_mapping`. An incomplete mapping means the comparison
  metrics were computed over a chain assignment that does not cover all chains, so the
  metrics cannot be trusted as a full-model comparison regardless of their values.
- Among the two complete models the tie-break chain never leaves key 2: `model_A`
  (TM-score 0.81) beats `model_C` (0.62) outright, so keys 3–7 are never consulted for
  that pair.

## How critical-region uncertainty affects the ranking

The critical region is residues 181–240. Its risk enters as key 6 — **after** mapping
completeness, TM-score, lDDT, RMSD, and aligned coverage — so it only separates models
that are otherwise tied. In this dataset no pair reaches key 6, so the critical-region
risk does not change the ordering; it is reported but is not the deciding factor.
Notably, `model_C` has the *lowest* critical-region risk (4.1 Å vs 4.8 Å for model_A)
yet still ranks below model_A, because global TM-score is compared first.
`model_B`'s high risk (8.5 Å) compounds its rejection but is not the stated reason for
it — the rejection is solely due to incomplete mapping.

Regarding *uncertainty* in the risk estimate itself: the rule defines risk as a
length-weighted mean over rows overlapping 181–240. Here each model contributes exactly
one row spanning the entire interval (181–240, length 60), so the weighting is trivial
and the reported risk equals that single row's value with no ambiguity. If the fixture
instead supplied partially overlapping rows (e.g., 181–200 and 201–240 with different
errors), the length-weighting would determine the aggregate, and any row not fully
covering the critical interval would introduce uncertainty about which residues' errors
are represented; the rule fixes that by weighting strictly by overlap length. No such
uncertainty exists in the present inputs.

## Reproducibility

`output/analysis.py` re-reads the two input CSVs, applies the tuple exactly as above,
and rewrites `output/model_ranking.csv`. Run: `python output/analysis.py`.

Scientific basis per the rule: TM-score is a length-normalized global structure
comparison measure (Zhang & Skolnick, 2004, DOI 10.1002/prot.20264); lDDT is a
superposition-free local distance-difference measure (Mariani et al., 2013,
DOI 10.1093/bioinformatics/btt473). RMSD and aligned coverage are retained only as
secondary diagnostics.
