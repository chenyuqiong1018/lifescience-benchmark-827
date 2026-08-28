# Structural Model Ranking Report

## Task and inputs

Rank all supplied structural models exactly per `inputs/SCORING_RULE.md`, using only the
files under `inputs/`:

- `inputs/model_metrics.csv` - per-model comparison metrics against one frozen reference:
  `tm_score`, `lddt`, `rmsd_a`, `aligned_residues`, `chain_mapping_complete`.
- `inputs/residue_errors.csv` - per-model per-residue-interval `mean_error_a` values.
- `inputs/SCORING_RULE.md` - the ranking rule.

Per the scoring rule, this fixture is a benchmark-informed local extension: the supplied
values are already-computed comparison metrics against one frozen reference. Accordingly,
this report makes **no coordinate-level, interface-geometry, or experimental claims**; it
only orders models using the provided numbers.

## Ranking rule (implemented exactly)

Models are ordered by the following stable tuple, evaluated in order:

1. complete chain mapping before incomplete mapping;
2. higher `tm_score`;
3. higher `lddt`;
4. lower `rmsd_a`;
5. more `aligned_residues`;
6. lower critical-region mean error;
7. lexical `model_id` as final tie break.

The critical region is residues 181-240. `critical_residue_risk` is the length-weighted
mean of `mean_error_a` over rows in `residue_errors.csv` overlapping that interval. In
this fixture each critical-region row covers the whole interval, so the reported risk
equals that row's value.

Output column semantics (as mandated):

- `global_score = tm_score` (no invented weighted composite);
- `interface_score = 1` for complete chain mapping, `0` otherwise (this is a
  mapping-completeness flag defined by the rule, not a measured interface quality);
- `critical_residue_risk` in angstrom;
- `decision`: `preferred` for rank 1; `alternate` for other mapping-complete models;
  `reject_incomplete_mapping` for incomplete mapping.

Implementation: `output/analysis.py` (reads only `inputs/`, writes
`output/model_ranking.csv`).

## Result

| rank | model_id | global_score (tm_score) | interface_score (mapping complete) | critical_residue_risk (A) | decision |
|------|----------|-------------------------|------------------------------------|---------------------------|----------|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

Supporting input metrics (from `model_metrics.csv`): model_A lddt 0.76, rmsd_a 2.4,
aligned_residues 238; model_B lddt 0.79, rmsd_a 2.0, aligned_residues 211; model_C
lddt 0.65, rmsd_a 3.8, aligned_residues 240.

Ranking derivation:

- Criterion 1 partitions the field: model_A and model_C have complete chain mapping,
  model_B does not, so A and C must both rank above B regardless of any other metric.
- Criterion 2 then separates the two complete models: model_A (tm_score 0.81) ranks
  above model_C (tm_score 0.62). Criteria 3-7 were not needed between them.
- model_B is last and receives `reject_incomplete_mapping`.
- Criteria 6 (critical-region error) and 7 (lexical id) were never reached as tie
  breakers in this dataset, but are implemented and would apply to exact ties.

## How chain-mapping completeness affects the ranking

Chain-mapping completeness is the first, dominant criterion of the rule, so it outranks
every quality metric. Its effect here is decisive:

- model_B has the best secondary quality numbers of the three on two diagnostics
  (lddt 0.79, rmsd_a 2.0) - better than both complete models on those fields - yet it
  ranks last and is rejected because its chain mapping is incomplete. An incompletely
  mapped model cannot be compared position-by-position with the reference across all
  chains, so the rule treats its other numbers as insufficient for selection.
- Completeness is exposed in the required schema as `interface_score` (1/0). model_A and
  model_C get 1; model_B gets 0.
- Among mapping-complete models only, the quality metrics decide: model_A becomes
  `preferred` and model_C becomes `alternate` purely on their relative metrics (here,
  tm_score at criterion 2).

## How critical-region uncertainty affects the ranking

The critical region (residues 181-240) is summarized per model as the length-weighted
mean of `mean_error_a` over overlapping rows, yielding a single uncertainty value in
angstrom:

- model_A: 4.8 A; model_C: 4.1 A; model_B: 8.5 A.

Its role in the ranking:

- It is criterion 6, applied only after mapping completeness, tm_score, lddt, rmsd_a and
  aligned_residues. It therefore acts as a tie-breaker that penalizes models whose
  functionally important region (181-240) is less certain, without overriding global
  similarity or coverage. In this dataset no earlier criteria tied, so critical-region
  risk did not change the order - but it is implemented and would reorder exact ties in
  favor of lower critical-region error.
- Independently of rank mechanics, the values corroborate the decisions: model_B carries
  the largest critical-region error (8.5 A), consistent with its rejection, although the
  binding reason for rejecting model_B is incomplete chain mapping, not its error value.
  model_C has the lowest critical-region risk (4.1 A) but remains `alternate` because its
  global similarity (tm_score 0.62) is lower than model_A's; the rule intentionally does
  not let local-region certainty override global similarity.

## Caveats

- All statements are restricted to the provided pre-computed comparison metrics against
  one frozen reference. No coordinate-level, interface-structure, or experimental
  properties are asserted, since none are present in the inputs.
- `interface_score` is the rule-defined 1/0 mapping-completeness flag, not an interface
  quality measurement.
- No composite score was invented: `global_score` is exactly `tm_score`.
