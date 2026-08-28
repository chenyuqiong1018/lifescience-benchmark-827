# Structural model ranking report

Generated strictly from `inputs/SCORING_RULE.md`, `inputs/model_metrics.csv`, and
`inputs/residue_errors.csv` by `output/analysis.py`. Per the rule, the supplied values
are pre-computed comparison metrics against one frozen reference; no coordinate-level,
interface, or experimental properties are claimed anywhere in this report.

## Final ranking (`output/model_ranking.csv`)

| rank | model_id | global_score (= tm_score) | interface_score | critical_residue_risk (Å) | decision |
|-----:|----------|--------------------------:|----------------:|--------------------------:|----------|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

Column definitions per the rule: `global_score` is `tm_score` itself (no weighted
composite was invented); `interface_score` is 1 for complete chain mapping and 0
otherwise; `critical_residue_risk` is the length-weighted mean of `mean_error_a` over
rows overlapping residues 181–240 (in this fixture each model has exactly one row
covering the whole interval, so the risk equals that row's value).

## How chain-mapping completeness affects the ranking

Completeness is the first key of the stable sorting tuple, so it dominates every
structural-quality metric. `model_B` has the best RMSD (2.0 Å), the highest lDDT (0.79),
and a higher TM-score (0.74) than `model_C` (0.62), yet its `chain_mapping_complete`
flag is `false`. It is therefore placed after both mapping-complete models and receives
`decision = reject_incomplete_mapping`, despite being second-best on TM-score. In other
words, an incomplete chain mapping is disqualifying: such a model can never be
`preferred` or `alternate`, regardless of its global-similarity metrics. The two
mapping-complete models, `model_A` and `model_C`, are then separated by the next key,
higher TM-score (0.81 vs 0.62), giving rank 1/`preferred` to `model_A` and
rank 2/`alternate` to `model_C`.

## How critical-region uncertainty affects the ranking

The critical region is residues 181–240. Its mean error (Å) is the 6th key of the
sorting tuple: it can only break ties among models that are identical on the five
stronger keys (mapping completeness, TM-score, lDDT, RMSD, aligned residues). In this
fixture no such tie occurs, so the critical-region risk does not change any rank.

It is still reported because it materially qualifies the preferred choice:

- `model_A` (preferred) carries a critical-region mean error of 4.8 Å — the highest
  local error among its own regions (1.4 Å for residues 1–180) — so its advantage is
  global, while uncertainty in residues 181–240 is comparatively large.
- `model_C` (alternate) has the lowest critical-region error (4.1 Å), making it the
  fallback of choice if downstream use emphasizes the 181–240 region.
- `model_B` combines incomplete mapping with the worst critical-region error (8.5 Å),
  reinforcing its rejection.

## Tie-breaks and secondary diagnostics

No tie-breaking beyond the first two keys was needed: the three models differ in
mapping completeness and TM-score alone. Had ties occurred, the rule's remaining keys
apply in order — higher lDDT, lower RMSD, more aligned residues, lower critical-region
mean error, and finally lexical `model_id`. RMSD and aligned-residue coverage are kept
as secondary diagnostics only and were never mixed into a unitless score.

## Scope limitations

All statements are limited to the metric values supplied in `inputs/`. No
coordinate-level interpretation, interface characterization, or experimental property
claims are made or implied.
