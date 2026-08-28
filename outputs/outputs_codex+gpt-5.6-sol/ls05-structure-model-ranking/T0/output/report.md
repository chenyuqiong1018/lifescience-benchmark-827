# Structural model ranking

## Frozen-rule result

| Rank | Model | Global score | Interface score | Critical-residue risk (Å) | Decision |
|---:|---|---:|---:|---:|---|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

## Input and rule audit

| Model | Chain mapping | TM-score | lDDT | RMSD (Å) | Aligned residues | Residues 181–240 risk (Å) |
|---|---|---:|---:|---:|---:|---:|
| model_A | complete | 0.81 | 0.76 | 2.4 | 238 | 4.8 |
| model_C | complete | 0.62 | 0.65 | 3.8 | 240 | 4.1 |
| model_B | incomplete | 0.74 | 0.79 | 2.0 | 211 | 8.5 |

The analysis found all required columns, no missing values, unique model identifiers, and exactly one error value covering every residue in the critical interval 181–240 for each model. Because the supplied rule is deterministic, no hypothesis test, p-value, or fitted composite score is warranted for these three models.

## Ranking interpretation

Chain-mapping completeness is the first and therefore decisive discriminator. It places model_A and model_C ahead of model_B, even though model_B's TM-score (0.74) is higher than model_C's (0.62). The prescribed `interface_score` records only this completeness state; it is not evidence of physical interface quality.

Among complete mappings, model_A is preferred because its higher TM-score (0.81 versus 0.62) is encountered before all remaining tie-breakers. Critical-region uncertainty is evaluated over residues 181–240. Model_C's risk is lower than model_A's (4.1 versus 4.8 Å), but risk is sixth in the tuple and cannot reverse the earlier TM-score ordering. Model_B combines incomplete mapping with the highest risk (8.5 Å) and is rejected for incomplete mapping.

Only supplied comparison metrics are interpreted. No coordinate-level geometry, physical interface, or experimental property is claimed.
