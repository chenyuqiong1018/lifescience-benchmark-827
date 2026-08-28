# Structural model ranking

## Result

| Rank | Model | Global score | Interface score | Critical-residue risk (Å) | Decision |
|---:|---|---:|---:|---:|---|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

## Interpretation

The ranking applies the supplied stable tuple exactly. Complete chain mapping is the first discriminator, so model_A and model_C rank ahead of model_B even though model_B has a higher TM-score than model_C. The interface score is only the prescribed mapping-completeness indicator: 1 for a complete mapping and 0 for an incomplete mapping. It is not a coordinate-derived interface-quality claim.

Within the mapping-complete group, model_A ranks first because its TM-score (0.81) exceeds model_C's (0.62). Model_C has lower critical-region uncertainty (4.1 Å versus 4.8 Å for model_A), but critical-region risk is the sixth tie-breaker and therefore does not override the earlier TM-score comparison. Model_B has both incomplete mapping and the highest critical-region uncertainty (8.5 Å), so its decision is `reject_incomplete_mapping`.

Critical-region risk is the length-weighted mean error over residues 181–240. Each supplied model has one row spanning that full interval, so the result equals that row's mean error. These values are already-computed comparison metrics; no coordinate-level, physical-interface, or experimental property is inferred.
