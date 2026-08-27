# Structural model ranking

## Final ranking

| Rank | Model | Global score | Interface score | Critical-residue risk (Å) | Decision |
|---:|---|---:|---:|---:|---|
| 1 | model_A | 0.81 | 1 | 4.8 | preferred |
| 2 | model_C | 0.62 | 1 | 4.1 | alternate |
| 3 | model_B | 0.74 | 0 | 8.5 | reject_incomplete_mapping |

## Audited input metrics

| Model | Chain mapping | TM-score | lDDT | RMSD (Å) | Aligned residues | Critical risk (Å) | Critical overlap (residues) |
|---|---|---:|---:|---:|---:|---:|---:|
| model_A | complete | 0.81 | 0.76 | 2.4 | 238 | 4.8 | 60 |
| model_C | complete | 0.62 | 0.65 | 3.8 | 240 | 4.1 | 60 |
| model_B | incomplete | 0.74 | 0.79 | 2.0 | 211 | 8.5 | 60 |

Both input schemas match exactly, all fields are populated, model identifiers are unique and consistent across files, and each model contributes 60 residue positions to the length-weighted mean over residues 181–240.

## Decision trace

1. **model_A versus model_C:** both have complete chain mappings, so the first tuple element ties. TM-score is the first differing element; 0.81 exceeds 0.62, making model_A preferred. model_C's lower critical-region risk (4.1 Å versus 4.8 Å) occurs only at the sixth discriminator and cannot reverse that decision.
2. **model_C versus model_B:** chain-mapping completeness is the first discriminator. model_C is therefore ahead despite model_B's higher TM-score (0.74 versus 0.62), higher lDDT, and lower RMSD. Incomplete mapping also fixes model_B's interface score at 0 and its decision at `reject_incomplete_mapping`.

Critical-region uncertainty is highest for model_B (8.5 Å), then model_A (4.8 Å), then model_C (4.1 Å). These risks can resolve a tie only after mapping completeness, TM-score, lDDT, RMSD, and aligned-residue count all tie.

## Evidence boundary

The supplied files contain comparison metrics but no coordinates or PDB identifier. Accordingly, no coordinate geometry, atom-level quality, physical interface, composition, visualization, or experimental property was calculated or claimed. `interface_score` is strictly the rule-prescribed binary chain-mapping indicator. With a deterministic frozen tuple and only three candidates, inferential hypothesis tests, p-values, effect sizes, and fitted composite scores would be unsupported and were not introduced.
