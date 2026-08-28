# Pocket Reliability Assessment

## Verdict

- Pocket reliable: **False**
- Prioritize mutations: **False**

The nominated pocket is **not reliable**. Both pocket intervals have pLDDT < 50 (coordinates should not be interpreted) and PAE > 10 A (uncertain relative placement). All pocket mutation candidates are therefore deferred pending experimental structure validation. Non-pocket candidates are out of scope.

## Method

Each mutation was joined to the residue interval in `confidence.csv` containing its residue number. Pocket candidates were ranked by higher pLDDT, then lower PAE, then lower absolute predicted ddG, then lexical mutation; non-pocket candidates follow. Reliability requires every pocket interval to satisfy pLDDT >= 70 and PAE <= 10 A.

## Ranked mutation priorities

| rank | mutation | region | pLDDT | PAE (A) | predicted ddG | pocket_support | confidence_penalty | decision |
|------|----------|--------|-------|---------|---------------|----------------|--------------------|----------|
| 1 | Y218F | pocket | 47 | 14.2 | 0.4 | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 2 | R244A | pocket | 43 | 17.8 | 2.8 | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 3 | L265W | pocket | 43 | 17.8 | 4.1 | unsupported_low_confidence | plddt_lt_50;pae_gt_10A | defer_structure_validation |
| 4 | S301A | core | 84 | 3.4 | 0.2 | cautious_support | none | out_of_scope_non_pocket |

## Confidence intervals

| residues | pLDDT | PAE to core (A) | pocket member |
|----------|-------|-----------------|---------------|
| 1-209 | 88 | 2.1 | false |
| 210-230 | 47 | 14.2 | true |
| 231-270 | 43 | 17.8 | true |
| 271-340 | 84 | 3.4 | false |

## Uncertainty and caveats

- pLDDT and PAE are per-model confidence estimates; low pLDDT (<50) and high PAE (>10 A) together mean the pocket geometry and its relative placement are too uncertain to trust for mutation design.
- All predicted ddG and activity-change values are **computational model hypotheses**. They are not measured binding affinities or enzymatic activities and must not be reported as such.
- No experimental binding or activity measurements are present in the inputs.

## References

- AlphaFold DB confidence guidance (pLDDT / PAE): https://alphafold.ebi.ac.uk/faq
- Terwilliger et al. (2024). AlphaFold predictions are hypotheses and do not replace experimental structure determination. Nat Methods. DOI 10.1038/s41592-023-02087-4
