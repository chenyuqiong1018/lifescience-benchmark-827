# Low-confidence pocket assessment

## Answer

The nominated pocket is **not reliable enough to prioritize mutations**. Both pocket intervals fail the pLDDT and PAE requirements, so all three pocket candidates are deferred for structure validation.

## Candidate-level uncertainty propagation

| Rank | Mutation | Region | Joined interval | pLDDT | PAE (Å) | pLDDT margin to 70 | PAE margin to ≤10 Å | Predicted ΔΔG hypothesis (kcal/mol) | Activity hypothesis | Decision |
|---:|---|---|---|---:|---:|---:|---:|---:|---|---|
| 1 | Y218F | pocket | 210–230 | 47 | 14.2 | -23 | -4.2 | 0.4 | increase | defer_structure_validation |
| 2 | R244A | pocket | 231–270 | 43 | 17.8 | -27 | -7.8 | 2.8 | increase | defer_structure_validation |
| 3 | L265W | pocket | 231–270 | 43 | 17.8 | -27 | -7.8 | 4.1 | unknown | defer_structure_validation |
| 4 | S301A | core | 271–340 | 84 | 3.4 | +14 | +6.6 | 0.2 | no_change | out_of_scope_non_pocket |

Negative margins indicate failure. Y218F is 23 pLDDT points below and 4.2 Å above the reliability limits. R244A and L265W are each 27 pLDDT points below and 7.8 Å above the limits. These are not borderline failures: every pocket candidate inherits two confidence penalties.

## Deterministic ranking trace

1. All pocket candidates precede the non-pocket S301A.
2. Y218F precedes R244A and L265W because its joined interval has higher pLDDT (47 versus 43), then lower PAE (14.2 versus 17.8 Å).
3. R244A and L265W tie on pLDDT and PAE; the lower absolute predicted ΔΔG hypothesis ranks R244A before L265W (2.8 versus 4.1 kcal/mol).

## Evidence-quality assessment

**What is supported:** the supplied interval join, pLDDT/PAE thresholds, deterministic ordering, and conclusion that the pocket fails the prescribed reliability rule.

**What is not available:** no PDB coordinates, UniProt identifier, sequence, ligand, pocket-prediction result, or atom-level quality data were supplied. Therefore AlphaFold download, fpocket/P2Rank, structural quality, visualization, and binding-site characterization tools were not run; manufacturing such results would exceed the evidence.

**Claim boundary:** predicted ΔΔG and predicted activity are hypotheses, not measured stability, binding, or activity effects. Their values cannot overcome the low-confidence structure. The proportionate recommendation is to validate the pocket structure or obtain a better-supported model before prioritizing mutagenesis, followed by experiments to establish actual effects.
