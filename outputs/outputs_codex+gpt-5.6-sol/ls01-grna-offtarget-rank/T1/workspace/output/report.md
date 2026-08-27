# Guide ranking report

This T1 run used the controlled skill set, freshly installed and opened in the specified order: `dna-rna-sequence-analysis`, `genome_annotation`, and `code_execution_analysis`. The skills informed local sequence validation, interpretation of coding-exon annotations, and reproducible computation. No remote tools or external data were used because the prompt expressly prohibits them.

Input validation confirmed six unique 20-nt DNA spacers and complete guide-ID coverage in the off-target annotations.

## Ranked result

1. **g02 — recommend; low risk; 0.67 activity.** Three-mismatch intronic hit; best activity after safety clearance.
2. **g06 — acceptable; low risk; 0.64 activity.** Three-mismatch intergenic hit.
3. **g04 — acceptable; low risk; 0.59 activity.** Four-mismatch intergenic hit.
4. **g03 — reject; high risk; 0.74 activity.** Two-mismatch coding-exon hit; activity cannot override the safety rule.
5. **g01 — reject; critical risk; 0.82 activity.** One-mismatch coding-exon hit, high expression risk, and one coding bystander; activity cannot override rejection.
6. **g05 — reject; critical risk; 0.78 activity.** One-mismatch coding-exon hit, high expression risk, and two coding bystanders; activity cannot override rejection.

The ranking is deterministic and unique. Every candidate appears exactly once.
