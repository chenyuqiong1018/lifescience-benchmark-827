# Guide ranking report

The two supplied CSV files were analyzed locally. The selected skills informed sequence validation, genome-annotation interpretation, and reproducible code structure. Their optional remote tools were not called because the task prohibits external data.

All six spacers are valid 20-nt DNA sequences, candidate IDs are unique, and the candidate and off-target annotation guide sets match exactly.

## Ranking

1. **g02 — recommend (low risk, activity 0.67).** Its nearest annotated hit has three mismatches and is intronic. It clears the safety filters and has the highest activity among eligible guides.
2. **g06 — acceptable (low risk, activity 0.64).** Its nearest annotated hit has three mismatches and is intergenic.
3. **g04 — acceptable (low risk, activity 0.59).** Its nearest annotated hit has four mismatches and is intergenic.
4. **g03 — reject (high risk, activity 0.74).** It has a two-mismatch coding-exon near match. The higher activity does not override the explicit safety rejection.
5. **g01 — reject (critical risk, activity 0.82).** It has a one-mismatch coding-exon hit, high expression risk, and one coding bystander. These safety findings outweigh its high activity.
6. **g05 — reject (critical risk, activity 0.78).** It has a one-mismatch coding-exon hit, high expression risk, and two coding bystanders. These safety findings outweigh its high activity.

The ordering is unique and follows the frozen rule: eligible guides first by decreasing activity; rejected high-risk guides before critical-risk guides; then decreasing activity and guide ID.
