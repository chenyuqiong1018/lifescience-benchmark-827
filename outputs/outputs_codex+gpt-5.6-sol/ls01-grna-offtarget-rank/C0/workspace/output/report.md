# Guide ranking report

All six candidate guides were ranked exactly once using only the supplied CSV files and the participant-visible frozen rule. No external data were fetched and the inputs were not modified.

## Result

`g02` is recommended. It clears the annotated rejection filters and has the highest on-target activity among the three non-rejected guides. `g06` and `g04` are acceptable in that order because their activity scores are 0.64 and 0.59 after both clear the safety filters.

`g03` is rejected as high risk because its annotated off-target is a coding-exon near match with two mismatches. Its higher activity (0.74) is an explicit trade-off and does not override the frozen safety rule. `g01` and `g05` are critical-risk rejects: both have a one-mismatch coding-exon hit and high expression risk; they also have one and two coding bystanders, respectively. Their high activity scores therefore do not make them acceptable.

## Deterministic order

1. g02 — recommend, low risk, activity 0.67
2. g06 — acceptable, low risk, activity 0.64
3. g04 — acceptable, low risk, activity 0.59
4. g03 — reject, high risk, activity 0.74
5. g01 — reject, critical risk, activity 0.82
6. g05 — reject, critical risk, activity 0.78
