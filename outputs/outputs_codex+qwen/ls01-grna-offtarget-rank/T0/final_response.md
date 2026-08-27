All deliverables are complete in `./output`, built solely from `./inputs` (inputs untouched, no external data).

## Results

I implemented the frozen deterministic policy from `inputs/RANKING_RULE.md` verbatim: reject on any coding-exon off-target ≤2 mismatches, any `expression_risk=high` hit, or `coding_bystander_count>0`; classify risk as critical/high/low; order kept guides by activity, rejected guides with `high` before `critical`.

| Rank | Guide | On-target | Risk | Decision | Why |
|---|---|---|---|---|---|
| 1 | g02 | 0.67 | low | **recommend** | intronic hit, 3 MM, low expr risk, 0 bystanders |
| 2 | g06 | 0.64 | low | acceptable | intergenic, 3 MM, low expr risk |
| 3 | g04 | 0.59 | low | acceptable | intergenic, 4 MM, low expr risk |
| 4 | g03 | 0.74 | high | reject | coding-exon hit with 2 MM |
| 5 | g01 | 0.82 | critical | reject | coding-exon 1 MM + high expr risk + 1 bystander |
| 6 | g05 | 0.78 | critical | reject | coding-exon 1 MM + high expr risk + 2 bystanders |

**Key trade-off, stated explicitly:** the two most active guides (g01 = 0.82, g05 = 0.78) are rejected because 1-mismatch coding-exon near matches in highly expressed loci are decisive safety evidence that activity cannot offset. The recommendation (g02) deliberately sacrifices ~0.15 activity versus g01 to eliminate all coding-exon/high-expression liability. Coding/exonic near matches and mismatch counts are treated as primary safety evidence in both the decisions and every per-guide rationale.

## Deliverables
- `output/ranked_guides.csv` — required columns; ranks 1–6 unique; all 6 guides appear exactly once (validated programmatically)
- `output/analysis.py` — reproducible implementation with built-in assertions (unique ranks, full guide coverage), re-runs idempotently
- `output/report.md` — method, per-guide assessment, explicit trade-off summary, verification checklist