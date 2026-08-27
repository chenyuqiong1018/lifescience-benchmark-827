All deliverables are complete and verified. Here's the final ranking produced under the frozen rule in `inputs/RANKING_RULE.md`:

| Rank | Guide | On-target score | Risk class | Decision | Why |
|-----:|-------|----:|----------|------------|-----|
| 1 | g02 | 0.67 | low | **recommend** | Intron hit, 3 mismatches, low expression risk, 0 coding bystanders |
| 2 | g06 | 0.64 | low | acceptable | Intergenic hit, 3 mm, low expression risk |
| 3 | g04 | 0.59 | low | acceptable | Intergenic hit, 4 mm, low expression risk |
| 4 | g03 | 0.74 | high | reject | Coding-exon hit with 2 mismatches (≤2 threshold) |
| 5 | g01 | 0.82 | critical | reject | Coding-exon hit with 1 mm + high expression risk + 1 bystander |
| 6 | g05 | 0.78 | critical | reject | Coding-exon hit with 1 mm + high expression risk + 2 bystanders |

**Deliverables (all in `./output/`):**
- **`ranked_guides.csv`** — 6 rows, unique ranks 1–6, every input guide exactly once, columns `rank,guide_id,on_target_score,risk_class,decision,rationale`
- **`analysis.py`** — deterministic script that reads only `./inputs/*.csv`, applies the frozen rule (with built-in assertions for unique ranks / full guide coverage), and regenerates the CSV
- **`report.md`** — full method, per-guide rationale, and explicit trade-off discussion

**Key trade-off (stated, not hidden):** the recommended guide g02 (0.67) is only the 4th most active of six. The three more active guides — g01 (0.82, the best in the set), g05 (0.78), g03 (0.74) — are all rejected because their off-targets are coding-exonic near-matches (1–2 mismatches), with high expression risk and coding bystanders for g01/g05. Per the frozen rule, coding/exonic near-matches and mismatch count act as hard safety evidence that vetoes activity: ~0.15 on-target activity is deliberately sacrificed (0.82 → 0.67) to avoid near-perfect matches in coding sequence. Among rejected guides, `high` (g03) ranks before `critical` (g01, g05), and g01 precedes g05 on higher activity.

No external data was fetched and inputs were not modified; work stayed entirely within the workspace.