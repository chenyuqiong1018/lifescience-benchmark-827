Done. All deliverables are in `output/`, inputs untouched, no external data used.

## Final ranking (`output/ranked_guides.csv`)

| rank | guide | score | risk | decision | why |
|---|---|---|---|---|---|
| 1 | g02 | 0.67 | low | **recommend** | Best activity among safety-passing guides; only off-target is intronic, 3 mismatches, low expression risk, 0 bystanders |
| 2 | g06 | 0.64 | low | acceptable | Intergenic off-target, 3 mm, low risk |
| 3 | g04 | 0.59 | low | acceptable | Intergenic off-target, 4 mm, low risk — weakest usable activity |
| 4 | g03 | 0.74 | high | reject | Coding-exon off-target at 2 mismatches |
| 5 | g01 | 0.82 | critical | reject | Coding-exon at **1 mm + high expression risk** + 1 bystander |
| 6 | g05 | 0.78 | critical | reject | Coding-exon at 1 mm + high expression risk + 2 bystanders |

## Key points

- **Followed the frozen rule in `inputs/RANKING_RULE.md`** exactly: rejection gates (coding exon ≤2 mm / high expression risk / bystanders >0), risk classes, non-rejected-first ordering, and within-rejected `high` before `critical`.
- **Trade-off stated, not hidden:** the two most active guides (g01 0.82, g05 0.78) and g03 (0.74, which beats every usable guide) are all rejected because coding/exonic near matches with low mismatch counts are decisive safety evidence — the ranking sacrifices up to ~0.23 on-target activity to avoid likely damaging coding edits. This is spelled out in each row's rationale and in §5 of the report.
- **Deliverables:** `output/ranked_guides.csv` (exact requested header), `output/analysis.py` (self-contained, deterministic — re-run reproduces the CSV byte-for-byte, hash-verified), `output/report.md` (evidence table, policy, trade-offs, validation notes).
- **Validated invariants:** unique sequential ranks 1–6, every input guide appears exactly once, scores match inputs, rank 1 = `recommend`, all rejected guides labeled `reject`.