# 任务卡： `ls01-grna-offtarget-rank`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS01-1｜CRISPR guide 活性与脱靶风险排名 — `ls01-grna-offtarget-rank`

**Formal status:** `ready` — Ranking rule, task-specific gold and accepted static oracle are frozen.

| 字段 | 内容 |
| --- | --- |
| ID | ls01-grna-offtarget-rank |
| Domain / sub-domain | molecular biology / CRISPR design |
| Level / time | L2, 40 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, A, O |
| 来源思想 | source: custom fixture |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/RANKING_RULE.md` — 970 bytes
- `inputs/candidates.csv` — 246 bytes
- `inputs/off_targets.csv` — 295 bytes

**Total:** 1,511 bytes (0.00 MiB).

### Prompt（运行时仅复制本框）

> Read the two supplied CSV files and rank all candidate guides by on-target activity and annotated off-target risk. Do not fetch external data or alter inputs. Write output/ranked_guides.csv with rank,guide_id,on_target_score,risk_class,decision,rationale, output/analysis.py, and output/report.md. Rankings must be unique; every input guide must appear once; explicitly treat coding/exonic near matches and mismatch count as safety evidence and state any trade-off rather than hiding it.

### Deliverables / Output contract

output/ranked_guides.csv

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ all guides exactly once

□ unique integer ranks 1..N

□ numeric values traceable to inputs

□ decisions/rationales nonempty

□ script reruns cleanly


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | All input guide IDs appear exactly once; required columns parse. |
| Core science | 40 | Guide ranks reproduce frozen activity/off-target risk rule, including coding near matches and mismatch weighting. |
| Direction / decision | 15 | Recommended/reject decisions agree with rank and frozen safety thresholds. |
| Summary consistency | 5 | Report states top choice and essential trade-off consistently with table. |
| Script / reproducibility | 10 | analysis.py exists, parses statically, uses relative inputs/output paths and recreates declared machine-readable artifacts. |

No scientific points may be emitted until a static oracle, tolerances and correct/empty/wrong controls are independently accepted 3/3. Missing core artifact gives zero deterministic points.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`dna-rna-sequence-analysis`, `genome_annotation`, `transcriptome_analysis`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

