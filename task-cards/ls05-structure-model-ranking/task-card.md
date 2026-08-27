# 任务卡： `ls05-structure-model-ranking`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS05-2｜结构模型置信度综合排名 — `ls05-structure-model-ranking`

**Formal status:** `ready_local_extension` — Frozen rule, gold, accepted oracle and 3/3 acceptance tests present; synthetic fixture, not upstream benchmark.

| 字段 | 内容 |
| --- | --- |
| ID | ls05-structure-model-ranking |
| Domain / sub-domain | structural biology / model confidence |
| Level / time | L2, 35 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, A, O |
| 来源思想 | benchmark-informed local extension (not an upstream benchmark item) |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/SCORING_RULE.md` — 1,518 bytes
- `inputs/model_metrics.csv` — 164 bytes
- `inputs/residue_errors.csv` — 162 bytes

**Total:** 1,844 bytes (0.00 MiB).

### Prompt（运行时仅复制本框）

> Using only the files in inputs/, rank every supplied structural model exactly according to inputs/SCORING_RULE.md. Write output/model_ranking.csv with rank,model_id,global_score,interface_score,critical_residue_risk,decision, output/analysis.py, and output/report.md. Explain how chain-mapping completeness and critical-region uncertainty affect the ranking. Do not claim coordinate-level, interface, or experimental properties that are not present in the inputs.

### Deliverables / Output contract

output/model_ranking.csv

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ every model once

□ unique ranks

□ input metrics preserved

□ rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | 10: every model occurs once, ranks are unique, and all required artifacts parse. |
| Core science | 40 | 40: exact rank tuple/order (20) plus nine frozen global/interface/critical-risk fields prorated (20). |
| Direction / decision | 15 | 15: three preferred/alternate/reject decisions prorated against frozen gold. |
| Summary consistency | 5 | 5: report names model_A and identifies model_B's incomplete mapping consistently. |
| Script / reproducibility | 10 | 10: analysis.py parses and contains no absolute user path; submission code is never imported by oracle. |

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

Skill: C0=`NONE`; T0=`AUTO`; T1=`protein_quality_assessment`, `protein_structure_analysis`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

