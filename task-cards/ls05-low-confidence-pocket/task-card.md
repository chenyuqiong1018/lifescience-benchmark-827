# 任务卡： `ls05-low-confidence-pocket`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS05-3｜低置信度口袋与突变优先级评估 — `ls05-low-confidence-pocket`

**Formal status:** `ready_local_extension` — Frozen pLDDT/PAE rule, gold, accepted oracle and 3/3 acceptance tests present; synthetic fixture, not upstream benchmark.

| 字段 | 内容 |
| --- | --- |
| ID | ls05-low-confidence-pocket |
| Domain / sub-domain | structural biology / pocket uncertainty |
| Level / time | L2, 35 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, A, O |
| 来源思想 | benchmark-informed local extension (not an upstream benchmark item) |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/SCORING_RULE.md` — 1,753 bytes
- `inputs/confidence.csv` — 142 bytes
- `inputs/mutation_candidates.csv` — 167 bytes

**Total:** 2,062 bytes (0.00 MiB).

### Prompt（运行时仅复制本框）

> Using only the files in inputs/, assess whether the nominated pocket is reliable enough to prioritize mutations, following inputs/SCORING_RULE.md exactly. Write output/mutation_priorities.csv with rank,mutation,pocket_support,confidence_penalty,decision, output/pocket_assessment.json, output/analysis.py, and output/report.md. Propagate pLDDT/PAE uncertainty and do not describe predicted ΔΔG/activity as measured effects.

### Deliverables / Output contract

output/mutation_priorities.csv

output/pocket_assessment.json

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ all candidates once

□ ranks unique

□ confidence evidence traceable

□ rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | 10: every mutation occurs once, ranks are unique, and all required artifacts parse. |
| Core science | 40 | 40: exact order (10), eight pocket-support/penalty fields prorated (20), and false/false calibrated pocket assessment with reason (10). |
| Direction / decision | 15 | 15: four defer/out-of-scope decisions prorated against frozen gold. |
| Summary consistency | 5 | 5: report states low confidence and deferral without prohibited experimental overclaim. |
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

Skill: C0=`NONE`; T0=`AUTO`; T1=`protein_quality_assessment`, `binding_site_characterization`, `alphafold_structure_pipeline`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

