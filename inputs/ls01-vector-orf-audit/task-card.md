# 任务卡： `ls01-vector-orf-audit`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS01-3｜表达载体 ORF 与克隆兼容性审计 — `ls01-vector-orf-audit`

**Formal status:** `ready` — ORF/tag audit rules, construct-level gold and accepted static oracle are frozen.

| 字段 | 内容 |
| --- | --- |
| ID | ls01-vector-orf-audit |
| Domain / sub-domain | molecular biology / construct QC |
| Level / time | L2, 35 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, A, O |
| 来源思想 | custom |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/AUDIT_RULE.md` — 704 bytes
- `inputs/constructs.csv` — 262 bytes

**Total:** 966 bytes (0.00 MiB).

### Prompt（运行时仅复制本框）

> Audit each construct for start/stop codons, reading frame, tag/linker compatibility and cloning flags represented in the input. Write output/construct_audit.csv with construct_id,frame_ok,start_ok,stop_ok,tag_ok,overall_status,issues, output/analysis.py, and output/report.md. Do not infer sequence features that are absent from the input.

### Deliverables / Output contract

output/construct_audit.csv

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ every construct once

□ boolean fields valid

□ issues agree with status

□ rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | Every construct is present once with required audit fields. |
| Core science | 40 | Start/stop, frame, tag/linker and cloning checks match frozen rules. |
| Direction / decision | 15 | Overall status and issue labels agree with component checks. |
| Summary consistency | 5 | Report names usable constructs and blocking defects consistently. |
| Script / reproducibility | 10 | Standard static rerunnable-script checks. |

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

Skill: C0=`NONE`; T0=`AUTO`; T1=`dna-rna-sequence-analysis`, `synthetic_biology_design`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

