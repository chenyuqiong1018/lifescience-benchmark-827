# 任务卡： `ls04-perturbseq-reference-map`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS04-2｜Perturb-seq 查询—参考映射 — `ls04-perturbseq-reference-map`

**Formal status:** `ready` — Feature alignment, pseudobulk transform, Hungarian assignment and mapping gold are frozen.

| 字段 | 内容 |
| --- | --- |
| ID | ls04-perturbseq-reference-map |
| Domain / sub-domain | single-cell / Perturb-seq mapping |
| Level / time | L3, 90 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, A, O |
| 来源思想 | CompBioBench perturb-seq-align-q1 |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/perturb.seq.align.q1.query.h5ad` — 19,163,116 bytes
- `inputs/perturb.seq.align.q1.ref.h5ad` — 41,142,620 bytes

**Total:** 60,305,736 bytes (57.51 MiB).

### Prompt（运行时仅复制本框）

> Map query perturbation groups to the labeled reference across the cell-type shift and identify the query guide IDs corresponding to PABPC1, NUDT21 and LEO1. Write output/guide_mapping.csv with target_gene,query_guide_id,score,runner_up_score,confidence, output/analysis.py, and output/report.md. Prevent target metadata leakage and quantify ambiguity.

### Deliverables / Output contract

output/guide_mapping.csv

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ exactly three target genes

□ unique guide call per target

□ finite scores

□ rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | Exactly PABPC1, NUDT21 and LEO1 occur once with unique query guide calls. |
| Core science | 40 | Three guide identities and mapping scores/ranking match frozen leak-free reference analysis. |
| Direction / decision | 15 | Confidence/ambiguity decisions agree with best and runner-up scores. |
| Summary consistency | 5 | Report repeats all three mappings consistently. |
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

Skill: C0=`NONE`; T0=`AUTO`; T1=`scvi-tools`, `scgpt`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

