# 任务卡： `ls10-neun-power-analysis`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls10-neun-power-analysis`

| 字段 | 内容 |
| --- | --- |
| ID | ls10-neun-power-analysis |
| Domain / sub-domain | biostatistics / power analysis |
| Level / time | L2, 40 min |
| Anchor / related capabilities | [biostatistics] |
| 来源思想 | BixBench / bix-19-q1/q2 |

### Inputs (authoritative packaged inventory)
- `inputs/NeuN_quantification.csv` — 218 bytes
- `inputs/README.md` — 770 bytes

**Total:** 988 bytes (0.00 MiB).

### Prompt

> Estimate the standardized mean difference (Cohen's d) between the two supplied groups and the required equal sample size per group for a two-sided independent t-test at alpha 0.05 and power 0.80. Write `output/power_result.json` with `group_labels,n_each,means,sds,pooled_sd,cohens_d,alpha,power,alternative,required_n_per_group`, `output/analysis.py`, and `output/report.md`. Round required n upward.

### Deliverables

one JSON object with group-keyed or label-aligned arrays; report; rerunnable script. Sample SD convention and signed-d order must be stated.

### Hard gates

both groups mapped correctly; finite means/SD/pooled SD/effect size; two-sided alpha/power specification exact; sample size rounded upward.

### DeterministicArtifactScore（0–80）

coverage/schema 10; means 8, SDs 8, pooled SD 6, absolute Cohen d 8 (`5e-3` tolerance), required n/group 10; specification/direction 15; report consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`biomarker_discovery`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

