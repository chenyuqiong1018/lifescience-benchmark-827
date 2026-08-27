# 任务卡： `ls10-treatment-response-model`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls10-treatment-response-model`

| 字段 | 内容 |
| --- | --- |
| ID | ls10-treatment-response-model |
| Domain / sub-domain | biostatistics / logistic regression |
| Level / time | L2, 45 min |
| Anchor / related capabilities | [biostatistics] |
| 来源思想 | BixBench / bix-51-q3/q4 |

### Inputs (authoritative packaged inventory)
- `inputs/README.md` — 806 bytes
- `inputs/data.xlsx` — 22,788 bytes

**Total:** 23,594 bytes (0.02 MiB).

### Prompt

> Fit a logistic regression for the binary treatment-response outcome using BMI, age and gender. Use complete cases, document outcome coding and gender reference level, and report the age log-odds coefficient and two-sided p-value. Write `output/model_coefficients.csv` with `term,estimate,std_error,z,p_value,odds_ratio`, `output/model_metadata.json`, `output/analysis.py`, and `output/report.md`.

### Deliverables

unique coefficient rows; metadata with formula, outcome coding, reference level, complete-case count and implementation/version; report; rerunnable script.

### Hard gates

specified model only; binary outcome coding and gender reference documented; age term unique and finite; coefficient, odds ratio and significance interpretation mutually consistent.

### DeterministicArtifactScore（0–80）

coverage/schema 10; age estimate, SE, two-sided p-value and odds ratio 10 each (`rel_tol=3e-3`, `abs_tol=5e-5`); age direction/decision 15; report consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.

## Release gate

A card enters the main result only if its reference submission passes 3/3 clean reruns, empty output and at least one format-correct scientific error fail 3/3, one domain reviewer and one grader reviewer accept it, and a timed calibration run can be frozen and rescored. On 2026-08-17, all ten cards have accepted static checkers; campaign-level reviewer and platform deviations remain governed by `docs/contracts/formal-eval-release-status-2026-08-17.md`.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`personalized_medicine`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

