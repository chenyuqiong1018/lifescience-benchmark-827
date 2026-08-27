# 任务卡： `ls06-eno1-significance-audit`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls06-eno1-significance-audit`

| 字段 | 内容 |
| --- | --- |
| ID | ls06-eno1-significance-audit |
| Domain / sub-domain | proteomics / multiple-testing audit |
| Level / time | L2, 30 min |
| Anchor / related capabilities | [tabular-analysis] |
| 来源思想 | BixBench / bix-37-q3 |

### Inputs (authoritative packaged inventory)
- `inputs/MeRIP_RNA_result.xlsx` — 1,155,180 bytes
- `inputs/Proteomic_data .xlsx` — 646,418 bytes
- `inputs/README.md` — 970 bytes

**Total:** 1,802,568 bytes (1.72 MiB).

### Prompt

> Retrieve ENO1's adjusted p-value from the supplied proteomics results and give a threshold-calibrated interpretation at FDR 0.05. Write `output/eno1_significance.json` with `gene,adjusted_p_value,fdr_threshold,significant,source_file,source_sheet`, `output/analysis.py`, and `output/report.md`. Do not relabel a raw p-value as adjusted.

### Deliverables

one JSON object; report; rerunnable script. `significant` is a JSON boolean.

### Hard gates

exact ENO1/source; adjusted rather than raw p-value; finite p in `[0,1]`; boolean agrees with FDR 0.05.

### DeterministicArtifactScore（0–80）

coverage/schema 10; adjusted p-value 40 (absolute tolerance `0.0005`); FDR decision 15; report consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`proteome_analysis`, `biomarker_discovery`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

