# 任务卡： `ls06-eno1-effect-size`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls06-eno1-effect-size`

| 字段 | 内容 |
| --- | --- |
| ID | ls06-eno1-effect-size |
| Domain / sub-domain | proteomics / effect-size calculation |
| Level / time | L2, 35 min |
| Anchor / related capabilities | [tabular-analysis] |
| 来源思想 | BixBench / bix-37-q1/q4 |

### Inputs (authoritative packaged inventory)
- `inputs/MeRIP_RNA_result.xlsx` — 1,155,180 bytes
- `inputs/Proteomic_data .xlsx` — 646,418 bytes
- `inputs/README.md` — 1,248 bytes

**Total:** 1,802,846 bytes (1.72 MiB).

### Prompt

> Using the supplied proteomics results, calculate ENO1 tumor-versus-normal fold change and log2 fold change. Write `output/eno1_effect.json` with `gene,tumor_value,normal_value,fold_change,log2_fold_change,source_file,source_sheet`, `output/analysis.py`, and `output/report.md`. State the fold-change direction and do not substitute the unrelated workbook.

### Deliverables

one JSON object with finite numeric values and source identifiers; UTF-8 Markdown report; rerunnable Python script. No additional artifact is required.

### Hard gates

exact ENO1/source row; all four core values within checker tolerance; fold-change/log2 direction internally consistent; source file and sheet traceable.

### DeterministicArtifactScore（0–80）

coverage/schema 10; Normal, Tumor, fold change and log2 fold change 10 each (raw values relative tolerance `5e-6`, fold change `2e-3`, log2 absolute tolerance `0.011`); direction 15; report consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`proteome_analysis`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

