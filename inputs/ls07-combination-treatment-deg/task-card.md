# 任务卡： `ls07-combination-treatment-deg`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls07-combination-treatment-deg`

| 字段 | 内容 |
| --- | --- |
| ID | ls07-combination-treatment-deg |
| Domain / sub-domain | transcriptomics / differential expression |
| Level / time | L3, 75 min |
| Anchor / related capabilities | [transcriptome-analysis] |
| 来源思想 | BixBench / bix-43-q3 |

### Inputs (authoritative packaged inventory)
- `inputs/README.md` — 1,193 bytes
- `inputs/counts_raw_unfiltered.csv` — 6,016,681 bytes
- `inputs/ensg_to_gene_name.tsv` — 2,360,408 bytes
- `inputs/sample_layout.csv` — 2,860 bytes

**Total:** 8,381,142 bytes (7.99 MiB).

### Prompt

> Use the files in `inputs/` to perform differential-expression analysis for `Cisplatin_IC50_CBD_IC50` versus `DMSO`. Use only the three replicates from each of those groups, with combination treatment as numerator and DMSO as denominator, and use `Group` as the only design term. Before fitting, retain genes for which at least one of the six selected samples has a raw count greater than 10. Use PyDESeq2 0.5.0 with `refit_cooks=True` and the standard `DeseqStats` contrast. A gene passes when `padj < 0.05`, `abs(log2FoldChange) > 0.5`, and `baseMean > 10`.
>
> Write all results under `output/`: `differential_expression.csv`, `summary.json`, `analysis.py`, and `report.md`. Do not use samples from other groups. Preserve unavailable adjusted p-values as null rather than converting them to zero. The report must be no more than 500 words and distinguish statistical association from causation.
### Deliverables

unique gene rows; JSON count and explicit contrast/design metadata; report; rerunnable script. CSV missing numeric values are empty.

### Hard gates

design and contrast recorded exactly; gene IDs unique; threshold rule uses strict inequalities exactly; summary count equals passing rows.

### DeterministicArtifactScore（0–80）

coverage/schema 10; frozen reference values and pass count 40 using the pinned DESeq2 environment/tolerances; direction/threshold decision 15; summary consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.
- Readiness closure: the prompt-authoritative PyDESeq2 0.5.0 analysis, full-row reference, 555-gene threshold result and 677/679 upstream discrepancy adjudication are frozen in the task-local oracle.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`biomarker_discovery`, `transcriptome_analysis`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

