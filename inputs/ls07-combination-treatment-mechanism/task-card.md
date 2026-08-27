# 任务卡： `ls07-combination-treatment-mechanism`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls07-combination-treatment-mechanism`

| 字段 | 内容 |
| --- | --- |
| ID | ls07-combination-treatment-mechanism |
| Domain / sub-domain | transcriptomics / pathway enrichment |
| Level / time | L3, 90 min |
| Anchor / related capabilities | [pathway-analysis] |
| 来源思想 | BixBench / bix-43-q5 |

### Inputs (authoritative packaged inventory)
- `inputs/README.md` — 1,692 bytes
- `inputs/Reactome_2022.background.txt` — 65,942 bytes
- `inputs/Reactome_2022.gmt` — 778,913 bytes
- `inputs/Reactome_2022.manifest.json` — 1,515 bytes
- `inputs/counts_raw_unfiltered.csv` — 6,016,681 bytes
- `inputs/ensg_to_gene_name.tsv` — 2,360,408 bytes
- `inputs/sample_layout.csv` — 2,860 bytes

**Total:** 9,228,011 bytes (8.80 MiB).

### Prompt

> Using the files in `inputs/`, perform the approved `Cisplatin_IC50_CBD_IC50` versus `DMSO` differential-expression analysis and enrichment with GSEApy 1.1.4 against the evaluator-supplied frozen `Reactome_2022` resource and supplied background universe. Identify the best-supported primary cellular mechanism. Do not download or substitute a current pathway library or identifier mapping.
>
> Write all results under `output/`: `pathway_enrichment.csv`, `mechanism_call.json`, `resource_manifest.json`, `analysis.py`, and `report.md`. The report must be no more than 500 words and must distinguish pathway enrichment from demonstrated causation.
### Deliverables

pathway table with declared tested universe/release; mechanism JSON referencing a table row; report; rerunnable script. Missing statistics are empty/null.

### Hard gates

exact pinned gene-set release and universe used; corrected enrichment statistics valid; mechanism call supported by a reported row; no causal overclaim.

### DeterministicArtifactScore（0–80）

coverage/schema 10; overlap/statistics/ranking and primary mechanism 40 against the pinned reference; evidence direction/restraint 15; summary consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.
- Readiness closure: the Reactome library/background, DE mapping, complete 1,818-pathway enrichment table and top-mechanism call are frozen in the task-local oracle; local results are explicitly separated from historical live-Enrichr anchors.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`go_term_analysis`, `string-ppi-enrichment`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

