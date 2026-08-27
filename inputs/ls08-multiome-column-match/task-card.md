# 任务卡： `ls08-multiome-column-match`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls08-multiome-column-match`

| 字段 | 内容 |
| --- | --- |
| ID | ls08-multiome-column-match |
| Domain / sub-domain | single-cell multi-omics / cross-modality matching |
| Level / time | L3, 75 min |
| Anchor / related capabilities | [multi-omics-integration] |
| 来源思想 | CompBioBench / multiome-match-atac-rna-q1 |

### Inputs (authoritative packaged inventory)
- `inputs/MATCHING_RULE.md` — 880 bytes
- `inputs/README.md` — 1,122 bytes
- `inputs/ensembl112_gene_coordinates.tsv` — 4,200,908 bytes
- `inputs/multiome.match.atac.rna.q1.atac.tsv.gz` — 15,259,355 bytes
- `inputs/multiome.match.atac.rna.q1.rna.tsv.gz` — 1,432,352 bytes

**Total:** 20,894,617 bytes (19.93 MiB).

### Prompt

> Recover the one-to-one matching between the eight permuted ATAC population columns and RNA populations. Write `output/column_mapping.csv` with `rna_population,atac_column,match_score,runner_up_score`, `output/score_matrix.csv`, `output/analysis.py`, and `output/report.md`. Enforce a bijection and explain the shared biological signal used.

### Deliverables

eight unique mapping rows; complete finite score matrix; report; rerunnable script. Score definition and preprocessing must be stated.

### Hard gates

all eight labels on each side appear exactly once; mapping is a bijection; all reported scores finite; mapping direction is RNA-to-ATAC.

### DeterministicArtifactScore（0–80）

coverage/schema 10; full score matrix and hidden permutation 40 under the frozen preprocessing/tolerance; bijection/direction 15; report consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.
- Readiness closure: feature mapping, top-2,000 RNA variance selection, cosine score matrix and Hungarian one-to-one permutation are frozen in the task-local rule and oracle.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`multiomics_integration`, `regulatory_region_analysis`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

