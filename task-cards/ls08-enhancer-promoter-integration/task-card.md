# 任务卡： `ls08-enhancer-promoter-integration`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## `ls08-enhancer-promoter-integration`

| 字段 | 内容 |
| --- | --- |
| ID | ls08-enhancer-promoter-integration |
| Domain / sub-domain | regulatory genomics / multimodal evidence integration |
| Level / time | L2, 45 min |
| Anchor / related capabilities | [regulatory-integration] |
| 来源思想 | CompBioBench / ep-interactions-q1 |

### Inputs (authoritative packaged inventory)
- `inputs/INTEGRATION_RULE.md` — 1,011 bytes
- `inputs/README.md` — 1,131 bytes
- `inputs/ep.interactions.q1.expr.csv` — 3,792 bytes
- `inputs/ep.interactions.q1.hic.csv` — 6,388 bytes

**Total:** 12,322 bytes (0.01 MiB).

### Prompt

> Integrate the supplied Hi-C and CRISPR-expression evidence for all candidate enhancer-promoter pairs and identify the least supported causal pair. Write `output/pair_evidence.csv` with `pair_id,contact_evidence,perturbation_effect,combined_support,rank`, `output/least_supported.json`, `output/analysis.py`, and `output/report.md`. Treat physical contact and perturbation evidence as distinct.

### Deliverables

every candidate pair exactly once, unique ranks, one least-supported call, report and rerunnable script; units/scales and tie policy stated.

### Hard gates

modalities joined by the true pair key; all candidates covered once; least-supported call equals minimum combined support under the frozen rule; contact is not described as perturbational proof.

### DeterministicArtifactScore（0–80）

coverage/schema 10; modality values, combined score and ranking 40; least-supported decision 15; report consistency 5; static/rerunnable script 10.

### JudgeScore（0–20）

| Dimension | 0 | 3 | 5 |
| --- | --- | --- | --- |
| Evidence | No task-specific evidence, fabricated/untraceable support, or contradiction with machine-readable artifacts. | At least one traceable task-specific input/result supports the main claim, but coverage, linkage or uncertainty is incomplete. | Every main claim links to the relevant input/result artifact; decisive measurements/rows and evidence-bound uncertainty are explicit. |
| Method | Missing or materially invalid method. | Broadly valid method with incomplete choices, direction, units or limitations. | Correct, auditable method with all consequential choices, direction, units and limitations stated. |
| Restraint | Unsupported causal/clinical/experimental claim. | Mostly calibrated but one important boundary is vague. | Claims stay within the supplied evidence and all important limitations are explicit. |
| Readability | Unusable or internally confusing report. | Understandable with notable ambiguity or clutter. | Concise, internally consistent and easy to audit against formal artifacts. |

The judge records main-claim count, supported-claim count, evidence locations and a short rationale. Judge points cannot repair a deterministic hard-gate failure.
- Readiness closure: background distance regression, robust residual z-score, guide aggregation, physical-evidence threshold and ranking/tie rule are frozen and independently recomputed by the checker.

### Ablation（不进入 Prompt）

Skill: C0=`NONE`; T0=`AUTO`; T1=`regulatory_region_analysis`, `region-gene-elements`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

