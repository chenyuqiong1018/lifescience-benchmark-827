# 任务卡： `ls04-differential-composition`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS04-1｜视网膜单细胞差异组成分析 — `ls04-differential-composition`

**Formal status:** `ready` — Marker panel, annotation rule, composition counts and depleted-cell gold are frozen.

| 字段 | 内容 |
| --- | --- |
| ID | ls04-differential-composition |
| Domain / sub-domain | single-cell / composition |
| Level / time | L3, 90 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, T, A, R, O |
| 来源思想 | CompBioBench differential-composition-q1 |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/ANNOTATION_RULE.md` — 724 bytes
- `inputs/MARKER_PANEL.tsv` — 673 bytes
- `inputs/differential.composition.q1.1.mtx.gz` — 29,076,744 bytes
- `inputs/differential.composition.q1.2.mtx.gz` — 31,011,024 bytes
- `inputs/differential.composition.q1.genes.txt.gz` — 223,316 bytes

**Total:** 60,312,481 bytes (57.52 MiB).

### Prompt（运行时仅复制本框）

> Compare the two retinal single-cell expression matrices and identify the cell population that is severely depleted in sample 2. Write output/composition.csv with sample,cell_type,n_cells,fraction, output/depleted_call.json, output/analysis.py, and output/report.md. Document QC, normalization, annotation evidence and uncertainty.

### Deliverables / Output contract

output/composition.csv

output/depleted_call.json

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ both samples

□ fractions valid/sum within tolerance

□ one call or explicit ambiguity

□ rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | Both samples and all frozen cell-type labels are represented; fractions parse. |
| Core science | 40 | Cell counts/fractions and depleted population match frozen annotation/composition gold. |
| Direction / decision | 15 | Depletion direction is sample-2 relative to sample-1 and agrees with fractions. |
| Summary consistency | 5 | Report states depleted population and magnitude consistently. |
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

