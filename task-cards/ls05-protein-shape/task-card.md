# 任务卡： `ls05-protein-shape`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS05-1｜PDB 蛋白几何字母形状识别 — `ls05-protein-shape`

**Formal status:** `ready_local_extension` — Geometry/PCA shape gold and view-metadata contract are frozen; report separately as a local-extension calibration task.

| 字段 | 内容 |
| --- | --- |
| ID | ls05-protein-shape |
| Domain / sub-domain | structural biology / geometry health check |
| Level / time | L1, 20 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, V, O |
| 来源思想 | CompBioBench protein-shape-q1 |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/protein.shape.q1.pdb` — 1,558,926 bytes

**Total:** 1,558,926 bytes (1.49 MiB).

### Prompt（运行时仅复制本框）

> Inspect the supplied PDB geometry and determine which one of B,D,F,H,J,L,N,P,R,T,V,X,Z it most resembles. Write output/shape_call.json with letter,confidence,orientation_notes and output/shape_view.png. Use only the supplied structure.

### Deliverables / Output contract

output/shape_call.json

output/shape_view.png.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

### Hard gates

□ allowed letter

□ valid nonempty PNG

□ confidence 0–1


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | Call JSON and PNG exist and parse; allowed letter vocabulary enforced. |
| Core science | 40 | Letter identity matches frozen visual/geometric gold. |
| Direction / decision | 15 | Orientation/confidence decision is valid and consistent with call. |
| Summary consistency | 5 | Orientation note concisely supports the same letter. |
| Script / reproducibility | 10 | No script points for L1 unless optional script meets standard; remaining points require valid reproducible view metadata. |

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

Skill: C0=`NONE`; T0=`AUTO`; T1=`protein_structure_analysis`, `protein_quality_assessment`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

