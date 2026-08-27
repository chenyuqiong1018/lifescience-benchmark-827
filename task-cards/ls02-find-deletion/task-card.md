# 任务卡： `ls02-find-deletion`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS02-2｜浅层双端测序大片段缺失定位 — `ls02-find-deletion`

**Formal status:** `ready` — GRCh38 chr22, breakpoint evidence, tolerances and accepted static oracle are frozen.

| 字段 | 内容 |
| --- | --- |
| ID | ls02-find-deletion |
| Domain / sub-domain | genomics / structural variation |
| Level / time | L3, 75 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, T, A, R, O |
| 来源思想 | CompBioBench find-deletion-q1 |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/find.deletion.r1.fq.gz` — 16,705,359 bytes
- `inputs/find.deletion.r2.fq.gz` — 16,701,241 bytes
- `inputs/reference/GRCh38_chr22.fa.gz` — 11,147,284 bytes
- `inputs/reference/README.md` — 295 bytes

**Total:** 44,554,179 bytes (42.49 MiB).

### Prompt（运行时仅复制本框）

> Locate the large deletion in the supplied shallow paired-end hg38 data. Write output/deletion.tsv with chrom,start_100kb,end_100kb,size_bp,supporting_signals, output/qc.json, output/analysis.py, and output/report.md. Round breakpoints to the nearest 100 kb and distinguish evidence from precision limits.

### Deliverables / Output contract

output/deletion.tsv

output/qc.json

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ one interval, start<end, hg38 coordinates, support nonempty, rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | Required interval/QC artifacts parse and contain one nominated event. |
| Core science | 40 | Chromosome, rounded breakpoints and size match frozen gold within declared tolerances; support fields are present. |
| Direction / decision | 15 | Deletion/no-deletion decision agrees with coverage/pair evidence. |
| Summary consistency | 5 | Report repeats interval and precision limit consistently. |
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

Skill: C0=`NONE`; T0=`AUTO`; T1=`chromosome_analysis`, `ucsc_genome_exploration`, `genome_annotation`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

