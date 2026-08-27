# 任务卡： `ls02-infer-genome-build`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS02-3｜VCF 基因组版本推断 — `ls02-infer-genome-build`

**Formal status:** `ready` — hg18/hg19/hg38 chr20 references and match-rate gold are frozen.

| 字段 | 内容 |
| --- | --- |
| ID | ls02-infer-genome-build |
| Domain / sub-domain | genomics / coordinate normalization |
| Level / time | L2, 40 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, T, A, O |
| 来源思想 | CompBioBench vcf-infer-build-q1 |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/references/README.md` — 563 bytes
- `inputs/references/hg18_chr20.fa.gz` — 19,496,195 bytes
- `inputs/references/hg19_chr20.fa.gz` — 18,072,551 bytes
- `inputs/references/hg38_chr20.fa.gz` — 18,840,364 bytes
- `inputs/references/reference_manifest.json` — 1,108 bytes
- `inputs/vcf.infer.build.q1.vcf.gz` — 3,788,749 bytes

**Total:** 60,199,530 bytes (57.41 MiB).

### Prompt（运行时仅复制本框）

> Determine whether the supplied chr20 VCF uses hg18, hg19, hg38 or T2T coordinates. Write output/build_call.json with build,confidence,n_variants_checked,n_ref_matches,n_ref_mismatches,evidence, output/analysis.py, and output/report.md. Base the call on reproducible allele/coordinate checks and do not treat chromosome naming alone as proof.

### Deliverables / Output contract

output/build_call.json

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ allowed build label

□ counts nonnegative

□ evidence present

□ rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | Build artifact parses and diagnostic counts are complete. |
| Core science | 40 | Reference-allele checks and final build match frozen diagnostic loci/gold. |
| Direction / decision | 15 | Confidence/call direction agrees with match-versus-mismatch evidence. |
| Summary consistency | 5 | Report states build and strongest evidence consistently. |
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

Skill: C0=`NONE`; T0=`AUTO`; T1=`chromosome_analysis`, `ucsc_genome_exploration`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

