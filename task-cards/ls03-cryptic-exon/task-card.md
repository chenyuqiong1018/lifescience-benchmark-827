# 任务卡： `ls03-cryptic-exon`

> Canonical participant-facing standalone card. The packaged-input inventory below matches the frozen task directory. Only the Prompt is pasted into a run; evaluator-only answers and outputs are never exposed.

## LS03-1｜高表达隐蔽外显子识别 — `ls03-cryptic-exon`

**Formal status:** `ready` — GRCh38 chr9, Ensembl 112 coding exons, junction supports and accepted static oracle are frozen.

| 字段 | 内容 |
| --- | --- |
| ID | ls03-cryptic-exon |
| Domain / sub-domain | transcriptomics / splicing |
| Level / time | L3, 90 min |
| Priority | P0 |
| Anchor / related capabilities | D / P, T, A, R, O |
| 来源思想 | CompBioBench cryptic-exon-q1 |
| Card version | standalone-v2 |

### Inputs (authoritative packaged inventory)
- `inputs/cryptic.exon.q1.fq.gz` — 16,920,968 bytes
- `inputs/reference/GRCh38_chr9.fa.gz` — 36,554,007 bytes
- `inputs/reference/README.md` — 505 bytes
- `inputs/reference/ensembl112_protein_coding_exons.tsv.gz` — 6,339,899 bytes

**Total:** 59,815,379 bytes (57.04 MiB).

### Prompt（运行时仅复制本框）

> Identify the protein-coding HGNC gene containing the highly expressed cryptic exon supported by two novel splice junctions. Write output/cryptic_exon.tsv with gene,chrom,start,end,left_junction_reads,right_junction_reads,expression_evidence, output/junctions.tsv, output/analysis.py, and output/report.md. Novelty must be assessed against the supplied annotation version.

### Deliverables / Output contract

output/cryptic_exon.tsv

output/junctions.tsv

output/analysis.py

output/report.md.

所有 CSV/TSV/JSON 必须可解析；ID 唯一；数值 finite；缺失值使用空字段或 null。

可重跑脚本必须只使用相对 inputs/ 与 output/ 路径；grader 不直接 import 未受信提交代码。

### Hard gates

□ one gene/interval

□ two flanking novel junctions

□ read support finite

□ rerun


### DeterministicArtifactScore（0–80, authoritative v2）

| Component | Points | Deterministic requirement |
| --- | ---: | --- |
| Coverage / schema | 10 | Gene, exon and both junction artifacts are complete and parseable. |
| Core science | 40 | HGNC gene, exon coordinates and two novel junctions/read supports match frozen gold/tolerances. |
| Direction / decision | 15 | Cryptic-exon decision agrees with novelty and bilateral junction support. |
| Summary consistency | 5 | Report repeats gene/exon and evidence consistently. |
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

Skill: C0=`NONE`; T0=`AUTO`; T1=`transcriptome_analysis`, `ensembl-sequence-retrieval`, `genome_annotation`.

MCP=`NONE`. No additional MCPs are configured; expected to reduce errors in numeric values and shapes, convergence assessment, logging, and environment reproducibility.

