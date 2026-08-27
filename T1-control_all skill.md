# T1-control_all Skill 指定方案

## 适用范围

- 候选池严格采用飞书 5.2「生命科学」模板中的 183 个 Skill：143 个生命科学 Skill + 40 个通用 Skill。
- 每题固定指定 1–5 个与核心方法、输入格式或交付物直接相关的 Skill；允许复用，但不为凑数加入无关 Skill。
- 同一题的所有 T1 重复运行必须使用完全相同的 Skill 集合和表内顺序。
- 本方案只改变 T1 的指定 Skill；不配置额外 MCP，task-card、inputs、模型和其他运行参数保持不变。

## 指定方案

| 序号 | Task ID | T1 指定 Skill |
|---:|---|---|
| 1 | `ls01-grna-offtarget-rank` | `dna-rna-sequence-analysis`; `genome_annotation`; `code_execution_analysis` |
| 2 | `ls01-primer-transcript-audit` | `dna-rna-sequence-analysis`; `ensembl-sequence-retrieval`; `transcriptome_analysis`; `code_execution_analysis` |
| 3 | `ls01-vector-orf-audit` | `dna-rna-sequence-analysis`; `synthetic_biology_design`; `code_execution_analysis` |
| 4 | `ls02-deleterious-mutation` | `genome_annotation`; `snp_functional_analysis`; `statistical-analysis`; `code_execution_analysis` |
| 5 | `ls02-find-deletion` | `chromosome_analysis`; `ucsc_genome_exploration`; `genome_annotation`; `code_execution_analysis` |
| 6 | `ls02-infer-genome-build` | `chromosome_analysis`; `ucsc_genome_exploration`; `code_execution_analysis` |
| 7 | `ls03-atac-sample-swap` | `regulatory_region_analysis`; `tissue_specific_analysis`; `exploratory-data-analysis`; `statistical-analysis`; `scientific-visualization` |
| 8 | `ls03-cryptic-exon` | `transcriptome_analysis`; `ensembl-sequence-retrieval`; `genome_annotation`; `code_execution_analysis` |
| 9 | `ls03-genome-coordinates` | `regulatory_region_analysis`; `region-gene-elements`; `statistical-analysis`; `scientific-visualization` |
| 10 | `ls04-differential-composition` | `scvi-tools`; `scgpt`; `statistical-analysis`; `scientific-visualization` |
| 11 | `ls04-perturbseq-reference-map` | `scvi-tools`; `scgpt`; `statistical-analysis`; `code_execution_analysis` |
| 12 | `ls04-spatial-deconvolution` | `scvi-tools`; `scgpt`; `statistical-analysis`; `code_execution_analysis`; `scientific-visualization` |
| 13 | `ls05-protein-shape` | `protein_structure_analysis`; `protein_quality_assessment`; `code_execution_analysis`; `matplotlib` |
| 14 | `ls05-structure-model-ranking` | `protein_quality_assessment`; `protein_structure_analysis`; `statistical-analysis`; `code_execution_analysis` |
| 15 | `ls05-low-confidence-pocket` | `protein_quality_assessment`; `binding_site_characterization`; `alphafold_structure_pipeline`; `scientific-critical-thinking`; `code_execution_analysis` |
| 16 | `ls06-eno1-effect-size` | `proteome_analysis`; `markitdown`; `exploratory-data-analysis`; `statistical-analysis`; `code_execution_analysis` |
| 17 | `ls06-eno1-significance-audit` | `proteome_analysis`; `biomarker_discovery`; `markitdown`; `statistical-analysis`; `code_execution_analysis` |
| 18 | `ls07-combination-treatment-deg` | `biomarker_discovery`; `transcriptome_analysis`; `statistical-analysis`; `code_execution_analysis` |
| 19 | `ls07-combination-treatment-mechanism` | `go_term_analysis`; `string-ppi-enrichment`; `scientific-critical-thinking`; `code_execution_analysis` |
| 20 | `ls08-enhancer-promoter-integration` | `regulatory_region_analysis`; `region-gene-elements`; `scientific-critical-thinking`; `code_execution_analysis` |
| 21 | `ls08-multiome-column-match` | `multiomics_integration`; `regulatory_region_analysis`; `exploratory-data-analysis`; `statistical-analysis`; `code_execution_analysis` |
| 22 | `ls09-opentrons-sop` | `synthetic_biology_design`; `protocol-generation-from-description`; `protocol-to-executable-json`; `unit_conversion_suite`; `code_execution_analysis` |
| 23 | `ls09-plate-dilution-recovery` | `bioassay_analysis`; `unit_conversion_suite`; `measurement-error-analysis`; `protocol-generation-from-description`; `code_execution_analysis` |
| 24 | `ls10-neun-power-analysis` | `biomarker_discovery`; `statistical-analysis`; `code_execution_analysis` |
| 25 | `ls10-treatment-response-model` | `personalized_medicine`; `exploratory-data-analysis`; `statistical-analysis`; `code_execution_analysis` |

## 运行约束

1. Skill 只提供方法与工具指导；若 Skill 默认流程与 task-card 冲突，以 task-card、冻结 inputs 和 oracle 规则为准。
2. 禁止用在线数据库或 Skill 示例数据替换随题材料；涉及固定版本、阈值、排序公式、字段名或输出结构时，严格执行题面合同。
3. 通用 Skill 的加入用于补足表格解析、统计检验、代码计算、误差审计、协议结构化和图形输出能力，不改变生命科学领域 Skill 的证据边界。
4. 记录实际加载的 Skill 名称及顺序，确保 T1 可复现，并可与 C0=`NONE`、T0=`AUTO` 公平比较。

## 核验摘要

- 题目数：25
- 每题 Skill 数：3–5
- 候选池：183（生命科学 143 + 通用 40）
- 表内所用生命科学 Skill：均来自序号 1–143
- 表内所用通用 Skill：均来自序号 144–183
- 额外 MCP：无
