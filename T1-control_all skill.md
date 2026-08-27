# T1 Control — All Skills

## 实验条件定义

- 条件名称：`T1-control_all-skill`
- 任务范围：物质科学评测集 30 题（`MS01-Q1`—`MS10-Q3`）。
- 固定候选池：物质科学模板中的全部 96 个 Skill，其中物质专用 48 个、通用 48 个。
- 允许范围：每道题只能从下方 96 个 Skill 中选择，不得调用候选池外能力。
- 选择方式：运行前依据题目 prompt、inputs 类型和交付物，从固定候选池选择最相关的 1—4 个 Skill；禁止把 96 个 Skill 全部注入上下文。
- 隔离要求：每次运行前清空已安装或已启用的 Skill，只安装并打开本题实际选中的 Skill；运行结束后再次清空。
- 审计要求：记录候选池版本、实际选择结果、选择理由、安装/启用/清理时间及失败回退。
- 无匹配回退：若没有合适 Skill，记录 `no_matching_skill` 并以无 Skill 方式完成，不得使用候选池外能力。

## 30 题适用范围

| 序号 | Task ID | T1 指定候选池 | 单题选择上限 |
| ---: | --- | --- | ---: |
| 1 | `MS01-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 2 | `MS01-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 3 | `MS01-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 4 | `MS02-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 5 | `MS02-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 6 | `MS02-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 7 | `MS03-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 8 | `MS03-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 9 | `MS03-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 10 | `MS04-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 11 | `MS04-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 12 | `MS04-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 13 | `MS05-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 14 | `MS05-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 15 | `MS05-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 16 | `MS06-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 17 | `MS06-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 18 | `MS06-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 19 | `MS07-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 20 | `MS07-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 21 | `MS07-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 22 | `MS08-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 23 | `MS08-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 24 | `MS08-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 25 | `MS09-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 26 | `MS09-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 27 | `MS09-Q3` | 全部 96 个物质＋通用 Skill | 4 |
| 28 | `MS10-Q1` | 全部 96 个物质＋通用 Skill | 4 |
| 29 | `MS10-Q2` | 全部 96 个物质＋通用 Skill | 4 |
| 30 | `MS10-Q3` | 全部 96 个物质＋通用 Skill | 4 |

## 96 Skill 候选池

| 序号 | Skill 名称 | 类型 |
| ---: | --- | --- |
| 1 | `admet_druglikeness_report` | 物质专用 |
| 2 | `aliphatic_ring_analysis` | 物质专用 |
| 3 | `buoyancy-acceleration-calculation` | 物质专用 |
| 4 | `capacitance-calculation` | 物质专用 |
| 5 | `cas_compound_lookup` | 物质专用 |
| 6 | `chembl-molecule-search` | 物质专用 |
| 7 | `chemical-mass-percent-calculation` | 物质专用 |
| 8 | `chemical-structure-analysis` | 物质专用 |
| 9 | `chemical_patent_analysis` | 物质专用 |
| 10 | `chemical_property_profiling` | 物质专用 |
| 11 | `chemical_safety_assessment` | 物质专用 |
| 12 | `chemical_structure_comparison` | 物质专用 |
| 13 | `combinatorial_chemistry` | 物质专用 |
| 14 | `compound-name-retrieval` | 物质专用 |
| 15 | `compound_database_crossref` | 物质专用 |
| 16 | `drugsda-admet` | 物质专用 |
| 17 | `drugsda-compound-retrieve` | 物质专用 |
| 18 | `drugsda-denovo-sampling` | 物质专用 |
| 19 | `drugsda-drug-likeness` | 物质专用 |
| 20 | `drugsda-linker-sampling` | 物质专用 |
| 21 | `drugsda-mol-properties` | 物质专用 |
| 22 | `drugsda-mol-similarity` | 物质专用 |
| 23 | `drugsda-mol2mol-sampling` | 物质专用 |
| 24 | `drugsda-rgroup-sampling` | 物质专用 |
| 25 | `electrical_circuit_analysis` | 物质专用 |
| 26 | `electromagnetic_analysis` | 物质专用 |
| 27 | `energy_conversion` | 通用 |
| 28 | `experimental_data_processing` | 通用 |
| 29 | `functional_group_profiling` | 物质专用 |
| 30 | `geometric-volume-calculation` | 通用 |
| 31 | `geometry_trigonometry` | 通用 |
| 32 | `lead_compound_optimization` | 物质专用 |
| 33 | `length_measurement` | 通用 |
| 34 | `material-density-volume-calculation` | 物质专用 |
| 35 | `mobility_analysis` | 物质专用 |
| 36 | `molecular-descriptors-calculation` | 物质专用 |
| 37 | `molecular-format-conversion` | 物质专用 |
| 38 | `molecular-properties-calculation` | 物质专用 |
| 39 | `molecular-property-profiling` | 物质专用 |
| 40 | `molecular-similarity-search` | 物质专用 |
| 41 | `molecular_fingerprint_analysis` | 物质专用 |
| 42 | `molecular_visualization_suite` | 物质专用 |
| 43 | `natural_product_analysis` | 物质专用 |
| 44 | `nuclear_physics` | 物质专用 |
| 45 | `optical-frequency-calculation` | 物质专用 |
| 46 | `optics_analysis` | 物质专用 |
| 47 | `polymer_property_analysis` | 物质专用 |
| 48 | `pubchem-smiles-search` | 物质专用 |
| 49 | `pubchem_deep_dive` | 物质专用 |
| 50 | `signal_processing` | 通用 |
| 51 | `smiles-to-cas-conversion` | 物质专用 |
| 52 | `smiles_comprehensive_analysis` | 物质专用 |
| 53 | `statistical_error_analysis` | 通用 |
| 54 | `substructure_activity_search` | 物质专用 |
| 55 | `thermal_analysis` | 物质专用 |
| 56 | `unit-conversion-nanoscale` | 通用 |
| 57 | `exploratory-data-analysis` | 通用 |
| 58 | `statistical-analysis` | 通用 |
| 59 | `markitdown` | 通用 |
| 60 | `openalex-database` | 通用 |
| 61 | `hypothesis-generation` | 通用 |
| 62 | `initialize-atlas-graph` | 通用 |
| 63 | `peer-review` | 通用 |
| 64 | `scientific-critical-thinking` | 通用 |
| 65 | `matplotlib` | 通用 |
| 66 | `scientific-visualization` | 通用 |
| 67 | `seaborn` | 通用 |
| 68 | `citation-management` | 通用 |
| 69 | `compute-env-setup` | 通用 |
| 70 | `customize` | 通用 |
| 71 | `figure-composer` | 通用 |
| 72 | `figure-style` | 通用 |
| 73 | `literature-review` | 通用 |
| 74 | `managed-model-endpoints` | 通用 |
| 75 | `paper-narrative` | 通用 |
| 76 | `pdf-explore` | 通用 |
| 77 | `product-self-knowledge` | 通用 |
| 78 | `remote-compute-modal` | 通用 |
| 79 | `remote-compute-ssh` | 通用 |
| 80 | `self-awareness` | 通用 |
| 81 | `skill-creator` | 通用 |
| 82 | `using-model-endpoint` | 通用 |
| 83 | `code_execution_analysis` | 通用 |
| 84 | `lab_protocol_from_literature` | 通用 |
| 85 | `measurement-error-analysis` | 通用 |
| 86 | `meta-analysis-execution` | 通用 |
| 87 | `protocol-extraction-from-pdf` | 通用 |
| 88 | `protocol-generation-from-description` | 通用 |
| 89 | `protocol-to-executable-json` | 通用 |
| 90 | `scientific-literature-search` | 通用 |
| 91 | `unit_conversion_suite` | 通用 |
| 92 | `web_literature_mining` | 通用 |
| 93 | `nature-figure` | 通用 |
| 94 | `drugsda-data-valid` | 通用 |
| 95 | `drugsda-file-transfer` | 通用 |
| 96 | `pubmed-article-search` | 通用 |

## 单题执行记录格式

| Task ID | 实际选择 Skill | 选择理由 | 安装成功 | 已启用 | 运行后已清理 |
| --- | --- | --- | --- | --- | --- |
| `MSxx-Qx` | `skill_a`; `skill_b` | 与 prompt、inputs 和交付物的对应关系 | 是/否 | 是/否 | 是/否 |

## 版本来源

- Skill 清单来源：飞书《项目模板》“物质科学模板（通用＋物质）”章节。
- 本文件固定 96 项全集作为 T1 可选范围；后续若清单变更，应更新文件并记录版本差异。

