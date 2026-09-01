```markdown
# report.md

## T1 条件与加载顺序
按指定顺序加载：  
`multiomics_integration` → `regulatory_region_analysis` → `exploratory-data-analysis` → `statistical-analysis` → `code_execution_analysis`

## 计算设置
- **唯一基因**：基因标识去重后保留唯一基因符号，避免重复特征干扰后续统计与匹配。  
- **strand-aware TSS**：以链感知的转录起始位点为锚点，按基因链方向定义调控区域。  
- **10kb bin**：将 TSS 附近调控信号按 10 kb 分箱汇总。  
- **log1p**：对计数/信号使用 `log1p(x) = log(1 + x)` 变换，以稳定方差。  
- **2000高方差基因**：选取 2000 个高方差基因构建跨组学特征矩阵。  
- **Pearson**：使用 Pearson 相关系数衡量 RNA population 与 ATAC column 的相似度。  
- **Hungarian双射**：使用 Hungarian 算法进行一对一最优匹配，保证每个 RNA population 唯一对应一个 ATAC column。

## 真实匹配结果
| rna_population | atac_column | match_score | runner_up_score | margin |
|---:|---:|---:|---:|---:|
| 0 | 5 | 0.546761 | 0.535126 | 0.011635 |
| 1 | 1 | 0.509316 | 0.488338 | 0.020978 |
| 2 | 4 | 0.526474 | 0.461484 | 0.064990 |
| 3 | 0 | 0.452778 | 0.389683 | 0.063095 |
| 4 | 6 | 0.382501 | 0.342562 | 0.039940 |
| 5 | 3 | 0.414109 | 0.334626 | 0.079483 |
| 6 | 7 | 0.389744 | 0.276556 | 0.113188 |
| 7 | 2 | 0.399550 | 0.308453 | 0.091097 |

## 结果摘要
- 8 个 `rna_population` 与 8 个 `atac_column` 形成 **完整双射**，`atac_column` 为 `0–7` 的唯一排列。  
- 所有匹配均为正相关，`match_score` 范围为 **0.3825–0.5468**。  
- 每个最优匹配均高于次优匹配，平均 margin 为 **0.0606**。  
- 最优总匹配分数为 **3.6212**。

## 结论
基于 strand-aware TSS、10 kb bin、log1p 变换、2000 高方差基因和 Pearson 相关，Hungarian 算法得到唯一且稳定的双射匹配。结果支持 RNA population 与对应 ATAC column 之间存在共享的调控/生物信号。
```