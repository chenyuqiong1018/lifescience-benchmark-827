# T0 多组学匹配报告

- **条件**: T0  
- **自动加载 Skill**: `multiomics_integration`, `regulatory_region_analysis`, `exploratory-data-analysis`, `statistical-analysis`, `code_execution_analysis`  
- **候选 Skill 总数**: 183 个，自动选择其中 5 个。

## 方法摘要

1. **唯一基因**：仅保留唯一基因标识，避免重复基因或多映射基因干扰后续相关性计算。  
2. **strand-aware TSS**：以基因链方向感知的 TSS 为锚点，使调控区信号与转录方向一致。  
3. **10kb bin**：在 TSS 附近按 10 kb 窗口聚合 ATAC / 调控区信号。  
4. **log1p 变换**：对信号强度做 `log1p` 变换，以稳定低值区间方差。  
5. **2000 高方差基因**：选择 2000 个高方差基因，保留具有生物学区分度的变异信号。  
6. **Pearson 相关**：计算每个 RNA population 与每个 ATAC column 之间的 Pearson 相关，作为匹配得分。  
7. **Hungarian 双射匹配**：使用 Hungarian 算法进行一一对应匹配，确保每个 RNA population 唯一匹配到一个 ATAC column，并最大化总匹配得分。

## 匹配结果

| rna_population | atac_column | match_score | runner_up_score | margin |
|---:|---:|---:|---:|---:|
| 0 | 5 | 0.547 | 0.535 | 0.012 |
| 1 | 1 | 0.509 | 0.488 | 0.021 |
| 2 | 4 | 0.526 | 0.461 | 0.065 |
| 3 | 0 | 0.453 | 0.390 | 0.063 |
| 4 | 6 | 0.383 | 0.343 | 0.040 |
| 5 | 3 | 0.414 | 0.335 | 0.079 |
| 6 | 7 | 0.390 | 0.277 | 0.113 |
| 7 | 2 | 0.400 | 0.308 | 0.091 |

- 所有 8 个 RNA population 均形成唯一匹配，无重复 ATAC column。  
- `match_score` 范围为 **0.383–0.547**，平均约 **0.453**。  
- 每个最优匹配均高于对应 `runner_up_score`，平均领先约 **0.061**。  
- 最强匹配为 **rna_population 0 → atac_column 5**，得分 **0.547**。  
- 领先幅度最大的匹配为 **rna_population 6 → atac_column 7**，领先 **0.113**。

## 生物信号解释

在 T0 条件下，RNA population 与 ATAC column 在 strand-aware TSS 附近 10 kb 调控窗口内呈现稳定正相关。经过唯一基因过滤、log1p 变换和 2000 高方差基因筛选后，Pearson 相关仍能支持 Hungarian 算法得到清晰的一一对应关系，说明 RNA 状态与局部染色质可及性之间存在可检测的配对特异性调控信号，而非随机匹配结果。