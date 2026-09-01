# ls08-multiome-column-match 最终报告  
## C0 条件，不使用任何 Skill

本报告基于 ls08-multiome-column-match 在 C0 条件下、未使用任何 Skill 的真实计算结果，总结 RNA population 与 ATAC column 之间的最终匹配关系，并说明所采用的计算策略、匹配得分含义、共享生物信号以及 runner-up 得分的解释。

---

## 1. 任务目标

该任务的目标是：在 C0 条件下，根据 RNA 与 ATAC 特征之间的相关性，将每个 RNA population 匹配到一个唯一的 ATAC column，并尽可能恢复两组多组学表示之间的对应关系。

匹配过程不依赖外部 Skill，仅使用输入数据本身的特征计算。最终采用 Hungarian 算法得到双射匹配，即每个 RNA population 对应一个 ATAC column，每个 ATAC column 也只被分配一次。

---

## 2. 计算方法

### 2.1 唯一基因

计算中采用唯一基因作为基因侧特征，即对基因进行去重和无歧义化处理，避免重复基因、同名基因或多映射基因对后续相关性计算造成干扰。每个基因在特征矩阵中只保留一个明确条目。

### 2.2 Strand-aware TSS

在构建基因相关的染色质可及性特征时，使用 strand-aware TSS，即根据基因所在链确定转录起始位点方向。

对于正链基因，TSS 通常取基因转录方向的起始位置；对于负链基因，TSS 方向需要按其转录方向处理。这样可以避免简单使用基因组坐标起点或终点造成的方向错误，使 ATAC 信号与基因调控区域的关系更加准确。

### 2.3 10kb ATAC bin

ATAC 信号按照 10kb bin 进行聚合。即以 10kb 长度的区间作为染色质可及性特征的基本汇总单位，将区间内的 ATAC 信号整理为可用于相关性计算的特征。

这种 binning 方式可以在保留局部染色质可及性信息的同时降低单点噪声，并使 RNA 与 ATAC 特征之间更容易进行稳定的相关性比较。

### 2.4 log1p 变换

对 RNA 和 ATAC 信号使用 log1p 变换，即：

log1p(x) = log(1 + x)

该变换用于压缩高信号值的尺度，降低极端计数或强峰值对相关性计算的主导影响，同时保留低信号区域的信息。log1p 对零值友好，不会引入 log(0) 的问题。

### 2.5 RNA 方差最高 2000 基因

RNA 侧特征并非使用全部基因，而是选择方差最高的 2000 个基因。

这些基因通常在不同 population 之间具有更强的表达差异，更能代表不同 RNA population 的生物学状态。使用高方差基因可以提高匹配信号的信噪比，减少低信息量基因对相关性估计的稀释。

### 2.6 Pearson 相关

对于每个 RNA population 与每个候选 ATAC column，计算二者之间的 Pearson 相关系数。

Pearson 相关用于衡量两个特征向量之间的线性共变关系。在本任务中，较高的 Pearson 相关表示该 RNA population 与该 ATAC column 在特征变化模式上更一致，提示二者可能来自同一生物学状态或同一组学对应关系。

### 2.7 Hungarian 双射匹配

在得到 RNA population 与 ATAC column 的相关性矩阵后，使用 Hungarian 算法进行全局最优匹配。

Hungarian 算法的目标是最大化总体匹配得分，同时满足双射约束：

- 每个 RNA population 只能匹配一个 ATAC column；
- 每个 ATAC column 只能被一个 RNA population 匹配；
- 最终结果是一组一一对应的匹配关系。

因此，最终匹配不仅考虑单个 RNA population 的局部最高相关，也考虑全局匹配的一致性和总得分。

---

## 3. 最终匹配结果

以下为 ls08-multiome-column-match 的真实最终匹配结果：

| RNA population | 匹配的 ATAC column | match_score | runner_up_score | match_score - runner_up_score |
|---:|---:|---:|---:|---:|
| 0 | 5 | 0.5467608355031154 | 0.5351257852754059 | 0.0116350502 |
| 1 | 1 | 0.5093155435867853 | 0.4883380242307774 | 0.0209775194 |
| 2 | 4 | 0.5264740567092872 | 0.4614841993560838 | 0.0649898574 |
| 3 | 0 | 0.4527781428024127 | 0.3896828715340393 | 0.0630952713 |
| 4 | 6 | 0.3825014972318415 | 0.3425616797593513 | 0.0399398175 |
| 5 | 3 | 0.4141085717464719 | 0.3346257692770320 | 0.0794828025 |
| 6 | 7 | 0.3897439385010290 | 0.2765558638985525 | 0.1131880746 |
| 7 | 2 | 0.3995499283923857 | 0.3084526204199870 | 0.0910973080 |

最终匹配关系为：

| RNA population | ATAC column |
|---:|---:|
| 0 | 5 |
| 1 | 1 |
| 2 | 4 |
| 3 | 0 |
| 4 | 6 |
| 5 | 3 |
| 6 | 7 |
| 7 | 2 |

可以看到，8 个 ATAC column 分别为 5、1、4、0、6、3、7、2，覆盖了 0 到 7 的全部 column，且每个 column 只出现一次，符合 Hungarian 双射匹配的要求。

---

## 4. 匹配得分总结

最终 match_score 的范围为：

- 最低：0.3825014972318415，对应 RNA population 4 匹配 ATAC column 6；
- 最高：0.5467608355031154，对应 RNA population 0 匹配 ATAC column 5。

所有最终 match_score 均为正相关，说明每个被匹配的 RNA population 与其对应 ATAC column 之间都存在一定程度的共享变化模式。

整体统计如下：

| 指标 | 数值 |
|---|---:|
| 平均 match_score | 0.4526540643 |
| 平均 runner_up_score | 0.3921033517 |
| 平均 match_score - runner_up_score | 0.0605507126 |

平均而言，最终匹配得分高于次优候选得分约 0.0606，说明 Hungarian 选择的最终匹配整体上优于备选匹配。

---

## 5. 共享生物信号的含义

match_score 使用 Pearson 相关来衡量 RNA population 与 ATAC column 之间的一致性。这里的正相关并不只是数学上的相似，而是反映了两种模态之间共享的生物学信号。

具体来说，如果某个 RNA population 与某个 ATAC column 匹配得分较高，说明二者在所选特征上具有相似的样本间变化趋势。也就是说，当某些基因表达在该 RNA population 中升高或降低时，对应 ATAC column 中的染色质可及性也倾向于呈现一致的变化。

这种共变关系通常提示：

- 二者可能来自相同的细胞群体或生物学状态；
- 基因表达变化与染色质开放程度变化之间存在调控层面的关联；
- RNA 所代表的转录状态与 ATAC 所代表的染色质状态捕获了相似的生物学轴；
- 该匹配结果具有生物学可解释性，而不仅是随机相似。

在本结果中，所有最终匹配得分均为正，且最高达到约 0.5468，说明 RNA 和 ATAC 之间确实存在可检测的共享生物信号。不同 population 的得分高低差异，反映了不同匹配对在 C0 条件下共享信号的强弱。

---

## 6. runner-up 得分的含义

runner_up_score 表示对应 RNA population 的次优候选匹配得分，可以理解为除最终 Hungarian 分配的 ATAC column 之外，最具竞争力的替代 ATAC column 所得到的相关分数。

runner_up_score 的意义在于衡量最终匹配的相对优势：

- 如果 match_score 明显高于 runner_up_score，说明最终匹配相对备选匹配更清晰，匹配置信度较高；
- 如果 match_score 与 runner_up_score 接近，说明最终匹配与某个替代匹配之间的差异较小，该 RNA population 的匹配相对更容易产生歧义；
- match_score - runner_up_score 可作为匹配稳定性的参考指标。

例如：

- RNA population 6 的 match_score 为 0.3897439385010290，runner_up_score 为 0.2765558638985525，差值约为 0.1132，是所有结果中最大的差值，说明 RNA population 6 匹配 ATAC column 7 的相对优势最明显。
- RNA population 0 的 match_score 为 0.5467608355031154，runner_up_score 为 0.5351257852754059，差值仅约 0.0116，说明虽然最终选择了 ATAC column 5，但其与次优候选之间的差距很小，匹配相对更接近歧义边界。
- RNA population 7 的差值约为 0.0911，RNA population 5 的差值约为 0.0795，说明这两个匹配也具有较强的相对优势。

因此，runner_up_score 并不是最终采用的匹配得分，而是用于评估最终匹配是否显著优于其他候选匹配的重要参考。

---

## 7. 各匹配结果解读

### RNA population 0 匹配 ATAC column 5

match_score 为 0.5467608355031154，是所有最终匹配中最高的，说明 RNA population 0 与 ATAC column 5 之间的相关性最强。

但其 runner_up_score 为 0.5351257852754059，与最终得分非常接近，差值仅约 0.0116。这说明 RNA population 0 虽然最终被分配到 ATAC column 5，但存在较强的竞争候选，匹配优势相对较小。

### RNA population 1 匹配 ATAC column 1

match_score 为 0.5093155435867853，runner_up_score 为 0.4883380242307774，差值约 0.0210。该匹配得分较高，但次优候选也相对接近，属于较强但存在一定竞争关系的匹配。

### RNA population 2 匹配 ATAC column 4

match_score 为 0.5264740567092872，runner_up_score 为 0.4614841993560838，差值约 0.0650。该匹配不仅得分较高，而且与次优候选之间有明显差距，属于较可靠的匹配。

### RNA population 3 匹配 ATAC column 0

match_score 为 0.4527781428024127，runner_up_score 为 0.3896828715340393，差值约 0.0631。该匹配具有中等强度的相关，并且相对次优候选有较明显优势。

### RNA population 4 匹配 ATAC column 6

match_score 为 0.3825014972318415，是所有最终匹配中最低的，但其 runner_up_score 为 0.3425616797593513，差值约 0.0399。这说明该匹配虽然绝对相关较弱，但最终选择仍优于次优候选。

### RNA population 5 匹配 ATAC column 3

match_score 为 0.4141085717464719，runner_up_score 为 0.3346257692770320，差值约 0.0795。该匹配具有较好的相对优势。

### RNA population 6 匹配 ATAC column 7

match_score 为 0.3897439385010290，runner_up_score 为 0.2765558638985525，差值约 0.1132。虽然绝对得分不是最高，但其相对次优候选的优势最大，说明该匹配在候选集合中较为明确。

### RNA population 7 匹配 ATAC column 2

match_score 为 0.3995499283923857，runner_up_score 为 0.3084526204199870，差值约 0.0911。该匹配也具有较明显的相对优势。

---

## 8. 结论

在 C0 条件且不使用任何 Skill 的情况下，ls08-multiome-column-match 通过唯一基因、strand-aware TSS、10kb ATAC bin、log1p 变换、RNA 方差最高 2000 基因、Pearson 相关以及 Hungarian 双射匹配，得到了 8 组 RNA population 与 ATAC column 的一一对应关系。

最终结果表明：

- 所有 RNA population 均成功匹配到唯一的 ATAC column；
- 所有最终 match_score 均为正，说明 RNA 与 ATAC 之间存在共享生物学信号；
- 平均 match_score 高于平均 runner_up_score，说明 Hungarian 全局匹配整体有效；
- RNA population 0 的最终匹配绝对得分最高，但与 runner-up 差距最小；
- RNA population 6 的最终匹配绝对得分不高，但相对 runner-up 优势最大；
- 整体结果支持 RNA population 与 ATAC column 之间存在可恢复的对应关系。

综上，该真实计算结果说明，在不使用任何 Skill 的 C0 条件下，基于高方差基因、TSS 相关 ATAC 特征、log1p 变换、Pearson 相关和 Hungarian 双射匹配，能够从 RNA 与 ATAC 特征中恢复出一组稳定且可解释的 column 对应关系。