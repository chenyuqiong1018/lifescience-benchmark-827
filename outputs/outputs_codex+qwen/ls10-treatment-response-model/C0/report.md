# ls10-treatment-response-model 报告

## 数据与编码

- 分析样本：完整病例，仅按 Efficacy、Age、Gender、BMI 纳入。
- 样本量：n = 80。
- 应答分布：response 41 例，nonresponse 39 例。
- Efficacy 编码：PR = 1，表示应答；SD/PD = 0，表示未应答。
- Gender 编码：Female 为参考类别。

## 模型设定

- 模型类型：Logistic 回归。
- 因变量：Efficacy，PR 相对 SD/PD。
- 预测变量：BMI、Age、Gender。
- 当前报告重点为 Age 的系数及其解释；BMI 和 Gender 作为协变量保留在模型中。

## Age 结果

| 指标 | 数值 |
| --- | ---: |
| Age coefficient | -0.0795084690 |
| SE | 0.0262977303 |
| z | -3.0233966 |
| two-sided p | 0.0024995441 |
| OR per year | 0.9235702 |

## 解释

在控制 BMI 和 Gender 后，Age 每增加 1 岁，获得 PR 相对 SD/PD 的优势乘以 0.9235702，即应答优势约降低 7.6%。

负系数表示：在相同 BMI 和相同 Gender 条件下，年龄越大，达到 PR 的可能性越低。two-sided p = 0.0024995441，说明 Age 与应答状态之间的关联具有统计学显著性。

该结果反映的是调整 BMI 和 Gender 后的年龄独立关联，不直接等同于因果效应。