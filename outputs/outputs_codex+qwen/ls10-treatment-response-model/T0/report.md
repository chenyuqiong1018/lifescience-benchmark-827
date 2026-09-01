# report.md

## 分析条件
- T0 条件下自动选择：personalized_medicine、exploratory-data-analysis、statistical-analysis、code_execution_analysis。
- 分析对象：完整病例，包含变量 Efficacy、Age、Gender、BMI。
- 结局编码：PR = 1，SD/PD = 0。
- Gender 以 Female 为参考组。
- 样本量：n = 80；报告计数为 41/39。

## 统计模型
- 采用 logistic 回归，模型包含 BMI、Age、Gender。
- 模型形式：Efficacy ~ BMI + Age + Gender。

## 主要结果
| 变量 | beta | SE | z | p | OR |
|---|---:|---:|---:|---:|---:|
| Age | -0.0795085 | 0.0262977 | -3.0234 | 0.00249954 | 0.92357 / year |

## 解释
在控制 BMI 和 Gender 后，Age 与疗效结局呈显著负相关。

Age 每增加 1 岁，达到 PR 相对于 SD/PD 的 odds 乘以 0.92357，约下降 7.6%。这表示在 BMI 和 Gender 相同的情况下，年龄越大，达到 PR 的相对可能性越低。

该结果应解释为控制 BMI 和 Gender 后的统计关联，不直接等同于因果效应。