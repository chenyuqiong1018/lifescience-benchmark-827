# report.md

## T1 条件

指定顺序加载：personalized_medicine、exploratory-data-analysis、statistical-analysis、code_execution_analysis。

## 数据与样本

- 分析类型：完整病例分析  
- 纳入变量：Efficacy、Age、Gender、BMI  
- 结局编码：PR = 1，SD/PD = 0  
- 样本量：n = 80  
- 性别分布：41 / 39  
- Gender 参考水平：Female  

## 统计模型

采用 logistic 回归模型，结局为 Efficacy，模型变量包含：

Efficacy ~ BMI + Age + Gender

## Age 结果

- beta = -0.0795085  
- SE = 0.0262977  
- z = -3.0234  
- p = 0.00249954  
- OR = 0.92357 / year  

## 解释

在控制 BMI 和 Gender 后，Age 与疗效结局存在显著负相关。年龄每增加 1 岁，达到 PR 相对于 SD/PD 的优势乘以 0.92357，约下降 7.6%。该结果表示在调整 BMI 和性别后的独立关联。