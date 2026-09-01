# report.md

## T1 条件
- 加载顺序：bioassay_analysis、unit_conversion_suite、measurement-error-analysis、protocol-generation-from-description、code_execution_analysis

## 孔板状态
- A1、A2、A3、B1：completed
- B2：取尖失败，且发生在吸液前；未移动液体
- 恢复范围：仅恢复 B2/B3；不重做已完成孔

## 恢复方案
- 目标：每孔 0.5 uM / 100 uL
- Source：A2，25 uM
- 操作：从 A2 取 2 uL，使用 P20；加入 98 uL media diluent R1，使用 P300
- 公式：C1V1 = C2V2  
  25 uM × 2 uL / 100 uL = 0.5 uM
- 库存：A2 source 与 R1 diluent 满足需求