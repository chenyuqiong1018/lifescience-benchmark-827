# report.md

## 条件

- 条件：C0
- 未使用 Skill
- 协议 API：2.16
- 模块：Magnetic Module Gen2

## 样本与转移

- 样本数：24
- 每样本净转移阶段：8
- 总净转移行数：192

## 吸头

- 移液器：P300 single，范围 20–300 µL
- 两盒吸头总数：192
- SOP 需求：144
- 结论：吸头数量足够

## 体积与孔容量

- lysis：80 µL
- beads：120 µL
- remove：250 µL
- wash 1 add/remove：180 µL / 180 µL
- wash 2 add/remove：180 µL / 180 µL
- elution：40 µL
- 最大孔体积：250 µL
- 结论：按给定事实，体积未超过最大孔体积

## 试剂

- 扣除 dead volume 后，各试剂均足够。

## 模拟结果

固定命令：

python -m opentrons.simulate output/protocol.py

实际返回：

- exit code：1
- 错误：ModuleNotFoundError: No module named 'opentrons'

## 结论

模拟失败。协议 not execution-ready。不得声称模拟成功或协议可执行。