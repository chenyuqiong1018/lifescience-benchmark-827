# report.md

## T0 条件自动选择与执行事实报告

### 候选选择
已从 183 个候选池中自动选择以下 5 个模块/任务：

- synthetic_biology_design
- protocol-generation-from-description
- protocol-to-executable-json
- unit_conversion_suite
- code_execution_analysis

### 数据规模与校验
- 样本数：24
- 阶段数：8
- 总行数：24 × 8 = 192
- 移液器：P300
- Tip 数量：144
- 试剂校验：通过
- 容量校验：通过

### 真实执行结果
固定模拟命令已真实执行，返回结果为：

- exit code: 1
- 错误信息：ModuleNotFoundError: No module named opentrons

### 结论
当前产物未通过真实运行环境验证，缺少 opentrons 依赖，不能视为可执行协议。

状态：not execution-ready。