# 文档索引

| 分类 | 文档 | 用途 |
| --- | --- | --- |
| 算法 | [算法说明书（10 页 PDF）](algorithm/BSpline2Analytic-algorithm.pdf) | 当前实现的数学方法、流程图、真实案例及工程边界。 |
| 算法 | [源文件与编译说明](algorithm/README.md) | LaTeX 源文件、编译入口和测试证据清单。 |
| 算法 | [精确路径设计](algorithm/IMPLEMENTATION_PLAN.md) | 有理样条识别、解析母面重建及原始设计。 |
| 算法 | [非有理转换实现说明](algorithm/APPROXIMATE_CONVERSION.md) | 参数优化、Jacobian、连续验收与三阶段计时。 |
| 测试 | [测试报告](testing/TEST_REPORT.md) | 各轮回归和专项验证记录。 |
| 测试 | [500 例测试集](testing/DATASET500.md) | 数据来源、配额、格式、生成与执行方式。 |
| 测试 | [2026-10-13 全量重测](testing/DATASET500_RETEST_2026-10-13.md) | 固定批次统计、慢例与误差上界。 |
| 需求 | [需求说明](requirements/REQUIREMENTS.md) | 输入约束、验收目标与阶段范围。 |
| 需求 | [原始需求 PDF](requirements/中船海舟样条转二次曲面.pdf) | 原始项目材料。 |

代码与构建入口见[仓库 README](../README.md)，测试代码见 [tests/](../tests/README.md)，项目 Skill 见[Skill 索引](../.agents/skills/README.md)。历史报告保留其运行批次和基线，不以目录迁移更新历史性能结论。
