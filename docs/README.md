# 文档索引

| 分类 | 文档 | 用途 |
| --- | --- | --- |
| 算法 | [算法说明书（10 页 PDF）](algorithm/BSpline2Analytic-algorithm.pdf) | 当前实现的数学方法、流程图、真实案例及工程边界。 |
| 算法 | [源文件与编译说明](algorithm/README.md) | LaTeX 源文件、编译入口和测试证据清单。 |
| 算法 | [精确路径设计](algorithm/IMPLEMENTATION_PLAN.md) | 有理样条识别、解析母面重建及原始设计。 |
| 算法 | [非有理转换实现说明](algorithm/APPROXIMATE_CONVERSION.md) | 参数优化、Jacobian、连续验收与三阶段计时。 |
| 测试 | [测试报告](testing/TEST_REPORT.md) | 各轮回归和专项验证记录。 |
| 测试 | [500 例测试集](testing/DATASET500.md) | 数据来源、配额、格式、生成与执行方式。 |
| 测试 | [独立参数对应验证器](testing/PAIRED_VALIDATION_2026-10-14.md) | 分段单调映射、连续误差界、604 项断言与新旧验证器配对计时。 |
| 测试 | [独立验证器准备阶段优化](testing/PAIRED_PREPARATION_2026-10-15.md) | 分离提取、控制点坐标变换、源面缓存及同批三模式比较。 |
| 测试 | [独立验证器位置误差界优化](testing/PAIRED_POSITION_2026-10-15.md) | 快速连续上界、直接升阶、可分离权重乘法及 500 例覆盖。 |
| 测试 | [独立混合验证器](testing/GENERAL_HYBRID_VALIDATION_2026-10-15.md) | 内部增强 Greville 采样、四边位置与切向连续认证，600+60 例三轮测量及采样合同边界。 |
| 测试 | [2026-10-13 全量重测](testing/DATASET500_RETEST_2026-10-13.md) | 固定批次统计、慢例与误差上界。 |
| 测试 | [一般兜底正例修复报告](testing/GENERAL_POSITIVE500_FIX_REPORT_2026-10-13.md) | 边界参数方向、球参数图与奇异端点修复；500 个正例及性能对比。 |
| 测试 | [一般兜底评分前去重报告](testing/GENERAL_POSITIVE500_DEDUP_REPORT_2026-10-13.md) | 重复候选免于评分；500 个正例的前后比较及评分次数。 |
| 测试 | [一般兜底分类加速报告](testing/GENERAL_POSITIVE500_STAGED_REPORT_2026-10-14.md) | 8 次初值尝试、分阶段评分与后备候选；500 正例的三配置对比。 |
| 测试 | [一般兜底圆环初值修复报告](testing/GENERAL_POSITIVE500_TORUS_REPORT_2026-10-14.md) | 四法线直接求轴与双候选比较；8 次预算下 500/500，修复原 8 个圆环失败例。 |
| 测试 | [一般兜底极端 320 例](testing/GENERAL_EXTREME320.md) | 尺度 1e-2～1e3、螺旋带、三角片、大世界坐标；独立构建与认证。 |
| 测试 | [一般兜底极端 320 例报告](testing/GENERAL_EXTREME320_REPORT_2026-10-14.md) | 251 例通过、38 例返回失败、31 例转换超时；完整结果与分阶段耗时。 |
| 需求 | [需求说明](requirements/REQUIREMENTS.md) | 输入约束、验收目标与阶段范围。 |
| 需求 | [原始需求 PDF](requirements/中船海舟样条转二次曲面.pdf) | 原始项目材料。 |

代码与构建入口见[仓库 README](../README.md)，测试代码见 [tests/](../tests/README.md)，项目 Skill 见[Skill 索引](../.agents/skills/README.md)。历史报告保留其运行批次和基线，不以目录迁移更新历史性能结论。
业务尺寸极端场景：[200 例构造与运行说明](testing/BUSINESS_EXTREME200.md)、[2026-10-14 测试报告](testing/BUSINESS_EXTREME200_REPORT_2026-10-14.md)。
历史极限压力集：[旧 200 例构造说明](testing/EXTREME200.md)。
[2026-10-13 极端测试报告](testing/EXTREME200_REPORT_2026-10-13.md)保存 200 条逐例结果、三阶段分组统计和超时记录。
