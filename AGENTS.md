# BSpline2Analytic 项目约定

处理本项目的代码修改时，读取并遵守 [bspline2analytic-git-sync Skill](.agents/skills/bspline2analytic-git-sync/SKILL.md)。

用户已持续授权：每次完成一轮代码修改及相关测试、构建或配套文档修改后，执行相应验证，提交到当前工作分支并推送至本项目 Git 远端，核验推送结果后再汇报完成。正常提交与推送无需重复询问；用户后续明确指示优先。

仓库为 https://github.com/KoinoSora/BSpline2Analytic，保持私有。SDK、许可证、凭据及构建产物不进入仓库。完整流程和失败处理以所链接的 Skill 为准。

## 一般兜底验证器约束

2026-10-15 用户明确要求：旧解析 UV 直线边界验证器 `ValidateGeneralHybrid` 废弃，保留实现和历史测试代码，但今后不再运行其单测、数据集测试或性能对照。历史报告和测量文件仅供查阅，不按其中旧命令重跑。当前开发和测试使用 `ValidateGeneralRationalBoundaries`；其测试也不得调用旧验证器作为对照。共享类型和边界提取代码仍可复用。本约束适用于上述一般兜底旧验证器；后续用户明确修改约束时再调整。

## 验证器约束（2026-10-15）

等参专用路径的旧投影验证器 `ValidateConversion`（`src/ConversionValidation.cpp`）已废弃。保留其实现、接口和历史测试证据，但今后不得直接或间接运行其测试、预热或性能对照；尚调用它的 `ConvertApproximate` 完整转换测试也暂停。新验证器使用 `ValidatePairedConversion` / `paired::ValidatePrepared`。此约束不表示正式转换入口已经完成迁移。

精确有理路径 `Recognize` 不属于上述旧验证器，保留测试；新验证器尚不支持有理源面时，应单列精确识别结果，其五模块耗时标为不适用，不能当作零耗时通过。数据来源证书、输出独立审计在算法计时外保留。一般路径的验证器沿用其单独约定。

## 文件分类

公开接口放在 `include/`，算法实现和内部头文件放在 `src/`，示例放在 `examples/`，测试代码与固定数据放在 `tests/`。算法、测试和需求文档分别放在 `docs/algorithm/`、`docs/testing/`、`docs/requirements/`；根目录保留构建入口与导航。路径调整时同步更新 CMake 和相对文档链接。

项目可复用 Skill 集中保存于 [.agents/skills/](.agents/skills/README.md)。编写算法说明书时，使用 [algorithm-documentation Skill](.agents/skills/algorithm-documentation/SKILL.md)。
