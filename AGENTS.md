# BSpline2Analytic 项目约定

处理本项目的代码修改时，读取并遵守 [bspline2analytic-git-sync Skill](.agents/skills/bspline2analytic-git-sync/SKILL.md)。

用户已持续授权：每次完成一轮代码修改及相关测试、构建或配套文档修改后，执行相应验证，提交到当前工作分支并推送至本项目 Git 远端，核验推送结果后再汇报完成。正常提交与推送无需重复询问；用户后续明确指示优先。

仓库为 https://github.com/KoinoSora/BSpline2Analytic，保持私有。SDK、许可证、凭据及构建产物不进入仓库。完整流程和失败处理以所链接的 Skill 为准。

## 一般兜底验证器约束

2026-10-15 用户明确要求：旧解析 UV 直线边界验证器 `ValidateGeneralHybrid` 废弃，保留实现和历史测试代码，但今后不再运行其单测、数据集测试或性能对照。历史报告和测量文件仅供查阅，不按其中旧命令重跑。当前开发和测试使用 `ValidateGeneralRationalBoundaries`；其测试也不得调用旧验证器作为对照。共享类型和边界提取代码仍可复用。本约束适用于上述一般兜底旧验证器；后续用户明确修改约束时再调整。

2026-10-16 用户明确取消当前一般兜底验证器 `ValidateGeneralRationalBoundaries` 的曲面法向夹角约束：不构造法向认证数据，不检查或报告法向夹角。保留内部采样距离、连续边界位置和边界曲线切向检查；上游分类所需的几何法向不属于本约束。

## 验证器约束（2026-10-15）

等参专用路径的旧投影验证器 `ValidateConversion`（`src/ConversionValidation.cpp`）已废弃。保留其实现、接口和历史测试证据，但今后不得直接或间接运行其测试、预热或性能对照；尚调用它的 `ConvertApproximate` 完整转换测试也暂停。新验证器使用 `ValidatePairedConversion` / `paired::ValidatePrepared`。此约束不表示正式转换入口已经完成迁移。

精确有理路径 `Recognize` 不属于上述旧验证器，保留测试；新验证器尚不支持有理源面时，应单列精确识别结果，其五模块耗时标为不适用，不能当作零耗时通过。数据来源证书、输出独立审计在算法计时外保留。一般路径的验证器沿用其单独约定。

第三个独立实验验证器为 `ValidateMotherDistanceConversion` / `mother::ValidatePrepared`：全域（包括边界）采用带浮点保护的母面距离界，保留四边区间对应位置检查。2026-10-16 最新要求：同时检查边界上的 U/V 参数切向，以同一子区间上的 `sin²(U夹角上界)+sin²(V夹角上界) < sin²(0.1°)` 联合预算验收；独立容差为 `jointTangentToleranceDeg`。不构造源面法向、不转入直接法向认证，不检查相邻面 G1。联合条件超限、无法认证、预算耗尽须分别明确说明，不能将联合条件未通过称为实际法向误差超限。U/V 和联合推导角度上界分别报告，不能混用历史沿边切向 CSV 字段。此变更不适用于其他验证器。允许与当前 `ValidatePairedConversion` 对照，但本轮用户要求修改后暂不运行测试；两个验收合同不同，不能混称为同一证书。正式转换入口尚未迁移，旧验证器禁测约束不变。

## 文件分类

公开接口放在 `include/`，算法实现和内部头文件放在 `src/`，示例放在 `examples/`，测试代码与固定数据放在 `tests/`。算法、测试和需求文档分别放在 `docs/algorithm/`、`docs/testing/`、`docs/requirements/`；根目录保留构建入口与导航。路径调整时同步更新 CMake 和相对文档链接。

项目可复用 Skill 集中保存于 [.agents/skills/](.agents/skills/README.md)。编写算法说明书时，使用 [algorithm-documentation Skill](.agents/skills/algorithm-documentation/SKILL.md)。
