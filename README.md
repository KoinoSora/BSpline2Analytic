# BSpline2Analytic

基于海星 SGK 1.7.11 的四边样条面片解析识别与转换模块，使用 C++17。提供精确 NURBS 识别、非有理等参线近似转换及一般情况兜底。

统一入口为 `SurfaceConversion.h` 中的 `Convert`：首先尝试适用的专用模块，失败后进入 `ConvertGeneral`；`attempts` 记录尝试、跳过及分类失败原因。精确有理模块和非有理专用模块仍按输入表示选择，两个原接口继续保留。边界退化、非直线或圆弧时专用分类失败，不表示整张曲面退化。共享源采样上下文按需构造边界、Greville 初值/评分集及优化样本。

## 仓库目录

| 目录/入口 | 内容 |
| --- | --- |
| [`include/`](include/) | 公开接口头文件，含 `GeneralConversion.h`。 |
| [`src/`](src/) | 算法实现与内部头文件。 |
| [`examples/`](examples/) | 使用统一转换接口的示例。 |
| [`tests/`](tests/README.md) | 回归、专项、性能测试及汇总脚本；固定输入在 `tests/data/`。 |
| [`docs/algorithm/`](docs/algorithm/README.md) | 算法 PDF、LaTeX 源文件、设计与实现说明。 |
| [`docs/testing/`](docs/testing/TEST_REPORT.md) | 测试报告、固定测试集说明与批次报告。 |
| [`docs/requirements/`](docs/requirements/REQUIREMENTS.md) | 需求说明及原始需求 PDF。 |
| [`.agents/skills/`](.agents/skills/README.md) | 项目 Git 同步与算法文档编写 Skill。 |
| `CMakeLists.txt`、`Build.ps1` | 根目录构建入口。 |

[文档索引](docs/README.md)列出各说明书。`build/` 仅保存本机构建和运行结果，不进入 Git。

## 当前功能

新增统一入口 [`ConvertIsoThenFallback`](include/SequentialConversion.h)：每个输入先尝试等参专用路径，失败后进入一般兜底候选验证与母面参数优化，验收使用当前母面/有理边界验证器。默认等参子预算 20 ms，两分支共享 50 ms 总预算；预算耗尽或真实墙钟超过上限均判失败。`Passed()` 为最终成功判据，`path`、`isoAttempted`、`fallbackAttempted` 和两个预算标志保留执行情况。该入口不调用下文仍保留的旧 `Convert` / `ConvertApproximate` 验证链。测试命令与计时口径见[统一非极端测试](tests/README.md#等参优先失败兜底的统一测试)。

```cpp
#include "SequentialConversion.h"
auto result = bspline2analytic::ConvertIsoThenFallback(inputSurface);
if (result.Passed()) {
    // iso_exact: result.exact；iso_mother: result.isoPatch / result.mother；
    // general_rational: result.fallback.patch / result.fallback。
    // 两条路径的验收合同不同，使用对应报告，不混称同一种全域证书。
}
// result.totalMs 包括分类、失败尝试、优化和验证，不包括调用方结果输出。
```

首次准备进行母面优化时，串联入口额外尝试一次五类母面的整体距离排除证书。证书从输入曲面实际区间点值构造，不使用案例标签；只有距离下界严格超过当前容差才提前拒绝。局部预算为 8 ms，计入原有 50 ms 总预算，未能证明则继续原优化流程。可用 `!result.timeBudgetExceeded && result.exclusion.Excluded()` 识别已证明的负例；这类结果没有输出几何，`Passed()` 仍为假。实现边界和 100 个负例的复测结果见[排除证书测试报告](docs/testing/NEGATIVE_EXCLUSION_2026-10-19.md)。

输入来自平面、圆柱面、圆锥面、球面或圆环面的精确解析转换，四条边界对应解析面参数域的两对等参线。允许节点插入、升阶、U/V 交换及参数反向；识别不依赖来源标签或固定控制点排列。

- `AnalyzeBoundaries`：提取四条边界，恢复直线或圆弧的几何参数和有向范围。
- `Recognize`：检查边界关系，恢复 SGK 解析曲面对象、矩形参数范围、U/V 方向对应和验证报告。
- 默认距离容差为 `1e-5`，单位跟随输入模型，通常为 mm。
- 默认进行边界和母面的 Bézier/Bernstein 误差界估计；节点区间采样用于快速检查。
- 失败结果区分无效输入、退化、病态、边界关系不一致、歧义、验证失败和不支持。

模块包含静态库、示例和测试工具。2026-10-16 一般兜底数据已重构为 [rational-domain-v2](docs/testing/GENERAL_RATIONAL_DATASET_V2.md)：600 主案例、60 专项、320 极端和 5 压力输入，正例约束定义于母面 NURBS 参数域。两方向次数均为 3–5，控制顶点总数不超过 200，算法几何容差保持 `1e-5`。原一般兜底数据已废弃并保留归档，旧测试入口关闭；`dataset500/case-*.b2a` 专用路径数据不受影响。本轮只生成数据、独立认证来源并编译入口，暂不运行兜底模块测试。以下旧一般兜底性能报告均为历史结果，不适用于新版输入。

候选评分前去重的最新结果见 [评分去重报告](docs/testing/GENERAL_POSITIVE500_DEDUP_REPORT_2026-10-13.md)，包含本轮重新测得的 500 个正例基线、修改后结果与实际评分次数。

8 次初值尝试与分阶段评分的历史对比见 [分类加速报告](docs/testing/GENERAL_POSITIVE500_STAGED_REPORT_2026-10-14.md)，当时 8 个圆环因没有圆环候选而失败。最新[圆环初值修复报告](docs/testing/GENERAL_POSITIVE500_TORUS_REPORT_2026-10-14.md)改用四法线直接求轴及双候选比较，在相同 8 次预算下达到 **500/500**，原失败 8 例全部修复；分类累计 0.419 秒，转换累计 37.097 秒。相对历史批次分类约减少 76.76%，该比较并非同进程配对基准。

非有理输入使用新增的 `ApproximateConversion.h`：`ConvertApproximate` 通过四边与少量内部等参线恢复五类候选，再进行参数优化和独立验收。默认距离严格小于 `1e-5`、边界切向及边界曲面法向角分别严格小于 `0.1°`。等参专用路径采用 `specialized-v2` 计时口径：`timing.classification` 包含分类计算与初始候选验证，`timing.optimization` 包含优化计算与优化后验证；四项明细为 `classificationCompute`、`classificationValidation`、`optimizationCompute`、`optimizationValidation`。初值合格时两项优化计时为零。独立专用接口的 `timing.validation` 为零；统一入口追加 PCurve 验证和一般兜底仍使用各自的计时。明细是阶段总计的组成部分，不能再次累加。默认启用解析 Jacobian（非光滑分支局部差分回退），可通过 `options.refinement.jacobianMode=JacobianMode::FiniteDifference` 切换全差分对照。算法及验证口径见 [非有理转换说明](docs/algorithm/APPROXIMATE_CONVERSION.md)。

独立实验验证器 `ValidatePairedConversion` 通过候选的有限矩形域、分段单调参数对应和正权重有理 Bernstein 误差界，验证非有理输入的位置、四边切向及法向。入口为 [`PairedConversionValidation.h`](include/PairedConversionValidation.h)，目前单独接入测试，现有转换流程仍使用原验证器。实现范围、状态含义与配对计时见[独立验证器测试报告](docs/testing/PAIRED_VALIDATION_2026-10-14.md)。

第三个独立入口 [`ValidateMotherDistanceConversion`](include/MotherDistanceValidation.h) 用带浮点保护的母面距离界检查整个参数域（包括边界），保留四边区间对应位置检查。2026-10-16 改为同时检查边界上的 U/V 参数切向，以同一子区间上的 `sin²(U夹角上界)+sin²(V夹角上界) < sin²(0.1°)` 联合预算验收；容差为独立的 `jointTangentToleranceDeg`。目标母面的两个参数方向解析正交，交换、反向和分段单调映射保持此性质。分别报告 U/V 上界和联合推导角度上界，不构造源面法向，不进行直接法向认证或失败回退，也不检查相邻面 G1。联合条件超限与无法认证分别报告 `JointTangentConditionExceeded`、`JointTangentConditionUnresolved`，节点预算耗尽报告 `JointTangentBudgetExceeded`；这些状态不能解释为实际法向角已超限。联合预算版本已完成 [16 例定向测试](docs/testing/MOTHER_JOINT_SELECTED_2026-10-16.md)，160 次验证与独立审计全部通过；9/16 例的分类加首次验证均值低于 10 ms，尚未重跑全 500 例。它不提供内部逐参数配对或有限面片双向 Hausdorff 证明。球面/圆环历史测试见[历史报告](docs/testing/MOTHER_DISTANCE_SPHERE_TORUS_2026-10-15.md)。当前仍独立接入测试，未替换正式转换流程。

母面验证器最新采用无样条叉乘的正比例切向残差认证，位置与联合角度容差不变。[16 例复测报告](docs/testing/MOTHER_RESIDUAL_SELECTED_2026-10-16.md)中，新版本 160 次验证与独立审计全部通过；同机比较的球面/圆环角度检查耗时分别下降 37.0%/35.9%。分类加首次验证仍有 4 个大圆环均值超过 10 ms。

该入口现使用逐片双精度 Bézier 提取、误差包络和独立源缓存；代表性慢例的改前/改后测量及高精度包络审计见[流式提取报告](docs/testing/MOTHER_STREAMING_2026-10-16.md)。

2026-10-15 起，等参旧验证器 `ValidateConversion` 已废弃，保留实现和历史证据但不再测试；调用它的完整转换回归暂停。`Build.ps1 -Dataset500Only` / `-SpherePolynomialOnly` 改为新验证器测试，有理输入另外测试精确 `Recognize`。当前命令、五模块计时和未覆盖范围见[测试说明](tests/README.md#等参专用路径计时)。下面涉及近似、平台、Jacobian、圆锥和极端完整转换的旧测试开关已禁用，不能按历史命令重跑。

独立验证器现已采用 U/V 两步提取和逐控制点坐标变换；可传入 `PairedValidationContext` 复用同一源面的准备结果。源面内容改变会自动失效，候选改变仍重新验证。首次验证及缓存命中的同批对比见[准备阶段优化报告](docs/testing/PAIRED_PREPARATION_2026-10-15.md)。

位置误差界现采用快速连续保守界、缓存直接升阶与可分离权重乘法；快速界不够紧时仍执行控制系数求界和必要的细分，容差保持不变。40 个普通非有理球面的同批位置计算加速 4.24 倍，详见[位置误差界优化报告](docs/testing/PAIRED_POSITION_2026-10-15.md)。

当前一般兜底验证器为独立入口 [`ValidateGeneralRationalBoundaries`](include/GeneralRationalBoundaryValidation.h)：线性约束定义于母面 NURBS 参数域，输出解析母面参数及四条三维 NURBS 边界，不输出解析 UV PCurve。内部保持增强 Greville 采样，边界采用有理误差样条位置/切向连续认证，并覆盖 SDK 输出系数，尚未替换正式入口。旧 `ValidateGeneralHybrid` 已废弃，保留代码但不再测试，也不再作为对照运行。范围和历史慢例测试见[新边界验证报告](docs/testing/GENERAL_RATIONAL_BOUNDARY_2026-10-15.md)，当前命令见[测试说明](tests/README.md)。

## 构建与运行

需要 Windows x64、Visual Studio 2022 C++ 工具链、CMake 3.22 及以上、SGK 1.7.11 SDK 和有效本地许可证。`cmake`、`ctest` 需在 PATH 中。

在本模块目录运行：

```powershell
# SDK 默认路径：本模块所在目录的 ../../../SGGK
.\Build.ps1 -Configuration Release
.\Build.ps1 -Configuration Debug

# 单独克隆本仓库时，指定 SDK 根目录（其中含 include、x64-win）
.\Build.ps1 -SdkDirectory 'D:\YourSdk\SGGK' -Configuration Release -Benchmark

# 一般兜底旧 GeneralOnly / GeneralPositiveOnly 开关已废弃。
# 新版数据生成和认证命令见 tests/README.md；本轮暂不运行兜底测试。
```

脚本默认构建库、示例及测试程序，并运行精确、近似、平台、Jacobian/五类扰动、Dataset500 及两组圆锥专项回归。`-Benchmark` 额外运行原有 1000 面片精确识别测试，记录写入 `build/benchmark-Release.csv`。`-ApproximateBenchmark` 运行 35 次近似转换计时（含优化成功和禁用优化的失败分支），记录写入 `build/approximate-benchmark-Release.csv`。数据生成、预热和日志不计入转换耗时。可用 `-BuildJobs` 指定构建并发数，默认 1。`-JacobianOnly` 单独运行解析导数及五类扰动对照测试，输出 `build/perturbation-Release.csv` 和附加 `.jacobian.csv`；`-PlateauOnly` 保留平台专项入口；`-Dataset500Only` 单独加载已保存的 500 个面片并输出 `build/dataset500-Release.csv`。`-ConeOnly` 只构建并运行 100 个圆锥固定案例及圆锥专项检查，输出 `build/dataset500-cone-Release.csv` 和 `build/cone249-Release.csv`，保留原 500 例结果。专项开关互斥。

直接运行程序时，需将 SDK 对应配置的 DLL 目录加入 PATH：Release 为 `x64-win/bin`，Debug 为 `x64-win/bind`，RelWithDebInfo 为 `x64-win/binw`。当前测试环境在 SDK 根目录加载许可证；测试脚本自动设置工作目录和 DLL 搜索路径。RelWithDebInfo 的目录映射已配置，本轮没有构建测试该配置。

模块也可通过 `add_subdirectory` 接入宿主工程，并链接 `BSpline2Analytic`。CMake 自动向调用方传递 `include/` 路径，原有 `#include "SurfaceConversion.h"` 等写法不变；手动集成时将头文件搜索目录设为模块的 `include/`。`src/` 只供库实现和内部测试使用。应显式传入 `SGGK_SDK_DIR`；现有 samples 顶层工程还需使用实际 SDK 的库名前后缀配置。独立模块构建是本轮验证过的入口。

## 接口示例

统一调用：

```cpp
#include "SurfaceConversion.h"

bspline2analytic::SurfaceConversionOptions options;
options.rational.distanceTolerance = 1e-5;
options.polynomial.validation.distanceTolerance = 1e-5;
options.general.fitting.validation.distanceTolerance = 1e-5;
auto result = bspline2analytic::Convert(inputSurface, options);
if (result.Succeeded()) {
    const auto& timing = result.Timing();
    const auto& finitePatch = *result.patch;
    auto analytic = finitePatch.MotherSurface();
    const auto& fourPCurves = finitePatch.PCurves();
    auto boundaryBound = result.pcurveValidation.boundaryHausdorffUpperBound;
    // fourPCurves[k].sdkCurve 是母面 SDK 参数空间中的有界一次样条。
    // start/end 是模块规范参数空间中的端点，两种参数约定不能混用。
}
```

`path` 标识最终尝试分支，`details` 为三种结果之一，保留该模块报告；`Timing()` 汇总所有尝试。统一入口在专用模块成功后同样恢复并验证线性 PCurve，成功才发布 `result.patch`。母面对象本身不包含裁剪；有限域由四条 PCurve 围成的简单参数四边形定义，允许退化为三角形。失败统一返回 `ConversionFailed`，并保留模块及候选的诊断信息。原独立接口仍遵循各自历史验收合同。

## 已废弃的独立混合验证器

`GeneralHybridValidation.h::ValidateGeneralHybrid` 于 2026-10-15 废弃。实现、共享类型和历史测试代码保留；其 CMake 测试目标和 CTest 注册已移除，旧批量/对照脚本停止执行。后续不运行该验证器的单测、数据集或性能比较，新验证器测试也不调用它。旧[独立混合验证器报告](docs/testing/GENERAL_HYBRID_VALIDATION_2026-10-15.md)仅作为历史记录查阅。

## 一般情况与线性 PCurve

`ConvertGeneral` 支持有理（正权重）与非有理连续样条，采用平面 3、球 4、圆柱 5、圆锥 6、圆环 7 个母面自由度生成少量候选。Greville 参数按完整带重数节点计算；缺少内部覆盖时补点。内部点交错划分为生成集和独立评分集，边界另行评估；无效法向不生成初值。球与圆环评分接近且边界表现接近时先尝试球，失败后保留独立圆环候选。圆环初值直接求解四条点法向直线的公共交线，每组最多构造两个候选；子午圆拟合补充初值点集中的最多 12 个分散点，经中心化与缩放后用 QR 求解，再交给独立点集评分。秩亏、无实轴线或非法半径的构造会被拒绝；固定 8 次尝试不提供任意输入上的搜索完备性保证。

求解器复用原阻尼最小二乘、约束修正、信赖域和停滞处理，通过 `RefinementProblem` 替换几何评价；一般母面导数当前使用中心差分及单侧回退。边界见证点反馈到母面距离优化；PCurve 恢复和有限面片验证在每轮后重新执行。边界不满足线性要求时可能失败，局部优化失败不表示不存在可行解。

默认每类尝试 8 组初值（失败也占一次，不按去重后的数量补抽）；候选去重后，从独立评分集选取分散的 16 个内部点粗评。每类每批保留最多 3 个完整评分候选，只有粗评最大距离已严格超过本类第三名完整分数时才暂缓补评；同分候选继续完成评分。补评复用已算过的内部点，再评估边界。首批全部转换失败后继续取后备候选，耗尽后返回失败。`seedAttemptsPerType` 可显式调整，`coarseScorePointCount=0` 可启用完整评分对照；最终连续验收容差不受这些选项影响。

四个源参数角点只投影一次，相邻边复用端点。投影参数沿边连续展开周期分支，生成二维直线段，不要求它是等参线；原边界参数到直线段的对应可非线性。球面另外尝试边界支持的坐标系。极点、锥顶及退化边使用邻边参数极限协调角度，不把相同三维点的不同参数角点直接合并。

采样初筛检查偏离、越界及回走。最终使用严格浮点模式的向外舍入区间运算与有理 Bernstein 提取，验证源面到有限域的距离，以及源边界与线性 PCurve 三维像的双向 Hausdorff 上界。边界位移、投影像偏离和端点覆盖都使用最终几何容差；参数距离不能直接套用长度容差。达到预算或遇到无法解析的参数奇异性时返回失败，不降级为采样验收。此实现不提供整个有限面片的双向 Hausdorff 或全局可逆重参数化证明。

`GeneralValidationReport::pcurveFailure` 区分 `NonlinearPCurve`、`PCurveOutOfRange`、`NonMonotonePCurve`、周期分支不一致、退化参数未解决、无效边界环、几何超差和连续验证未解决。`GeneralConversionResult::attempts` 保留各候选信息。一般圆环当前限定为规则环形圆环（大半径大于小半径）；自交、尖点及 Lemon 类型不会被该兜底路径接受。

非有理样条：

```cpp
#include "ApproximateConversion.h"
auto result = bspline2analytic::ConvertApproximate(inputSurface);
if (result.status == bspline2analytic::ConversionStatus::Success) {
    auto analytic = result.patch->surface;
    auto domain = result.patch->sdkDomain;
    auto distanceBound = result.validation.geometry.upperBound;
    auto classifyMs = result.timing.classification.elapsedMs;
    auto optimizeMs = result.timing.optimization.elapsedMs;
    auto validateMs = result.timing.validation.elapsedMs;
}
// 失败时 patch 可能仅是诊断候选，不能作为验收合格结果。
```

精确有理路径：

宿主程序按 SGK 要求管理 `sggk::init()` / `sggk::fini()`，确保 SDK 对象在 `fini()` 前销毁。

```cpp
#include "BSpline2Analytic.h"

bspline2analytic::RecognitionOptions options;
options.distanceTolerance = 1e-5;
auto result = bspline2analytic::Recognize(inputSurface, options);

if (result.status == bspline2analytic::Status::Success) {
    auto analytic = result.surface;       // 原生 SGK 解析曲面
    auto domain = *result.domain;         // 展开周期后的矩形参数域
    auto type = result.parameters.type;
    auto bound = result.validation.surfaceDistanceBound;
} else {
    // 使用 result.status 和 result.message 处理失败。
}
```

四边顺序为 `Umin, Umax, Vmin, Vmax`。`mapping.analyticCorners` 的顺序为 `(umin,vmin), (umax,vmin), (umin,vmax), (umax,vmax)`。`swappedUV` 表示输入 U 对应解析 V；`reversedU/V` 分别描述输入 U/V 增大时，对应解析参数是否减小。这些是方向信息，不能作为 NURBS 参数到解析角度的线性映射；点的对应可通过解析对象的 `CalcParam` 求取，并按返回域处理周期分支。

返回的原点、轴向和参数零点采用等价局部坐标系，不要求重现原设计坐标系。`surface` 是解析母面，`domain` 是独立的面片范围；接口不创建或替换拓扑面。

## 验证及适用限制

`BernsteinBound` 表示对输入面片到解析母面的全参数域距离进行 Bernstein 系数界估计，并加入浮点保护项。它不是使用定向舍入实现的严格区间算术证书，也不是裁剪面片之间的双向 Hausdorff 证明。边界范围和方向通过采样追踪、角度展开及一致性检查恢复；不声称已证明任意重参数化的全域单调性。

设置 `requireBound=false` 会省略母面的 Bernstein 界验证，结果明确标为 `Sampled`，`surfaceDistanceBound` 保持无穷大；边界本身仍执行误差界检查。达到细分深度或工作量限制时返回失败，不将采样通过冒充全域验证通过。

普通非退化面片是当前主范围。周期接缝与单周范围已覆盖；球极点、锥顶等坍缩边返回 `Degenerate`。圆锥按同一锥瓣恢复；Apple/Lemon 圆环已有规则区域测试，但跨奇异点、跨锥顶或多周范围未作为支持承诺。高度病态的极短圆弧可能报告无法恢复或无法验证。

上述 `BernsteinBound` 和 `requireBound` 说明属于原 `Recognize` 接口。新增近似路径使用向外舍入区间/Bernstein 连续界，验证有限解析面片距离及边界切法向，并提供 `Passed/Exceeded/Unresolved`；没有允许采样代替验收的选项。其法向比较仅沿边界，参数及 SGK 对象需保持同步。

当前不处理非等参裁剪、任意自由曲面的可转换性搜索和拓扑重建。浅锥、短圆弧、奇异边界仍没有近似路径的专项保证；相邻面 G1、内存泄漏审计和下游运算加速待后续专项验证。新验证器目前耗时明显高于纯类型判断，不能把类型判断时间当作完整转换时间。

## 文档及仓库约定

- [当前实现算法说明书（10 页 PDF）](docs/algorithm/BSpline2Analytic-algorithm.pdf)：总体与子模块流程图、数学推导、真实案例和工程化边界；[可编辑源文件、编译与证据清单](docs/algorithm/README.md)。
- [需求说明](docs/requirements/REQUIREMENTS.md)
- [算法与代码实施方案](docs/algorithm/IMPLEMENTATION_PLAN.md)
- [测试报告](docs/testing/TEST_REPORT.md)
- [500 个固定面片的数据格式、生成与测试](docs/testing/DATASET500.md)
- [非有理转换、参数对应与三阶段计时](docs/algorithm/APPROXIMATE_CONVERSION.md)
- [原始需求 PDF](docs/requirements/中船海舟样条转二次曲面.pdf)
- [开发与 Git 同步约定](AGENTS.md)

本仓库仅管理模块代码与文档。SGK SDK、运行库、许可证及构建产物不上传。SGK 及原始需求材料的权利归各自权利人；本项目不另行授予第三方材料的使用或分发许可。

## 极端场景专项

当前专项使用[业务尺寸极端集 200 例](docs/testing/BUSINESS_EXTREME200.md)，主要半径与长度在 0.01–1000，包含短弧、点边、大坐标与参数表示变化；[旧极限压力集](docs/testing/EXTREME200.md)保留作历史参考。运行 `.\Build.ps1 -Configuration Release -ExtremeBoundaryOnly` 只执行新数据及 7 个模块诊断，需要 Python 3.8 或以上。单例进程隔离，外部看护上限 300 秒；超时不伪造阶段耗时。结果记录分类、优化、验证，以及每组 `max(分类+优化)`。
新矩阵结果见[业务尺寸极端测试报告](docs/testing/BUSINESS_EXTREME200_REPORT_2026-10-14.md)。旧数据结果见[历史极端测试报告](docs/testing/EXTREME200_REPORT_2026-10-13.md)；两批编号或均值不直接比较。

等参专用路径的新四项计时与固定 500 例复测见[2026-10-14 报告](docs/testing/SPECIALIZED500_TIMING_REPORT_2026-10-14.md)：包含分类计算、分类验证、优化计算及优化验证，以及退化状态失败和未修改基线核对。

非有理验证新增五步互斥计时，50 个非有理球面三轮专项测量见[验证分步计时报告](docs/testing/SPHERE_VALIDATION_TIMING_2026-10-14.md)。
