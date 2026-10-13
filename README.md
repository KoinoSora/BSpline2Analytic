# BSpline2Analytic

基于海星 SGK 1.7.11 的四边样条面片解析识别与转换模块，使用 C++17。提供精确 NURBS 识别及非有理样条近似转换两条路径。

统一入口为 `SurfaceConversion.h` 中的 `Convert`：输入 `IsRational()` 为真时调用 `Recognize`，否则调用 `ConvertApproximate`。不使用案例来源标签分流，不在失败后静默切换路径。精确解析转换的平面通常权重为 1，因此也进入非有理分支。两条原接口继续保留。

## 当前功能

输入来自平面、圆柱面、圆锥面、球面或圆环面的精确解析转换，四条边界对应解析面参数域的两对等参线。允许节点插入、升阶、U/V 交换及参数反向；识别不依赖来源标签或固定控制点排列。

- `AnalyzeBoundaries`：提取四条边界，恢复直线或圆弧的几何参数和有向范围。
- `Recognize`：检查边界关系，恢复 SGK 解析曲面对象、矩形参数范围、U/V 方向对应和验证报告。
- 默认距离容差为 `1e-5`，单位跟随输入模型，通常为 mm。
- 默认进行边界和母面的 Bézier/Bernstein 误差界估计；节点区间采样用于快速检查。
- 失败结果区分无效输入、退化、病态、边界关系不一致、歧义、验证失败和不支持。

模块包含静态库、示例、回归测试和批量性能测试。本轮已验证 Release，历史 Debug 结果及各轮范围见 [测试报告](TEST_REPORT.md)。

非有理输入使用新增的 `ApproximateConversion.h`：`ConvertApproximate` 通过四边与少量内部等参线恢复五类候选，再进行参数优化和独立验收。默认距离严格小于 `1e-5`、边界切向及边界曲面法向角分别严格小于 `0.1°`。`timing.classification / optimization / validation` 分别累计类型判断、求解、验收耗时及调用次数，初值合格时优化计时为零。默认启用解析 Jacobian（非光滑分支局部差分回退），可通过 `options.refinement.jacobianMode=JacobianMode::FiniteDifference` 切换全差分对照。算法及验证口径见 [非有理转换说明](APPROXIMATE_CONVERSION.md)。

## 构建与运行

需要 Windows x64、Visual Studio 2022 C++ 工具链、CMake 3.22 及以上、SGK 1.7.11 SDK 和有效本地许可证。`cmake`、`ctest` 需在 PATH 中。

在本模块目录运行：

```powershell
# SDK 默认路径：本模块所在目录的 ../../../SGGK
.\Build.ps1 -Configuration Release
.\Build.ps1 -Configuration Debug

# 单独克隆本仓库时，指定 SDK 根目录（其中含 include、x64-win）
.\Build.ps1 -SdkDirectory 'D:\YourSdk\SGGK' -Configuration Release -Benchmark
```

脚本默认构建库、示例及测试程序，并运行精确、近似、平台、Jacobian/五类扰动、Dataset500 及两组圆锥专项回归。`-Benchmark` 额外运行原有 1000 面片精确识别测试，记录写入 `build/benchmark-Release.csv`。`-ApproximateBenchmark` 运行 35 次近似转换计时（含优化成功和禁用优化的失败分支），记录写入 `build/approximate-benchmark-Release.csv`。数据生成、预热和日志不计入转换耗时。可用 `-BuildJobs` 指定构建并发数，默认 1。`-JacobianOnly` 单独运行解析导数及五类扰动对照测试，输出 `build/perturbation-Release.csv` 和附加 `.jacobian.csv`；`-PlateauOnly` 保留平台专项入口；`-Dataset500Only` 单独加载已保存的 500 个面片并输出 `build/dataset500-Release.csv`。`-ConeOnly` 只构建并运行 100 个圆锥固定案例及圆锥专项检查，输出 `build/dataset500-cone-Release.csv` 和 `build/cone249-Release.csv`，保留原 500 例结果。专项开关互斥。

直接运行程序时，需将 SDK 对应配置的 DLL 目录加入 PATH：Release 为 `x64-win/bin`，Debug 为 `x64-win/bind`，RelWithDebInfo 为 `x64-win/binw`。当前测试环境在 SDK 根目录加载许可证；测试脚本自动设置工作目录和 DLL 搜索路径。RelWithDebInfo 的目录映射已配置，本轮没有构建测试该配置。

模块也可通过 `add_subdirectory` 接入宿主工程，并链接 `BSpline2Analytic`。应显式传入 `SGGK_SDK_DIR`；现有 samples 顶层工程还需使用实际 SDK 的库名前后缀配置。独立模块构建是本轮验证过的入口。

## 接口示例

统一调用：

```cpp
#include "SurfaceConversion.h"

bspline2analytic::SurfaceConversionOptions options;
options.rational.distanceTolerance = 1e-5;
options.polynomial.validation.distanceTolerance = 1e-5;
auto result = bspline2analytic::Convert(inputSurface, options);
if (result.Succeeded()) {
    const auto& timing = result.Timing();
    if (const auto* exact = std::get_if<bspline2analytic::RecognitionResult>(&result.details)) {
        auto analytic = exact->surface;
        auto domain = exact->domain;
    } else {
        const auto& approximate = std::get<bspline2analytic::ConversionResult>(result.details);
        auto analytic = approximate.patch->surface;
        auto domain = approximate.patch->sdkDomain;
    }
}
```

`path` 标识实际执行分支，`details` 保留原始状态及验证报告。成功语义依照被调用接口及其选项；统一入口不会将精确路径的母面 Bernstein 界升级为近似路径的有限面片区间证书。`RecognitionResult::timing` 同样提供三阶段时间，其中优化为零，分类包含边界识别和参数域恢复，验证记录母面候选的验证。

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
- [需求说明](REQUIREMENTS.md)
- [算法与代码实施方案](IMPLEMENTATION_PLAN.md)
- [测试报告](TEST_REPORT.md)
- [500 个固定面片的数据格式、生成与测试](DATASET500.md)
- [非有理转换、参数对应与三阶段计时](APPROXIMATE_CONVERSION.md)
- [原始需求 PDF](中船海舟样条转二次曲面.pdf)
- [开发与 Git 同步约定](AGENTS.md)

本仓库仅管理模块代码与文档。SGK SDK、运行库、许可证及构建产物不上传。SGK 及原始需求材料的权利归各自权利人；本项目不另行授予第三方材料的使用或分发许可。
