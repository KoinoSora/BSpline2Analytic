# BSpline2Analytic

基于海星 SGK 1.7.11 的四边 NURBS 面片解析识别与转换模块，使用 C++17。

## 当前功能

输入来自平面、圆柱面、圆锥面、球面或圆环面的精确解析转换，四条边界对应解析面参数域的两对等参线。允许节点插入、升阶、U/V 交换及参数反向；识别不依赖来源标签或固定控制点排列。

- `AnalyzeBoundaries`：提取四条边界，恢复直线或圆弧的几何参数和有向范围。
- `Recognize`：检查边界关系，恢复 SGK 解析曲面对象、矩形参数范围、U/V 方向对应和验证报告。
- 默认距离容差为 `1e-5`，单位跟随输入模型，通常为 mm。
- 默认进行边界和母面的 Bézier/Bernstein 误差界估计；节点区间采样用于快速检查。
- 失败结果区分无效输入、退化、病态、边界关系不一致、歧义、验证失败和不支持。

模块包含静态库、示例、回归测试和批量性能测试。Release 与 Debug 均已验证，结果见 [测试报告](TEST_REPORT.md)。

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

脚本构建 `BSpline2Analytic`、`BSpline2AnalyticSample`、`BSpline2AnalyticTests`，并运行回归测试。`-Benchmark` 额外运行 1000 个不同面片的测试，原始记录写入 `build/benchmark-Release.csv`。数据生成和变换不计入识别耗时，计时包含完整 `Recognize` 流程。

直接运行程序时，需将 SDK 对应配置的 DLL 目录加入 PATH：Release 为 `x64-win/bin`，Debug 为 `x64-win/bind`，RelWithDebInfo 为 `x64-win/binw`。当前测试环境在 SDK 根目录加载许可证；测试脚本自动设置工作目录和 DLL 搜索路径。RelWithDebInfo 的目录映射已配置，本轮没有构建测试该配置。

模块也可通过 `add_subdirectory` 接入宿主工程，并链接 `BSpline2Analytic`。应显式传入 `SGGK_SDK_DIR`；现有 samples 顶层工程还需使用实际 SDK 的库名前后缀配置。独立模块构建是本轮验证过的入口。

## 接口示例

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

当前不处理一般非有理样条的近似恢复、非等参裁剪和拓扑重建。相邻面 G1、内存泄漏审计、严格区间误差证明及下游运算加速仍待后续专项验证。

## 文档及仓库约定

- [需求说明](REQUIREMENTS.md)
- [算法与代码实施方案](IMPLEMENTATION_PLAN.md)
- [测试报告](TEST_REPORT.md)
- [原始需求 PDF](中船海舟样条转二次曲面.pdf)
- [开发与 Git 同步约定](AGENTS.md)

本仓库仅管理模块代码与文档。SGK SDK、运行库、许可证及构建产物不上传。SGK 及原始需求材料的权利归各自权利人；本项目不另行授予第三方材料的使用或分发许可。
