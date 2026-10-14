# 独立参数对应验证器：实现范围与测试结果

本轮新增独立入口 `ValidatePairedConversion`，实现“解析面片有理化 → 参数对应 → 连续误差界”的方案，并单独接入 CTest。基线提交为 `4f351c0fb82111891ebb0ae43ae423f5b9dfe8ba`，构建环境为 Windows x64、MSVC 19.42、Release、SGK 1.7.11。

## 代码范围和调用

| 文件 | 本轮职责 |
| --- | --- |
| [`include/PairedConversionValidation.h`](../../include/PairedConversionValidation.h) | 独立公开接口、状态、预算、映射和计时报告。 |
| [`src/PairedRationalChart.cpp`](../../src/PairedRationalChart.cpp) | 解析面片的正权重有理 Bézier 分片、区间裁剪、局部坐标变换及私有 Bernstein 乘积系数缓存。 |
| [`src/PairedConversionValidation.cpp`](../../src/PairedConversionValidation.cpp) | 仿射/分段单调映射、共同分区、位置与边界角度连续认证。 |
| [`src/PairedValidationInternal.h`](../../src/PairedValidationInternal.h) | 新模块内部数据结构。 |
| `src/ApproximateGeometry.cpp`、`src/ApproximateInternal.h` | 保存分类采样的源边参数；把原 `RecoverApproxPatch` 的边界参数展开提取为可复用函数，恢复逻辑保持原有数值行为。 |
| `tests/PairedValidationTests.cpp`、`tests/SummarizePairedValidation.py`、`CMakeLists.txt` | 专项、固定案例对照、CSV 审计及测试注册。 |

新模块复用基础区间算术、Bernstein 提取、输入检查和边界参数展开；不调用旧 `ValidateConversion`，不调用优化器，也不自动回退。现有 `ConvertApproximate`、`Recognize`、一般兜底、优化目标和最终发布流程均未替换。测试程序为了配对比较和构造优化后候选，会显式调用旧验证器及优化器。

```cpp
#include "PairedConversionValidation.h"

bspline2analytic::PairedValidationOptions options;
options.validation.distanceTolerance = 1e-5;
options.validation.boundaryTangentToleranceDeg = 0.1;
options.validation.boundaryNormalToleranceDeg = 0.1;
auto report = bspline2analytic::ValidatePairedConversion(input, candidate, options);
if (report.Passed()) {
    auto bound = report.certificate.geometry.upperBound;
    auto timeMs = report.certificate.timing.totalMs;
}
```

`candidate` 使用分类恢复或优化后的**当前**几何参数、矩形域、U/V 交换/反向及法向符号；验证期间不重新恢复并覆盖这个域。每次调用重新提取源面，避免跨候选缓存失效。单次调用内部复用提取结果和乘积系数。

## 实现要点和证书含义

1. 平面、圆柱、圆锥、球和规则环形圆环解析构造成正权重有理 Bézier 分片；角跨度按不大于四分之一圆拆分。这是精确有理表示的区间包围，不是采样拟合。当前不额外物化一个 SDK NURBS 对象。
2. 先尝试源参数到有理图归一化参数的仿射对应；失败后利用与 `RecoverApproxPatch` 相同的边界参数展开，构造两对相对边共同支持的分段线性单调函数 `f(u)`、`g(v)`。先将解析角参数换算为真正的有理弧参数，再形成映射，不能直接把角度当作 NURBS 参数。映射采样只用来提议对应关系。
3. 将源节点、映射折点和有理图分片边界拉回源域形成共同分区。每格内映射为仿射，区间 de Casteljau 裁剪给出复合有理面 `R=A/W`。
4. 对非有理源面 `P` 构造 `E=P*W-A`，统一 Bernstein 次数。正权重下，`E/W` 的范数上界为各有理控制向量范数上界的最大值。无法直接通过时对子格细分；中点区间下界只能证明**当前对应**超差，不能证明不存在其他对应。
5. 四条外边分别比较源切向与有理导数分子 `A'W-AW'`；源法向 `P_u × P_v` 与解析法向分子的方向同时检查。要求点积严格正，并用叉积/点积比的连续上界判断角度。所有位置及四边角度均通过才返回成功。

映射连续、严格单调、覆盖整个候选矩形，因此位置界覆盖整个成对参数域；不声称该映射就是最近点投影。边界的切向和法向按此对应检查，和旧投影验证器的配对语义不同。映射折点及源内部节点采用分片单侧导数检查。

`NotCertified`、`BudgetExceeded`、`Unsupported` 都不表示几何上不存在可行解。两条相对边需要不同的速度函数时，可分离映射可能过于保守；专项中有一个旧验证器通过而新验证器返回 `NotCertified` 的平面案例，防止把这种失败误报为几何不可行。提前退出时未覆盖全域的上界置为无穷，不发布局部成功证书。

本版输入范围为连续、非有理、每方向次数 1～32 的 B 样条。球极点、锥顶、退化参数迹、多圈参数域、自交/尖点圆环及有理输入不在本版范围。默认每方向最多 512 个映射段、映射递归深度 12、每次认证最多 16384 个共同单元；还遵守原 `ValidationOptions` 的细分/节点预算。没有降低容差或采样放行分支。

## 非有理球面：三轮配对计时

使用固定 50 个非有理球面，每轮均重新加载并核验来源证书。对同一个真实分类结果的首个恢复候选，交替执行新旧验证器；每类预热一次，初始化、数据读取、来源证书、预热、额外采样核查及 CSV 写入均在验证计时之外。未固定 CPU 核或保证机器独占。

每轮 40 个普通案例全部通过两个验证器，10 个退化案例分为 `NoCandidate` 5、`ClassificationFailed` 2、`Unsupported` 3。三个运行的状态、误差界、映射次数及节点工作量完全一致。均未进入优化器。

| 范围 | 旧验证平均 ms/例 | 新验证平均 ms/例 | 倍率 |
| --- | ---: | ---: | ---: |
| 第 1 轮，40 个双方通过案例 | 339.647 | 32.591 | 10.42 |
| 第 2 轮，40 个双方通过案例 | 342.502 | 33.059 | 10.36 |
| 第 3 轮，40 个双方通过案例 | 339.978 | 32.594 | 10.43 |
| 三轮合计 120 次双方通过 | **340.709** | **32.748** | **10.40** |

验证耗时减少 **90.39%**。这个结果仅比较验证器，不代表整个优化流程有相同比例的加速。三轮平均分类计算为 1.393 ms/例；按“初值验证属于分类”的口径，本对照中的“分类计算 + 首候选验证”为旧 342.102、新 34.141 ms/例。它不是完整 `ConvertApproximate` 的端到端计时。优化计算/优化验证在本对照中未执行，不能据此给出优化器性能结论。

| 新验证的独占阶段 | 平均 ms/次 | 占验证内部总计 |
| --- | ---: | ---: |
| 准备、源提取、有理图及共同分区 | 21.767 | 66.47% |
| 位置误差界计算 | 8.807 | 26.89% |
| 位置递归细分 | 0.000 | 0.00% |
| 四边切向、法向及其递归检查 | 2.102 | 6.42% |
| 报告整理 | 0.071 | 0.22% |

五项为独占计时；`nurbs_ms`、`mapping_ms` 是准备阶段的子项，不能再次相加。有理图构造平均 0.107 ms。40 个普通球面全部在第一种仿射映射下通过，没有运行自适应分段映射，也没有位置细分；本批加速主要来自低次数有理误差式及边界一维计算，不能归因于分段映射。下一步热点是源提取及准备阶段。

最大位置上界 `1.410211624e-6`，最大切向角上界 `0.000799918°`，最大法向角上界 `0.001052697°`，均严格低于 `1e-5 / 0.1° / 0.1°`。

逐例结果和三轮均值：[原始数据目录](paired-validation-2026-10-14/)、[sphere-case-means.csv](paired-validation-2026-10-14/sphere-case-means.csv)、[summary.json](paired-validation-2026-10-14/summary.json)。

## 完整 Dataset500 分流与专项

读取全部 500 个原始输入，核验来源证书。280 个普通非有理案例均由真实分类首候选通过新旧两种验证器；200 个有理输入显式检查新接口返回 `Unsupported`；其余 20 个非有理退化输入分别为 6 个 `Unsupported`、9 个分类失败、5 个无候选。因此这不是“新接口认证 500/500”。

| 来源类型 | 双方通过数 | 旧验证平均 ms | 新验证平均 ms | 倍率 |
| --- | ---: | ---: | ---: | ---: |
| 平面 | 100 | 34.301 | 6.634 | 5.17 |
| 圆柱 | 50 | 153.398 | 15.120 | 10.15 |
| 圆锥 | 40 | 159.374 | 15.970 | 9.98 |
| 球 | 40 | 339.486 | 32.451 | 10.46 |
| 圆环 | 50 | 451.321 | 31.486 | 14.33 |
| 合计 | **280** | **191.502** | **17.609** | **10.88** |

该表为一轮全量诊断；三轮球面统计另表给出。全部 280 例都在第一种映射下通过。最大三项上界为 `3.130945186e-6`、`0.003627392°`、`0.004193019°`。

专项 **604 项断言通过**，覆盖：精确平面、独立的位置/切向/法向超差、非线性参数平面的 512 段映射、7 种交换/反向组合、相对边速度不兼容、候选几何与域改变、实际优化后再次验证、反法向、无效容差、预算耗尽、不支持的舍入环境、五类解析图及多角度分片，以及跨周期接缝的多分片圆柱完整证书。非线性参数平面先仿射失败，再分段映射通过。

最终 CTest 再次执行 `PairedSphere`、`PairedDataset500`，要求分别认证 40/280 个普通案例，防止仅检查已成功案例而遗漏通过率回归。该次结果保存在 `sphere-final.csv`、`coverage-final.csv` 和 `coverage-final.log`；它们不混入上述性能统计。

## 旧路径回归与既有失败

- 原精确回归：3847 项断言通过。
- 原两条专用接口重跑 500 例：460 `Success`、25 `Degenerate`、4 `Unresolved`、6 `InvalidInput`、5 `NoCandidate`，保留既有 15 个异常预期结果，CTest 退出码为 8。与原 `specialized500-2026-10-14/results.csv` 逐例比较，状态、三项误差界、优化迭代数、原始候选数、恢复数、验证节点数共八项均无差异。
- 原 `Cone` 和 `Approximate` 测试仍失败，分别为 `Classification cache or deferred recovery regressed`、`sphere-antipodal-meridians extended chart failed`。从修改前提交单独 `git archive` 并在另一构建目录重新编译，这两个相同失败均复现。基线日志随本报告保存，未为使本轮通过而修改原断言或放宽标准。

因此新增三个定向测试通过，但不能声称全仓库测试全部通过。上述旧问题留在原有模块范围内。

## 复现

在仓库根目录执行，`SGGK_SDK_DIR` 指向本机授权 SDK：

```powershell
cmake -S . -B build-paired -G 'Visual Studio 17 2022' -A x64 -DSGGK_SDK_DIR='D:/YourSdk/SGGK'
cmake --build build-paired --config Release --target BSpline2AnalyticPairedValidationTests
ctest --test-dir build-paired -C Release -R '^BSpline2Analytic\.Paired(Validation|Sphere|Dataset500)$' --output-on-failure
```

CTest 设置 SDK 工作目录和 DLL 搜索路径。直接执行程序时需自行设置这两项；入口为 `--unit`、`--sphere <data-dir> <csv>`、`--dataset <data-dir> <csv>`。三轮对比时将球面输出分别保存为 `sphere-run1.csv` 至 `sphere-run3.csv`，全量输出为 `dataset500.csv`。

```powershell
python tests/SummarizePairedValidation.py docs/testing/paired-validation-2026-10-14 --baseline docs/testing/specialized500-2026-10-14/results.csv
```

数据、源文件及测试程序摘要见 [provenance.json](paired-validation-2026-10-14/provenance.json)。SDK、许可证和编译产物不提交。
