# 独立验证器准备阶段优化

基线为 `3985a1bc218ef7231cae6dc4a42a7fefb3d1aa82`，本轮只修改独立 `ValidatePairedConversion` 的源面准备与相关测试。共享 `IntervalBounds.cpp::ExtractCells`、原投影验证器、分类和优化算法、一般兜底及正式转换流程均未修改。保持原容差、向外舍入、严格浮点设置和连续证书要求。

## 实现

- 新增私有 [`PairedSourcePreparation.cpp`](../../src/PairedSourcePreparation.cpp)：每个有效 U/V 节点区间的提取矩阵各计算一次，控制网先沿 U、再沿 V 变换。同一 U 区间的中间结果复用于全部 V 区间。非均匀节点、重节点和源域端点的局部裁剪保留。
- [`PairedRationalChart.cpp`](../../src/PairedRationalChart.cpp) 每个候选只构造一次区间坐标架，直接变换各 Bernstein 控制点，消除平移常数的通用升阶及多项式临时对象。
- 新增可选 `PairedValidationContext`，保存世界坐标下的源 Bézier 分片。每次调用精确比较次数、节点、重数、源参数域和控制点；使用相同对象地址不能跳过这些比较。源数据改变会重新提取，候选几何、候选域、方向、映射和最终证书始终重新计算。
- 增加 `sourceExtractionMs`、`sourceCacheMs`、`localTransformMs`、`sourceCacheHits`、`uExtractionMatrices`、`vExtractionMatrices`。三个耗时均属于准备阶段内部子项，不能再次加到验证总计上。

直接调用原签名仍为独立的首次验证，不保留跨调用缓存。需要多候选或优化后复验时，可在调用方复用上下文：

```cpp
#include "PairedConversionValidation.h"

bspline2analytic::PairedValidationContext context;
auto first = bspline2analytic::ValidatePairedConversion(source, initialCandidate, context);
// 调用方更新候选；本轮没有更换正式转换流程中的验证器。
auto next = bspline2analytic::ValidatePairedConversion(source, refinedCandidate, context);
```

上下文仅保存一个源面，可移动、不可复制；不同并发调用方使用各自的上下文。构造上下文不进行几何准备，首次调用中完成提取并计时。缓存命中仍检查输入、预算及浮点环境。源码变换与乘加顺序调整可能改变区间上界的末尾数字，不要求与旧版逐位相等；同一新版的缓存/不缓存结果要求完全一致。

## 同批三模式比较

构建环境为 Windows x64、MSVC 19.42、SGK 1.7.11、Release，区间相关文件采用 `/fp:strict`。比较程序从 Git 提取基线的三个私有文件，只重命名命名空间和公开函数以避免符号冲突，并增加提取/坐标变换两个子计时器。基线与新版链接同一个未修改的基础区间库。

使用 40 个普通非有理球面的真实分类首候选，每例三次测量，轮换三种模式的执行顺序。来源证书、分类、预热、结果审计和日志不计入验证耗时；不运行优化器。未固定 CPU 核，也不保证机器独占。以下数据全部来自同批比较，不与上一日报告的绝对耗时直接相除。

| 平均耗时，ms/例 | 修改前独立验证器 | 新版首次验证 | 新版缓存命中 |
| --- | ---: | ---: | ---: |
| 整个验证，外部计时 | 41.194 | **17.217** | **14.191** |
| 准备阶段 | 27.658 | **3.912** | **0.677** |
| 源面提取 | 24.880 | **3.281** | 0.000 |
| 局部坐标变换 | 2.562 | **0.417** | 0.438 |
| 缓存核对与复制 | 0.000 | 0.000 | 0.020 |
| 位置误差界计算 | 10.846 | 10.650 | 10.827 |
| 四边检查 | 2.600 | 2.568 | 2.601 |

首次验证的准备阶段加速 **7.07 倍**，耗时减少 **85.86%**；整个验证加速 **2.39 倍**。缓存命中相对修改前整个验证加速 **2.90 倍**。首次填充缓存的整次验证另外实测平均 **18.311 ms**，并未隐藏这次成本；上下文的空对象构造不在验证计时内。缓存命中性能只适用于同一源面后续复验，不能代表首次验证或整个优化流程。

| 源控制点数 | 修改前整次验证 ms | 新版首次验证 ms | 缓存命中 ms |
| --- | ---: | ---: | ---: |
| 64 | 3.216 | 2.741 | 2.383 |
| 100 | 21.406 | 9.968 | 7.957 |
| 121 | 37.660 | 15.997 | 13.141 |
| 144 | 50.946 | 20.298 | 16.964 |
| 196 | 92.744 | 37.083 | 30.509 |

每组 8 个案例，各重复三次。新版首次验证每次恰好提取一次；缓存命中每次提取次数为零。三轮结果稳定，同一新版缓存与不缓存的映射、误差界和工作量完全一致。三项最大上界分别为 `1.410211571e-6`、`0.000799916961°`、`0.001052695512°`，仍严格低于 `1e-5 / 0.1° / 0.1°`。

原始对比：[comparison.csv](paired-preparation-2026-10-15/comparison.csv)；自动汇总：[comparison.summary.json](paired-preparation-2026-10-15/comparison.summary.json)。源码摘要和库摘要记录在汇总及 [provenance.json](paired-preparation-2026-10-15/provenance.json) 中。

## 正确性与覆盖

本轮四项定向 CTest 全部通过：

| 测试 | 结果 |
| --- | --- |
| `PairedValidation` | **47491 项断言通过**。包含原负例、分段映射、实际优化后使用缓存复验，以及新旧提取系数包围交集、SDK 点值交叉检查。 |
| `PairedSphere` | 50 个非有理球面全部读取；40 个普通案例通过，10 个退化案例按原状态分流。 |
| `PairedDataset500` | 全部 500 例读取并核验来源；280 个普通非有理案例通过，200 个有理输入仍明确不支持，20 个非有理退化案例单独记录。 |
| `Regression` | 原精确路径 **3847 项断言通过**。 |

完整 500 例的新接口状态与基线逐例相同：280 `Passed`、206 `Unsupported`、9 `ClassificationFailed`、5 `NoCandidate`。其中 `Unsupported` 包括 200 个有理输入与 6 个非有理退化案例。没有把“接口按预期不支持”称为几何认证成功。

新增专项覆盖非均匀节点、重节点、U/V 不同次数、反向/交换、多分片和周期接缝；缓存测试覆盖控制点原地修改、节点插入、候选位移、候选域变化、降低预算、上下文移动及移动后的重新使用。源坐标变换后不污染缓存中的世界坐标分片。实际优化后的候选测试确认复验命中源面缓存，并重新获得证书。

覆盖结果见 [dataset500.csv](paired-preparation-2026-10-15/dataset500.csv)、[sphere.csv](paired-preparation-2026-10-15/sphere.csv)、[coverage-summary.json](paired-preparation-2026-10-15/coverage-summary.json) 和 [coverage.log](paired-preparation-2026-10-15/coverage.log)。本轮未运行全仓库测试，也未处理上一报告中原 `Cone`、`Approximate` 测试的既有失败。

## 复现

```powershell
cmake --build build-specialized-timing --config Release --target BSpline2AnalyticPairedValidationTests BSpline2AnalyticTests
ctest --test-dir build-specialized-timing -C Release -R '^BSpline2Analytic\.(PairedValidation|PairedSphere|PairedDataset500|Regression)$' --output-on-failure
python tests/BenchmarkPairedPreparation.py --baseline 3985a1b --build-dir build-specialized-timing --output comparison.csv --runs 3
```

基准脚本读取指定构建目录中的 SDK 配置，独立比较构建放在其 `paired-preparation-comparison/` 下，不切换工作分支、不修改正式源码。当前库需先完成 Release 构建。原始区间算术、共享源面提取及正式转换路径保持原样，后续切换正式验证入口属于另一项集成工作。
