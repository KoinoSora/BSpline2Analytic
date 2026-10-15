# 独立验证器位置误差界优化

基线为 `d345b30ea7aaa1725f4de140c0c7893af3dacfea`。本轮只修改独立 `ValidatePairedConversion` 的位置误差界及其目标分片准备，并添加测试和计时。共享 `IntervalBounds.cpp`、源面提取、原投影验证器、分类和优化算法、一般兜底及正式转换流程均未修改。保持严格浮点、向外舍入、正权重和原容差要求。

## 实现与认证含义

新增私有 [`PairedPositionBounds.cpp`](../../src/PairedPositionBounds.cpp)。对公共分片上的多项式源面 `S` 和有理目标 `N/W`，仍先构造误差分子 `E = WS - N`。

1. 先利用 Bernstein 系数凸包求保守上界。设三个坐标的系数绝对值上界为 `Mx, My, Mz`，权重系数最小下界为 `wmin > 0`，则整个分片满足 `||E/W|| <= sqrt(Mx² + My² + Mz²) / wmin`。严格小于容差才通过；失败只表示这个上界不够紧，继续原有有理控制系数求界及必要的递归细分。无效或非正权重不能由快速界通过。
2. 对齐次数采用缓存的直接升阶系数，沿 U/V 分开变换，代替逐级升阶。常量方向直接复制系数。误差相减复用乘积缓冲区，紧界求解中已同次的坐标直接引用，三个坐标共享权重倒数。
3. 三个坐标共享乘法系数；由解析构造得到的 `W(u,v)=Wu(u)Wv(v)` 保留两个一维因子，乘法沿 U/V 分开进行。反向、交换和分段映射同步变换因子；不通过浮点矩阵的近似秩猜测可分离性。通用二维乘法仍作为后备路径。
4. 只有接触外边界的公共分片才裁剪目标法向量。四条外边界的切向及法向认证仍全部执行。

快速界是覆盖连续域的保守上界，不是采样判断。返回的位置上界允许比旧版更宽；这不表示真实几何误差增大。算法仍要求完整覆盖和严格小于容差。同一新版的源面缓存命中与不缓存结果保持完全一致。

`PairedValidationReport` 新增 `fastBoundPasses`、`tightBoundChecks` 和 `surfaceComposeMs`、`surfaceProductMs`、`surfaceDifferenceMs`、`surfaceFastBoundMs`、`surfaceTightBoundMs`。五个耗时是 `certificate.timing.surfaceBounds` 的内部子项，不能再次加到验证总计上；正常完成的节点满足 `fastBoundPasses + tightBoundChecks == surfaceNodes`。节点计数累计所有映射尝试与细分节点。

## 同批性能对比

环境为 Windows x64、MSVC 19.42、SGK 1.7.11、Release，区间文件采用 `/fp:strict`。使用 40 个普通非有理球面的真实分类首候选，每例每种模式测三次，轮换执行顺序。分类、来源证书、预热、结果检查和日志不计入验证耗时；不运行优化器。未固定 CPU 核，不保证机器独占。

比较程序从 Git 提取基线的四个私有文件，重命名私有命名空间和公开验证函数；保留旧版数值实现，仅省去比较中未调用的旧版公开缓存包装，以避免重复符号。两版使用同一个未修改的基础区间库及相同编译配置。下表全部来自最终版本的同批比较，不能与上一轮报告的绝对时间直接相除。

| 平均耗时，ms/例 | 修改前独立验证器 | 新版首次验证 | 新版缓存命中 |
| --- | ---: | ---: | ---: |
| 整次验证，外部计时 | 16.757 | **8.849** | **5.577** |
| 准备阶段 | 3.862 | 3.877 | 0.620 |
| 位置误差界 | 10.322 | **2.433** | 2.459 |
| 四边检查 | 2.491 | 2.457 | 2.416 |

位置误差界加速 **4.24 倍**，耗时减少 **76.43%**；首次整次验证加速 **1.89 倍**。首次填充源面缓存的整次验证另计，平均 **9.034 ms**；空上下文构造不在验证计时内。缓存命中性能只代表后续复验，不能代替首次成本，也不代表整个优化器的加速倍数。

新版首次验证的位置误差界内部平均耗时：

| 操作 | ms/例 |
| --- | ---: |
| 源分片裁剪、目标复合及必要的法向量裁剪 | 0.433 |
| 权重乘法 | 1.203 |
| 次数对齐与误差分子相减 | 0.774 |
| 快速连续上界 | 0.020 |
| 较紧的控制系数上界 | 0.000 |

三轮共 1875 个分片节点全部由快速界通过，40 例均只需一次映射尝试，无位置细分。位置最大上界为 `1.652644492e-6`，仍严格小于 `1e-5`；切向、法向最大上界分别为 `0.000799916961°`、`0.001052695512°`，低于 `0.1°`。新旧位置上界不要求逐位相同；同一新版的三轮结果及缓存/不缓存映射、误差界和工作量完全一致。

原始数据：[comparison.csv](paired-position-2026-10-15/comparison.csv)；自动汇总：[comparison.summary.json](paired-position-2026-10-15/comparison.summary.json)。基线源码摘要和当前库摘要均已记录。

## 正确性与覆盖

四项定向 CTest 全部通过；单元测试单独运行，其余三项在同一覆盖批次运行。

| 测试 | 结果 |
| --- | --- |
| `PairedValidation` | **72783 项断言通过**。包含独立逐级升阶参照、已知多项式系数、五类母面的八种参数方向、权重裁剪、通用与可分离乘法交叉核对，以及原有负例和缓存测试。 |
| `PairedSphere` | 50 个非有理球面全部读取；40 个普通案例通过，10 个退化案例保持原分流状态。 |
| `PairedDataset500` | 全部 500 例读取并核验来源证书；280 个普通非有理案例通过。 |
| `Regression` | 原精确路径 **3847 项断言通过**。 |

新增端到端专项验证了：快速界过宽但紧界可通过；两种界均过宽但细分后可通过；距离为 `9.9e-6` 时通过、等于 `1e-5` 或为 `1.01e-5` 时不通过。另有非正权重、内点凸起、边界切向/法向超限、分段映射、交换/反向、节点预算等负例与边界情况。

完整 500 例的状态与 `d345b30` 逐例一致：280 `Passed`、206 `Unsupported`、9 `ClassificationFailed`、5 `NoCandidate`。`Unsupported` 包含 200 个有理输入和 6 个非有理退化案例；接口仍只支持非有理输入。280 个成功案例共 11211 个位置节点由快速界通过，三项最大上界分别为 `3.649806305e-6`、`0.003627390382°`、`0.004193016408°`，均低于原容差。它们的边界切向和法向上界与基线逐例完全相同。

覆盖 CSV：[sphere.csv](paired-position-2026-10-15/sphere.csv)、[dataset500.csv](paired-position-2026-10-15/dataset500.csv)；审计：[coverage-summary.json](paired-position-2026-10-15/coverage-summary.json)；日志：[unit.log](paired-position-2026-10-15/unit.log)、[coverage.log](paired-position-2026-10-15/coverage.log)。源代码、500 个输入及输出摘要见 [provenance.json](paired-position-2026-10-15/provenance.json)。CSV 汇总脚本通过历史、准备阶段和当前三种字段集合检查；复制/缩放时间仅用于解析器测试，不作为性能数据。

本轮未运行全仓库测试，也未处理历史 `Cone`、`Approximate` 回归测试的既有失败。性能结论限于上述独立验证器和测试输入。

## 复现

```powershell
cmake --build build-specialized-timing --config Release --target BSpline2AnalyticPairedValidationTests BSpline2AnalyticTests
ctest --test-dir build-specialized-timing -C Release -R '^BSpline2Analytic\.(PairedValidation|PairedSphere|PairedDataset500|Regression)$' --output-on-failure
python tests/BenchmarkPairedPreparation.py --baseline d345b30 --build-dir build-specialized-timing --output comparison.csv --runs 3
```

脚本沿用准备阶段的比较工具，现也支持准备阶段优化后的基线。临时基线代码和二进制放在忽略的构建目录，不切换分支，不修改正式源码。后续将独立验证器接入正式转换流程属于另一项集成工作。
