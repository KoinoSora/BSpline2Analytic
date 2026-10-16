# 等参优先、失败兜底：非极端全量测试（2026-10-16）

新增 `SequentialConversion.h::ConvertIsoThenFallback`，完成 1,160 个固定输入、每例三次，共 3,480 次测试。所有输入均先进入等参专用路径，失败后进入一般兜底；没有向算法传入来源母面、参数图或预期标签。985 个普通正例中，920 个三次均成功且独立输出审计通过。无解负例没有被误接受，但全部 100 例均因预算退出，不能算作正确拒绝。**本轮回归未全部通过。**

## 正确性及 10 ms 判断

“小于 10 ms”以每例三次实测整体耗时的最大值严格小于 10 ms 判定。算法耗时包含整个 API 的执行与返回，不含日志、CSV、来源认证或输出审计。

| 数据集 | 输入数 | 普通正例通过 | 无解负例正确拒绝 | 诊断数据成功且审计通过 | 全部成功且审计通过 | 三次均 <10 ms | 单次最大 / ms |
|---|---:|---:|---:|---:|---:|---:|---:|
| 等参专用 Dataset500 | 500 | 460/460 | — | 40/40 | 500 | 485/500 | 24.0191 |
| 新版兜底 main | 600 | 440/500 | 0/100，均预算失败 | — | 440 | 425/600 | 46.4167 |
| 新版兜底 extra | 60 | 20/25 | — | 13/35 | 33 | 37/60 | 46.4460 |
| 合计 | 1,160 | 920/985 | 0/100 | 53/75 | 973 | 947/1,160 | 46.4460 |

全部 3,480 次真实墙钟均不超过 50 ms，真实超时为 0。183 例、549 次因协作预算退出而判计算失败；提前退出不意味着几何无解。成功且审计通过并三次均小于 10 ms 的共有 **943 例**；只统计普通正负例的“正确且 <10 ms”为 894 例。所有成功输出的独立审计均通过，误接受已认证无解负例为 0。

| 数据集 | 进入兜底的案例数 | 总预算失败案例数 | 三次算法累计 / ms | 平均每次 / ms |
|---|---:|---:|---:|---:|
| 等参专用 | 43 | 0 | 5,042.6815 | 3.3618 |
| 兜底 main | 440 | 160 | 27,263.9692 | 15.1466 |
| 兜底 extra | 53 | 23 | 3,489.8333 | 19.3880 |
| 合计 | 536 | 183 | 35,796.4840 | 10.2863 |

536 个进入兜底的案例中，349 个最终成功且审计通过。等参阶段自身触及子预算共 21 次，这些尝试仍可继续进入兜底，不直接视为总预算失败。上述累计只加算法实测值，不加进程启动、数据加载、来源认证、输出审计或结果输出时间。

## 尚未通过的输入

- 普通正例共 65 个预算失败，均为球面：main 的 `204–215`、`224–235`、`244–255`、`264–275`、`284–295`，以及 extra 的 `30–34`。main 中平面、圆柱、圆锥、圆环各 100/100 成功，球面为 40/100。
- main 的 100 个无解负例均预算失败，没有错误成功，也没有得到及时、完整的拒绝结果。
- 75 个诊断输入中，53 个成功，18 个预算失败，4 个非超时拒绝（extra `0–3`）。这些结果保留原标签，不混入普通正确率。
- 549 次总预算退出中，523 次停在兜底参数优化，16 次停在 Jacobian，8 次在候选切换，分类边界评分及边界验收各 1 次。这里描述的是退出位置，不把它当成几何失败根因。

完整清单：[逐例结果（含中文结果列）](results/sequential-non-extreme-20261016/case-results.csv)、[未通过与超时案例](results/sequential-non-extreme-20261016/failures.csv)、[未达到 10 ms 的 213 例](results/sequential-non-extreme-20261016/over-10ms.csv)、[逐次预算退出](results/sequential-non-extreme-20261016/timeouts.csv)。

## 串联与预算实现

1. 有理输入先调用精确识别；非有理输入使用两子区间起步的增量等参搜索，再调用当前母面验证器，保留失败候选重试。
2. 等参成功立即返回。失败则执行一般候选分类，调用当前有理边界验证器；候选未认证时复用母面残差优化器，以当前样本进行一轮母面参数优化，改善后重新验收，再继续其他候选。
3. 两条路径共用一次起点和 50 ms 总期限；等参有 20 ms 子预算，为兜底留出时间。沿用 `ExecutionDeadline` 的 10% 清理余量（最多 5 ms），因此默认总预算在约 45 ms 协作取消、等参约 18 ms 取消。没有重新给兜底分配 50 ms，也没有截断记录的耗时。
4. 将取消检查接入等参采样、搜索、精确母面界、母面位置细分、边界位置及联合切向检查；取消异常不会被误吞成普通几何错误。兜底分类、验证和优化使用同一剩余预算。
5. API 内完成时以及测试外层 API 返回后都检查真实墙钟；超过 50 ms 强制失败。标准 Windows 环境的协作取消不保证进程必定在第 50 ms 被抢占，实际超出仍必须原样记录并判失败。本轮没有发生真实超出。

没有调用已废弃的 `ValidateConversion` 或 `ValidateGeneralHybrid`，也没有运行旧 general 数据或任何极端/压力输入。仅编译旧实现不构成执行旧验证器。

等参非有理路径保持 `1e-5` 母面/边界位置约束及严格 `0.1°` U/V 切向联合预算；兜底保持内部增强采样、连续边界位置及沿边切向认证，不构造曲面法向认证。两个分支不是相同的证书合同，兜底成功不等价于全域连续母面或双向 Hausdorff 证明。

## 数据与审计口径

覆盖 Dataset500 的全部 `case-*.b2a`，以及 `general-rational-v2/main` 600 例和 `extra` 60 例。独立 200/200/320 极端数据和 5 个压力输入不运行。V2 两方向次数 3–5、控制点不超过 200、几何容差 `1e-5` 的限制不变。

统一默认配置意味着 extra 中 `invalid_options` 和 `reject_budget` 的文件在本轮用作几何输入覆盖，没有注入历史非法选项或单节点预算；已知母面拒绝和非线性边界探针也没有向算法提供已知母面。它们与内部采样探针、Dataset500 的 40 个坍缩边输入共同单列诊断，不伪称复现历史特殊配置断言。尤其“偏离已知母面”不代表不存在另一可接受母面。

来源认证在计时前完成。成功输出的独立位置、边界和切向审计在计时后完成；审计可以否定成功输出，但不能把失败转为成功，也不反馈到候选生成。每次计时使用新上下文，无逐例预热。原始 CSV 保留 API 内时间与外部完整调用时间，最终汇总使用后者 `total_ms`；子阶段时间是其组成部分，不重复相加。

运行环境为 Windows x64、MSVC 19.42、Release，同一工作站串行逐例子进程，未固定 CPU 亲和性、未独占机器。批次期间另有一次 CMake 配置和两项很短的路由单测，未隐藏可能的系统调度影响；本报告保留所有实测值，不挑选最快值，也不将本次数据宣传为无干扰性能保证。

## 验证及复现

构建成功；`MotherDistanceUnit`、`GeneralRationalBoundaryUnit`、`SequentialShortCircuit`、`SequentialFallback` 均通过。路由单测覆盖成功直接返回、强制等参失败后进入兜底、非法/零总预算、极小预算取消且不发布成功输出。全量运行及记录复核无执行/完整性错误；运行器返回 **1**，准确表示存在预算失败，并非全部回归通过。

```powershell
cmake -S . -B build-sequential-20261016 -A x64
cmake --build build-sequential-20261016 --config Release --target BSpline2AnalyticSequentialTests BSpline2AnalyticMotherDistanceTests BSpline2AnalyticGeneralRationalBoundaryTests -- /m:1 /nr:false
ctest --test-dir build-sequential-20261016 -C Release -R 'BSpline2Analytic\.(SequentialShortCircuit|SequentialFallback|MotherDistanceUnit|GeneralRationalBoundaryUnit)$' --output-on-failure --no-tests=error
python tests/RunSequentialNonExtreme.py --exe build-sequential-20261016/Release/BSpline2AnalyticSequentialTests.exe --sdk ../../../SGGK --data tests/data --output build-sequential-20261016/full --runs 3
```

输出目录必须全新；`--resume` 只接受相同可执行文件、输入哈希、重复次数和配置，已完整的案例只复核记录、不重跑算法。

证据：[原始 3,480 行](results/sequential-non-extreme-20261016/runs.csv)、[逐例汇总](results/sequential-non-extreme-20261016/cases.csv)、[分组汇总](results/sequential-non-extreme-20261016/summary.csv)、[输入与可执行文件 SHA256](results/sequential-non-extreme-20261016/run.json)、[源码哈希及构建信息](results/sequential-non-extreme-20261016/provenance.json)、[路由单测日志](results/sequential-non-extreme-20261016/sequential-unit.log)。SDK、许可证和二进制不入库。
