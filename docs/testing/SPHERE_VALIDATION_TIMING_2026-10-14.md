# 非有理球面验证五步计时（2026-10-14）

本轮仅给非有理等参验证增加分步计时，未调整距离公式、区间运算、容差、预算或验收条件。基线为 `ff17321`。仅运行固定 Dataset500 的 50 个非有理球面（奇数 ID 301–399），包含 40 个普通案例及 10 个退化案例；没有运行其他类型或有理路径测试。

**主要结果：40 个普通案例的验证平均为 433.034 ms/例，边界验证占 85.15%。每轮的 625 个曲面小片和 564 个边界段均在初始分片上完成判定，没有额外细分。**

## 五步计时定义

`ConversionValidationReport::timing` 保存本次验证的五项互斥墙钟耗时；`ConversionStatistics::initialValidationTiming` 和 `refinementValidationTiming` 分别累计所有初始候选、所有优化后验证，避免覆盖前面的候选或轮次。CSV 新增 `initial_*`、`refined_*` 两组列和 `validation_timing_schema=exclusive-validation-v1`，原 `specialized-v2` 四项阶段口径保持不变。

| 步骤 | 字段后缀 | 计时范围 |
| --- | --- | --- |
| 1. 检查与准备 | prepare | 浮点环境/输入/候选检查、ExtractCells、BoundaryMaps、局部坐标角点、边界端点认证及初始化辅助数据 |
| 2. 曲面距离上界 | surface_bounds | 每个曲面节点的局部坐标变换、整片 SurfaceDistance、中心点求值及其距离界 |
| 3. 曲面判定与细分 | subdivide | 曲面遍历、预算/容差/深度分支、Split、递归调度；排除第 2、5 步 |
| 4. 边界验证 | boundary | 边界坐标变换、偏导/切向/法向构造、整段及中心点 BoundEdge、单调性判定与边界细分；排除单独计时的报告更新 |
| 5. 报告整理 | report | 约束初始状态、采样上界汇总、Record 中的状态/误差界/见证区域记录、异常处理、最终消息及正常退出清理 |

第 3 步只包含曲面递归；边界递归属于第 4 步。验证器生成反馈信息，实际 AddWitnessSamples 补点和求解仍属于优化计算，不计入第 5 步。部分边界见证位置的局部准备及对应关系判定留在第 4 步。

计时器切换活动步骤时暂停外层时间，递归不重复累计；五项之和等于内部 `totalMs`。内部总计不包含验证调用前后的返回赋值、统计累加等微小外围时间，因此略小于原 `classify_validate_ms + optimize_validate_ms`。本批每轮平均差值为 0.0374 ms。计时器自身开销包含在测量中，微秒量级步骤只用于判断数量级，不作为精确优化收益。

## 测量方法与结果

Windows x64、MSVC 19.42.34435.0、SGK 1.7.11、Release，独立构建目录 `build-specialized-timing`。三轮串行测试，每轮每例正式转换一次，首个普通球面预热一次。SDK 初始化、读文件、来源重建、独立输入证书、预热与输出日志不计入转换耗时。未绑定 CPU 或独占机器；本轮是分步诊断，不与上轮墙钟时间直接解释为加速或退化。

下表“每轮累计”为三轮累计除以 3；占比按各步骤累计除以验证内部总计计算，未对逐例百分比简单平均。单位为毫秒。

| 步骤 | 全部 50 例每轮累计 | 全部占比 | 普通 40 例每轮累计 | 普通案例均值 | 普通占比 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1. 检查与准备 | 1130.1837 | 6.5090% | 1088.1272 | 27.2032 | 6.2820% |
| 2. 曲面距离上界 | 1482.2252 | 8.5365% | 1482.2252 | 37.0556 | 8.5572% |
| 3. 曲面判定与细分 | 0.2349 | 0.0014% | 0.2349 | 0.0059 | 0.0014% |
| 4. 边界验证 | 14749.2048 | 84.9441% | 14749.2048 | 368.7301 | 85.1504% |
| 5. 报告整理 | 1.5660 | 0.0090% | 1.5612 | 0.0390 | 0.0090% |
| 验证合计 | 17363.4145 | 100% | 17321.3532 | 433.0338 | 100% |

40 个普通案例每轮各验证一次，全部在初始候选验证通过；优化验证及优化迭代均为零。10 个退化案例中，每轮只有 3 个进入验证，均止于准备阶段，累计内部验证平均 42.0613 ms；其余 7 个未进入验证，五步耗时均为零。三个进入验证的退化案例随后在优化采样准备失败，每轮累计优化计算平均 0.4877 ms，仍无优化迭代和优化后验证。

| 测量轮次 | 普通 40 例验证累计（ms） | 边界验证占比 |
| --- | ---: | ---: |
| 1 | 17121.6401 | 85.4775% |
| 2 | 18023.6467 | 84.8500% |
| 3 | 16818.7728 | 85.1394% |

| 控制点档 | 普通案例数 | 三轮平均验证耗时（ms/例） | 边界验证占比 |
| --- | ---: | ---: | ---: |
| 50–79 | 8 | 106.9852 | 95.2377% |
| 80–109 | 8 | 301.3968 | 89.8461% |
| 110–139 | 8 | 440.7717 | 87.4881% |
| 140–169 | 8 | 542.0866 | 84.4196% |
| 170–200 | 8 | 773.9287 | 81.1079% |

三轮平均最慢案例为 ID 383（196 个控制点）：验证 956.3855 ms，其中准备 76.3575 ms、曲面距离界 123.0650 ms、曲面调度 0.0253 ms、边界验证 756.8482 ms、报告整理 0.0895 ms。

## 与代码实现对应的结论

每轮普通案例的曲面节点数为 625，与输入非零节点矩形数之和相等；边界访问数为 1189 − 625 = 564，与初始边界段数相等。不存在预算耗尽，所有普通案例成功，因此本批没有额外二分；第 3 步微小的非零时间来自遍历和条件判定。

当前优先检查对象是 `ConversionValidation.cpp` 的 boundary 分支及 `BoundEdge`：每个边界段重新进行局部变换和二维求导，构造边界切向、叉积法向，然后分别对整段和中心点计算距离及两个方向夹角。其内部调用多项式点积/叉积、乘法及区间反三角函数。本轮已定位到这一整体模块，尚未测量这些内部操作各自占比，不能将 85% 全归因于某一种运算。

## 正确性与修改范围

- Release 的 `BSpline2AnalyticDataset500Tests` 构建成功，日志未发现编译警告/错误。
- 三轮每轮完整输出 50 行，所有来源证书和固定数据检查通过；50 个输入 SHA-256 与上一轮一致。
- 五项耗时守恒、非负值、初始/优化后调用次数对应、内部不超过外层计时均通过检查；本批覆盖正常验证、准备阶段异常和未进入验证的返回路径。
- 与上轮完整 500 例记录逐例核对状态、几何界、切向界、法向界、访问片数、验证调用次数及优化迭代数，三轮差异均为 0。
- 每轮状态为 Success 40、Degenerate 1、NoCandidate 5、Unresolved 1、InvalidInput 3。退化案例仍预期 Degenerate，故存在 9 个既有状态异常：317、319、337、339、357、359、377、379、397。三个 CTest 退出码均为 8，汇总脚本退出码为 1，不能宣称全部通过。
- 未运行其他曲面测试、完整 500 例、Debug 或必须优化的专项；本批未覆盖新增计时在实际细分、优化后验证上的非零执行。
- 生产代码仅变更近似验证计时、近似转换统计累计及公开统计结构。未改有理验证、一般验证、数值公式、数据及容差。新增公开成员后调用方需重新编译并链接。

## 复现和证据

从仓库根目录执行（已配置 SDK 的 Release 构建）：

```powershell
cmake --build build-specialized-timing --config Release --target BSpline2AnalyticDataset500Tests --parallel 1
ctest --test-dir build-specialized-timing -C Release -R '^BSpline2Analytic.SpherePolynomialDataset$' --output-on-failure
```

每轮结束后将 `build-specialized-timing/dataset500-sphere-polynomial-Release.csv` 分别保存为 `run1.csv`、`run2.csv`、`run3.csv`，然后执行：

```powershell
python tests/SummarizeSphereValidationTiming.py run1.csv run2.csv run3.csv --baseline docs/testing/specialized500-2026-10-14/results.csv --output summary
```

本机继承环境中同时存在 PATH/Path，构建与 CTest 由 Python 启动并以 `{k.upper(): v for k, v in os.environ.items()}` 规范化环境键，避免 MSBuild 重复键错误。预期状态异常导致非零退出时 CSV 仍完整写出，汇总会保留该失败结论。

- 三轮逐例数据：[run1.csv](sphere-validation-timing-2026-10-14/run1.csv)、[run2.csv](sphere-validation-timing-2026-10-14/run2.csv)、[run3.csv](sphere-validation-timing-2026-10-14/run3.csv)。
- [逐例三轮平均](sphere-validation-timing-2026-10-14/case-means.csv)、[完整汇总](sphere-validation-timing-2026-10-14/summary.json)。
- [环境及源码/二进制指纹](sphere-validation-timing-2026-10-14/provenance.json)、[50 个输入指纹](sphere-validation-timing-2026-10-14/input-sha256.csv)。
- 测试日志：[run1.log](sphere-validation-timing-2026-10-14/run1.log)、[run2.log](sphere-validation-timing-2026-10-14/run2.log)、[run3.log](sphere-validation-timing-2026-10-14/run3.log)。
