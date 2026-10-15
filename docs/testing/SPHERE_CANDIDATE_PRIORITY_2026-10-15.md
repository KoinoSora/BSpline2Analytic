# 球面候选优先级调整与专项回归（2026-10-15）

## 修改范围

基线为 `61c3883daf0e0d5566520f759342693452b4ba73`，本轮仅调整等参专用分类器 `RecoverRankedCandidates` 的近球面圆环后置规则。存在最佳采样评分小于位置容差的球面候选时，将同时满足 `abs(majorRadius) <= 5 * distanceTolerance` 和 `abs(majorRadius) < radius` 的圆环类型后置。默认排序半径阈值为 `5e-5`，位置认证仍为 `1e-5`，切向及法向认证仍为 `0.1°`。候选数量、几何参数、类型内部顺序和后续候选回退保持不变。

`PairedValidationTests` 的球面模式新增常规输入首候选必须为球面的回归断言；`RunSpecializedDataset500.py --sphere-polynomial-only` 选择完整的 50 个非有理球面，保留退化预期，常规分母为 40。没有修改测试输入、优化器、两套验证算法、一般路径或公共接口。

## 测试方法

按用户要求仅执行球面测试。使用基线 Git 快照叠加本轮三个代码文件，独立构建 Release 库与两个测试程序，排除共享工作区中的其他未提交开发内容。构建配置为 MSVC 19.42、SDK 1.7.11；未固定 CPU 亲和性，机器非独占。

- 独立配对测试：`BSpline2Analytic.PairedSphere`，50 个非有理球面，包含 40 个常规和 10 个退化输入。按原模式比较旧投影验证器和新配对验证器；每类型先预热一个常规输入。真实首候选参与验证，来源证书和原误差要求保留。
- 完整生产转换：逐例子进程、无预热，保留原 20 秒外部看门狗。该时限包含读盘及来源认证；算法时间取内部四项计时器。
- 50 个输入均完成来源认证。归档中保存输入、修改文件和可执行文件哈希，以及生产验证内部五项互斥计时。

## 结果

| 检查 | 结果 |
| --- | --- |
| Release 构建 | 两个球面测试所需目标构建成功 |
| 独立配对球面 CTest | 通过，退出码 0；40/40 常规输入的新旧验证均 Passed，首候选均为 Sphere |
| 完整生产转换：常规输入 | 40/40 Success，0 超时；每例 1 次初始验证，0 次优化计算／优化验证 |
| 完整生产转换：退化输入 | 10/10 NoCandidate，仍未满足原 Degenerate 预期 |
| 完整生产转换运行器 | 退出码 1，原因仅为上述既有退化状态不符 |

此前超时的 **309、333、341、343、355、363、365、367、381、389** 全部恢复成功。与同一低次数据集的基线 CSV 对照，只有这十例的生产状态从 ProcessTimeout 变为 Success，其余四十例状态不变。

既有退化问题仍为 **317、319、337、339、357、359、377、379、397、399**。未修改断言接受这些返回，不能将专项 CTest 通过解读为完整 50 例生产回归全部通过。

## 耗时

下表仅统计 40 个常规输入，单位为毫秒／例。

| 指标 | 均值 | 最大值 |
| --- | ---: | ---: |
| 生产分类计算 | 3.004 | 5.299 |
| 生产分类验证（旧投影验证器） | 268.656 | 509.258 |
| 生产优化计算／优化验证 | 0 / 0 | 0 / 0 |
| 生产四项合计 | 271.660 | 514.297 |
| 独立新配对验证 | 8.436 | 19.250 |

此前十个超时输入的新配对验证均值为 8.351 ms、最大为 15.946 ms。生产四项合计低于 10 ms 为 0/40；仅新配对验证低于 10 ms 为 28/40。新验证器尚未替换生产入口，因此本轮解决的是候选排序导致的异常耗时和首候选不支持问题，不能宣称整个算法达到 10 ms。

两个测试入口的预热、运行方式不同，不据此给出严格加速比。四阶段及五个验证模块的均值、最大值见[结构化汇总](evidence/sphere-candidate-priority-20261015/summary.json)。

## 证据与复现

- [完整生产 CSV](evidence/sphere-candidate-priority-20261015/specialized-sphere.csv)、[内部验证计时明细](evidence/sphere-candidate-priority-20261015/specialized-sphere-details.csv)
- [独立配对 CSV](evidence/sphere-candidate-priority-20261015/paired-sphere.csv)、[CTest 日志](evidence/sphere-candidate-priority-20261015/paired-ctest.log)
- [生产运行日志](evidence/sphere-candidate-priority-20261015/specialized-run.log)、[构建日志](evidence/sphere-candidate-priority-20261015/build.log)、[来源与哈希](evidence/sphere-candidate-priority-20261015/provenance.json)

从仓库根目录构建 Release 的 `BSpline2AnalyticPairedValidationTests` 与 `BSpline2AnalyticDataset500Tests` 后：

```text
ctest --test-dir <build> -C Release -R "^BSpline2Analytic\.PairedSphere$" --output-on-failure
python tests/RunSpecializedDataset500.py --exe <build>/Release/BSpline2AnalyticDataset500Tests.exe --data tests/data/dataset500 --sdk <SDK根目录> --output <新的结果路径.csv> --timeout 20 --sphere-polynomial-only
```
