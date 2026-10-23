# 圆锥近锥顶端点因子：实现与定向验证

2026-10-23，Windows x64、MSVC 19.42、Release，数值文件 `/fp:strict`。基线为 `195b8fd1e7fbe2a36f1e468afb2a7c26a079bbd3`。本轮只改变母面距离验证器的位置后备；完整验收仍检查原边界位置和严格的 0.1° U/V 联合切向预算，位置容差仍为 1e-5。

## 实现

`MotherDistanceBounds.cpp` 优先使用原 Bernstein 距离界；原界通过时不构造因子。仅圆锥失败原单元保留归一化局部坐标，在 `MotherConeFactor.cpp` 尝试非负线性因子。源面提取、分类、优化、一般兜底和输入文件不变。

以端点平均高度比选择一个方向，只有小端/大端不超过 0.25 才尝试。该阈值仅控制工作量，不参与接受证明。每个原单元至多尝试两个因子：高度比提出的线性因子，以及必要时的零端点因子。归一化大端为 1，逆向递推商多项式，避免除以小端值。

令名义局部坐标为 P，源面提取及坐标变换误差为 eta_src。双精度递推只生成 Q，随后用向外舍入区间算术重构全部控制系数的 R=P-gQ；不假设小系数严格为零。线性乘积使用邻接系数公式，复杂度与控制点数成正比。因子 g 非负，圆锥母面 K 含锥顶且正缩放不变，因此

```
distance(X,K) <= eta_src + maxNorm(R) + max(g) * distanceBound(Q,K).
```

对 Q 重新构造二次隐式残差及分母，使用缓存权重的 `paired::Products`。Q 的高度非负时，可以用正分母商界及平方根残差界；到锥顶的范数上界是通用补充。负锥瓣不会借隐式残差为零通过。最终上界向上舍入恢复世界尺度，并严格与原容差比较。扰动不除以 g，也不修改源几何。

选定分解的隐式残差、径向平方、高度、重构余项和因子随节点一起细分，不在子节点重复分解或构造平方。因子保持原参数关系，例如右半子片的 g=t 变成 (1+tau)/2。分解状态属于当前候选，不写入源缓存。节点数、深度和截止时间沿用原预算。

新增报告字段：`coneFactorAttempts`、`coneFactorCells`、`coneFactorRemainderSkips`、`coneFactorAcceptedNodes`、`coneFactorBuildMs`、`coneFactorPerturbation`。构建计时包含在 `residualBuildMs`，不能重复累加；扰动是已构建分解在访问节点上的最大世界坐标上界，包含源包络。零次尝试表示未触发。内部位置接口的禁用开关仅用于同一验证器的测试对照，公开选项没有放宽验收的开关。

整体结果仍沿用已有的带浮点保护母面距离合同：新后备采用区间计算，不代表原快速分支升级为全流程严格区间证书，也不提供有限面片双向 Hausdorff 证明。`CheckPatch` 对真正包含锥顶的参数图仍返回 Unsupported；位置模块的锥顶白盒通过不等于完整转换支持退化边。

## 正确性验证

- `MotherConeFactor`：28,600 项检查通过。包含 3/5 次、严格零端点/1e-12/1e-4 正端点、U/V 反向和交换、空间旋转、1000 倍尺度及百万量级平移的 96 个组合，并独立采样检查上界。
- 容差内端点扰动被补回；超差扰动和负锥瓣不通过；候选变更不能复用旧证书。覆盖源包络占满预算、节点预算、截止时间取消，以及需要细分的非恒定误差。
- 非零端点的因子细分保留原参数；在最大深度内不能证明时仍未认证。原方法在近锥顶探针上达到深度 12 后失败，新方法在根节点通过。
- 圆弧截面的近锥顶面片通过完整 `ValidateMotherDistanceConversion`，确实触发因子；修改有限域后仍失败，未绕过边界检查。锥顶参数图仍 Unsupported。
- 原 `MotherDistanceUnit`：1,585 项检查通过。未运行退休的 `ValidateConversion` / `ValidateGeneralHybrid`，也未运行已废弃一般数据。
- 最终 CTest 的 `MotherConeFactor`、`MotherDistanceUnit`、`SequentialShortCircuit`、`SequentialFallback` 共 4/4 通过，记录见 [ctest.log](results/cone-factor-20261023/ctest.log)。

独立合成面使用多项式勾股参数化，完整验证器测试另用五次正弦/四次余弦截面。大坐标采样审计在实际存储控制点中心化后独立执行 de Casteljau，保留求值舍入余量；直接在世界坐标调用 SDK 求值的舍入量可能大于本例几何上界。采样只用于回归核查，不替代连续界推导。

## 耗时对照

编译修改前程序后保存其可执行文件，再编译新版本。55 个固定输入包括 Dataset500 全部 40 个常规非有理圆锥及既有选例中的其他 15 例。按案例交替先运行的版本，每版每例 5 次首次验证和 5 次源缓存复验，合计 **1,100/1,100 次验证和独立审计通过**。来源证书、SDK 初始化、读取、采样审计及输出在计时之外；首次时间包含分类与候选验证，未运行优化器或完整串联转换。原始数据及哈希见 [results](results/cone-factor-20261023/manifest.json)。

首次“分类＋验证”均值，单位 ms：

| 类型 | 案例数 | 基线 | 新版 | 变化 |
| --- | ---: | ---: | ---: | ---: |
| 平面 | 1 | 1.7405 | 1.6656 | -4.30% |
| 圆柱 | 1 | 2.2200 | 1.8164 | -18.18% |
| 圆锥 | 40 | 4.7587 | 4.7123 | -0.98% |
| 球面 | 6 | 6.3642 | 6.0863 | -4.37% |
| 圆环 | 7 | 8.8666 | 9.2203 | +3.99% |

普通圆锥整体验证未观察到均值回退，但残差构造子项从 0.3052 增至 0.3274 ms（+7.27%）。其他类型也有明显计时变化；本批测量不能证明各类输入固定加速或普遍满足 1%–2% 回退目标。没有为了得到更好数字重复正式批次。

每版 275 次首次测量中，基线 17 次、新版 19 次不低于 10 ms，最大值分别为 13.5673、15.3441 ms；没有声称全部达到 10 ms。此批未测完整串联转换，不能据此作其 50 ms 全量性能结论。

以下为单个合成面、**已准备源缓存的位置模块**，每种状态 100 次、逐次交替开关，共 800 次。两端都使用当前母面实现，禁用后备等价于保留原位置计算路径；不是完整转换耗时，也不是不同完整验证合同之间的比较。

| 探针 | 原位置路径 ms | 启用因子 ms | 结果 | 节点数变化 |
| --- | ---: | ---: | --- | --- |
| 普通截锥 | 0.01596 | 0.01645 | 两者 100/100 通过 | 1 → 1 |
| 近锥顶 | 0.04539 | 0.20446 | 原 0/100，新 100/100 通过 | 13 → 1 |
| 近锥顶＋2e-6 端点扰动 | 0.04576 | 0.18523 | 原 0/100，新 100/100 通过 | 13 → 1 |
| 近锥顶＋3e-5 超差扰动 | 0.04648 | 0.08348 | 两者均未认证 | 13 → 13 |

新分支降低了近锥顶探针的节点数，却比原来的快速失败多花约 0.14–0.16 ms；主要成本为区间重构和二次残差构造。超差位置探针因两次余项预算检查也增加成本。公开验证器有先行采样拒绝，这个白盒失败耗时不能直接当作完整负例耗时。当前收益是恢复认证，尚无证据宣称它能减少真实业务数据的候选重试或完整转换总时间。

## 复现

```powershell
cmake -S . -B build-cone-factor -G 'Visual Studio 17 2022' -A x64
cmake --build build-cone-factor --config Release --target BSpline2AnalyticMotherConeFactorTests BSpline2AnalyticMotherDistanceTests BSpline2AnalyticMotherDistanceSelectedTests BSpline2AnalyticSequentialTests -- /m:1 /nr:false
ctest --test-dir build-cone-factor -C Release -R '^BSpline2Analytic\.(MotherConeFactor|MotherDistanceUnit|SequentialShortCircuit|SequentialFallback)$' --output-on-failure
```

`BenchmarkMotherConeFactor.py --help` 列出对照参数；必须提供先前保存的基线程序、提交和 SHA-256，输出目录必须为新目录。完整命令记录在 [manifest.json](results/cone-factor-20261023/manifest.json)。[measurements.csv](results/cone-factor-20261023/measurements.csv) 保留全部首次/复验数据，[summary.csv](results/cone-factor-20261023/summary.csv) 为分类汇总；[position-probes.csv](results/cone-factor-20261023/position-probes.csv) 和 [position-summary.csv](results/cone-factor-20261023/position-summary.csv) 为位置白盒探针。[runs.log](results/cone-factor-20261023/runs.log) 包含命令、退出码与单测计数。SDK、许可证及构建产物不入库。
