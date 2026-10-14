# 500 个固定面片

数据保存在 `tests/data/dataset500`，含 500 个 `.b2a` 文本文件和 `manifest.csv`。不是运行时随机选择容易通过的输入：固定种子生成后保存，测试加载同一批文件，重新验证数据条件，再按实际有理性直接调用 `Recognize` 或 `ConvertApproximate`。全部几何单位跟随 SDK 模型，不做单位换算。

## 配额

| 类型 | 精确普通 | 精确退化 | 非有理普通 | 非有理退化 | 合计 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 平面 | 50 | 0 | 50 | 0 | 100 |
| 圆柱 | 50 | 0 | 50 | 0 | 100 |
| 圆锥 | 40 | 10 | 40 | 10 | 100 |
| 球面 | 40 | 10 | 40 | 10 | 100 |
| 圆环 | 50 | 0 | 50 | 0 | 100 |
| 合计 | 230 | 20 | 230 | 20 | 500 |

控制点数以最终输入控制网的 `nu*nv` 统计，包含退化边上重合的顶点；权重不另计数量。50–79、80–109、110–139、140–169、170–200 五档各 100 个；每档、每种曲面内精确及非有理各 10 个。包含恰好 50 和 200 个控制点的面片。每个文件还保存一个较小的来源参考面，它不是额外测试输入，也不适用输入控制点配额。

250 个解析参考面各产生一对案例：偶数编号为精确转换，紧随的奇数编号为非有理拟合。`pair=id/2`；类型索引 0–4 分别为平面、圆柱、圆锥、球面、圆环。每类连续 100 个文件。

精确平面的权重通常全部为 1，因此 250/250 是来源配额，而非实际接口分支配额。实际预期分流为 200 个有理识别、300 个非有理转换。

## 构造方式

1. 用固定整数种子改变位置、轴向、尺寸、半锥角、圆环半径比和矩形参数域。普通圆弧跨度约 20.6–35 度，球极点片段约 28.6 度；不声称覆盖全部大跨度、短圆弧或浅锥情形。圆环采用普通 `R>r>0`。
2. SDK `ToBSpline` 从解析面生成等参矩形 NURBS。将其样条参数域仿射规范化到 `[0,1]^2`，不把样条参数当作解析角度。
3. 精确路径使用升阶、均匀/非均匀节点插入及节点重数变化达到控制点档；节点位置与次数不用于识别标签。
4. 非有理路径在同一参数域以实际来源曲面点做张量 Bernstein 插值，使用带列选主元的 QR，尝试满足控制点预算的次数。平面用三次；其他类型尝试五、七、九次，圆柱/圆锥轴向用三次。不是删除 NURBS 权重或直接把采样点当控制点。
5. 插入节点补足控制点档；对两种输入及参考面同时执行 U/V 反向、交换。元数据 `variant` 的低三位分别对应 U 反向、V 反向、UV 交换。
6. 退化案例是锥顶或球极点。将解析极限位置的浮点残差归并到一个点，并对最终输入的整排退化边控制点强制相等。拟合及节点处理之后重新验证，不以“短边”替代点边。

当前参考几何在既定预算内均成功生成，没有根据转换器成功与否更换案例。生成器对预算内无法通过配对误差证明的情况报错；不自动扩大控制点上限或放宽容差。

## 数据入库证明与转换验收

`PairedSurfaceBounds` 是测试专用的独立同参数曲面对比器，不调用识别器或使用其分类结果证明输入合格。对两张曲面 `X1/W1` 和 `X2/W2`，直接由控制网进行向外舍入的 Bézier 提取，在相交节点矩形上验证：

```text
E = X1*W2 - X2*W1
distance <= upper(norm(E)) / lower(W1*W2)
```

要求正权重，必要时细分；工作量或深度不足时不能通过。拟合选择阶段目标上界 `<8e-7`；最终输出及文件重载后必须通过 `<1e-6`。精确转换案例对保存来源的数值偏差必须通过 `<1e-9`。不同停止阈值会导致返回上界不同，它们不是实际最大误差的测量值。

参考面和最终面使用同一完整参数域，因此同参数距离上界足以约束两张有限面片的几何距离，包括对应边界；不需要假定投影参数线性。

该入库检查不代替转换验收。转换时仅给 `Recognize` / `ConvertApproximate` 输入面，不提供已知类型、真实参数或参考面。普通非有理分支仍要求有限面片距离 `<1e-5`，边界切向和法向分别 `<0.1` 度。精确分支保留原 `Recognize` 的验证口径。退化案例目前期望 `Degenerate`，不表示已支持锥顶或球极点的成功转换；退化边角度不填零、不计入普通误差统计。

本数据的拟合误差小于转换距离容差，很多初始候选可直接通过，优化时间可能为零。原平台扰动及五类扰动测试继续负责检验必须优化才能通过的情形。

## 运行

```powershell
# 从模块目录运行，设置 SDK 和许可证后：
.\Build.ps1 -Configuration Release -Dataset500Only

# 已构建时单独运行：
ctest --test-dir build -C Release -R '^BSpline2Analytic\.Dataset500$' --output-on-failure

# 可选：使用 Python 标准库汇总 CSV（输出到 stdout）
python tests/SummarizeDataset500.py build/dataset500-Release.csv

# 仅圆锥：加载 200–299，并运行圆锥停止行为、候选回退及扰动检查
.\Build.ps1 -Configuration Release -ConeOnly
python tests/SummarizeDataset500.py --cone-only build/dataset500-cone-Release.csv

# 仅非有理球面：301、303、…、399，共 50 个（40 个普通、10 个退化）
.\Build.ps1 -Configuration Release -SpherePolynomialOnly
python tests/SummarizeDataset500.py --sphere-polynomial-only build/dataset500-sphere-polynomial-Release.csv
```

全部 CTest 默认也包含此测试。准备数据、读取文件、独立配对验证、每类型/实际分支首次预热及日志均在转换计时之外。每个输入正式转换一次，均值和 P95 是不同输入间的统计，不是同一案例多次重复测量；没有绑定 CPU 或排除系统负载影响。退化输入的提前返回时间单列，不能混入普通转换速度。

CSV 保留每个案例的 `classify_ms` 与 `optimize_ms`，汇总脚本在“类型/实际分支”和“控制点档/实际分支”两组表中同时输出 `max_i(classify_ms[i] + optimize_ms[i])` 及对应案例文件名。先对同一案例的两个阶段求和，再取组内最大值，不将两个阶段各自的最大值相加；不含验收或计时外围开销。并列时列出编号最小的案例。退化分支另行输出相同指标。完整的旧 CSV 可直接重新汇总，无需重新运行转换。

`-ConeOnly` 只运行 `ConeDataset` 和 `Cone` 两组：前者读取 100 个已有圆锥输入（80 个普通、20 个退化），后者独立检查案例 249 的旧初值强制优化、新初值零优化、失败候选回退及圆锥平台扰动。前者写入 `build/dataset500-cone-Release.csv`；后者预热一次后对 249 测量五次，写入 `build/cone249-Release.csv`。不覆盖原 500 例 CSV。也可直接使用 `BSpline2AnalyticDataset500Tests --cone-only <数据目录> <结果.csv>`。

`-SpherePolynomialOnly` 只运行 `SpherePolynomialDataset`，跳过涉及其他类型的 `Probe`，但保留该子集全部来源、控制点档、退化配额、真实分支及成功验收检查。独立输出 `build/dataset500-sphere-polynomial-Release.csv`，不覆盖全量结果。命令行也可使用 `BSpline2AnalyticDataset500Tests --sphere-polynomial-only <数据目录> <结果.csv>`。汇总脚本在此模式下要求编号完整且类型、来源、实际分支全部匹配；不将缺行 CSV 视为完整结果。

验证内部现在还输出准备、曲面距离界、曲面判定/细分、边界验证、报告整理五项互斥计时，分别累计初始候选和优化后调用。专项脚本 `SummarizeSphereValidationTiming.py` 接受一轮或多轮该球面 CSV，检查五项之和和基线结果一致性。五步详细定义及本轮 50 个非有理球面的三轮结果见[验证计时报告](SPHERE_VALIDATION_TIMING_2026-10-14.md)。

历史完整重测见 [2026-10-13 报告](DATASET500_RETEST_2026-10-13.md)，独立 CSV 为 `build/dataset500-20261013-current-Release.csv`。

新增的 CSV 尾部字段记录分类缓存构建次数、缓存中的曲面点及偏导求值次数、原始候选数、参数域恢复次数和优化迭代数。有理路径这些近似分类计数为零；`class_cached_points` 只统计分类缓存的 `CalcPoint` 调用，不含等参曲线拟合、优化采样和验收求值。旧 CSV 无需补列即可继续汇总。

手动运行程序需在 SDK 根目录运行，并将 SDK 对应配置 DLL 目录加入 PATH。生成器命令为 `BSpline2AnalyticDataset500Tests.exe --generate <输出目录>`。建议另选空目录，不覆盖基准输入；重新生成可能受 SDK/编译器浮点差异影响，必须重新验证。

只读诊断命令为 `BSpline2AnalyticDataset500Tests.exe --inspect <case-XXX.b2a>`：对已知参考母面和算法生成候选分别给出采样诊断及至多 500 分片的验证结果，帮助定位问题。该模式使用已知参考面仅作对照，不计为正式转换结果，不能用其 `Unresolved` 推断默认预算下也失败。

## B2A_CASE_V1 文本格式

文件用空白分隔数值，double 以 17 位有效数字保存：

1. 标记 `B2A_CASE_V1`。
2. `id pair type bucket variant seed refitted collapsedMask storedBound`。
3. 解析面原点、轴向、X 方向各三个坐标；然后是 `radius majorRadius halfAngle`，仅使用对应类型的有效参数。
4. 原始解析参数域 `u0 u1 v0 v1`，方向变换前的 SDK 域。
5. 来源样条和输入样条各一段：`degreeU degreeV nu nv hasWeights`；U 节点数及每个 `knot multiplicity`；V 同理；随后 U 外层、V 内层的控制点 `x y z [weight]`。

`collapsedMask` 低四位依次对应 Umin、Umax、Vmin、Vmax。`storedBound` 是生成时记录；测试不会只信任它，而会重新计算。`manifest.csv` 是便于浏览的索引，不是独立证书。

## 本批发现并修复的问题

- 圆锥有限域验收：角度界原先用两个独立坐标范围相除，丢失共同高度缩放相关性。现用正分母 Bernstein 系数比值的凸包界，并保留非正分母时不可判定的返回。阈值、预算没有放宽。
- 球面候选排序：近似圆心的小偏差可能生成大半径接近零的圆环，它偶然更低的采样评分会导致长时间尝试奇异极限。当已有采样距离低于最终容差的球面候选时，将 `abs(majorRadius)<=distanceTolerance` 的圆环排到后面；所有候选仍保留，接受与否仍取决于完整验收。
- 圆环有限域验收：管截面角度原先先取独立的 `(rho-R,z)` 区间再旋转，丢失圆弧相关性。令 `L=z*cos(ref)+R*sin(ref)`，在分母不跨零时用 `(L²-rho²*sin²(ref))/(L+rho*sin(ref))` 计算旋转后的有向横坐标界，分子保持 Bernstein 多项式；分母可能跨零时回退到保守直接区间式。两端都靠近时选择较近的边界作为参考角，让稳定表达同时适用于上、下边界。数据、误差阈值和预算均未改变。

最终计时及实际通过情况见 [TEST_REPORT.md](TEST_REPORT.md)。

## 专用路径四项计时（2026-10-14）

500 例直接调用等参专用接口，不执行一般兜底或统一入口附加 PCurve 验证。有理 200 例、非有理 300 例；其中 40 个坍缩边案例应返回 `Degenerate`，不计为转换成功。固定数据不重新生成。

`timing_schema=specialized-v2`：`classify_ms = classify_compute_ms + classify_validate_ms`，`optimize_ms = optimize_compute_ms + optimize_validate_ms`，评价耗时为两者之和 `classify_optimize_ms`。`validate_ms` 在这两个独立接口中为零。初始候选的全部验证尝试归入分类验证；每次求解后的验证归入优化验证，包括失败、零迭代及重试。每项额外输出 calls；它不是优化迭代次数。`jacobian_ms` 是优化计算内部子项，不重复加入总计。

来源重建、输入证书、SDK 初始化、每类型/实际接口首次非退化案例的一次预热和日志均不计入转换时间。每个固定输入正式转换一次，P95 为跨案例 nearest-rank。四项结果汇总：`python tests/SummarizeSpecializedDataset500.py <结果.csv> <汇总目录>`；脚本核验 500 个 ID、分支配额、四项与总计关系，并分别统计普通和退化案例。

本口径的完整复测结果见 [2026-10-14 四项计时报告](SPECIALIZED500_TIMING_REPORT_2026-10-14.md)。
