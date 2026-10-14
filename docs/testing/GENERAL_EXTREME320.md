# 一般兜底极端输入：GeneralExtreme320

本批独立保存于 `tests/data/general-extreme320/`，只调用 `ConvertGeneral`。原有一般 600 例、额外 60 例、旧 500 例和其他极端数据均不执行。文件为 `extreme-000.b2a` 至 `extreme-319.b2a`，复用 `B2A_GENERAL_V1` 格式及只读加载器；来源真值仅供测试认证，不传入转换器。

本轮完整结果见 [2026-10-14 测试报告](GENERAL_EXTREME320_REPORT_2026-10-14.md)及[逐例 CSV](GENERAL_EXTREME320_RESULTS_2026-10-14.csv)。

退化 UV 三角域补丁后的专项结果见 [球面三角片 48 例报告](GENERAL_SPHERE_TRIANGLE48_REPORT_2026-10-14.md)及[逐例 CSV](GENERAL_SPHERE_TRIANGLE48_RESULTS_2026-10-14.csv)，该轮仅重测 ID 128–175。

## 矩阵与尺寸

| ID（含端点） | 组 | 输入数 | 配置 |
| --- | --- | ---: | --- |
| 000–063 | `radius_arc` | 64 | 圆柱、球、圆锥、圆环；名义半径尺度 100/1000；短边约 0.01 或弓高/容差为 0.1、1、10。 |
| 064–095 | `cylinder_helix` | 32 | 半径 0.1/1/10/100；方位角跨度 0.3/1.2/3/6 rad；左右旋及跨接缝。 |
| 096–127 | `slender_cylinder_helix` | 32 | 中心线长度 10/100/500/1000，宽度不低于 0.01，最大长宽比 1e5。 |
| 128–175 | `sphere_triangle` | 48 | 普通 UV 三角片、极点三角片、细长三角片各 16 个；南北极、接缝和参数方向。 |
| 176–207 | `cone_helix` | 32 | 半锥角 0.03/0.08/0.25/0.6 rad，起始高度 1/5/20/50；螺旋带及跨接缝。 |
| 208–239 | `slender_cone` | 32 | 半锥角 0.02/0.05/0.15/0.4 rad，沿母线或螺旋方向的细长带；横向名义宽度 0.0101。 |
| 240–271 | `torus_extreme` | 32 | 大小半径比最高 1e5，内侧狭颈、双接缝、小角度法线集中。 |
| 272–319 | `world_coordinates` | 48 | 6 个代表面片，各平移 0/1e3/1e6/1e9，方向为 `(d,-0.73d,0.41d)`。 |

偶数 ID 为有理输入，奇数为同一解析参数配置的非有理输入，共 160 对。大坐标组包含平面、普通/细长圆柱螺旋带、球面三角片、圆锥螺旋带及圆环。其中球三角片的平移副本也保留退化边；整批共 56 个带点退化边的输入。大坐标配对属于相同局部几何的数值表示测试，不宣称每个输入都具有独立形状。

半径、整体长度和宽度限定在 `[1e-2,1e3]`；圆锥有限片的轴向高度与截面半径也检查此范围（允许规定的点退化端）。角度、归一化源参数 `[0,1]^2`、曲率弓高和世界坐标不套用长度单位限制。三角片允许在指定边界收缩到零宽度，内部必须保持二维。

大半径组采用稳定关系 `h=2R*sin(theta/4)^2`。短边子组按最小截面半径反算跨度，使最短边目标为 `0.010001`，在业务下界之上留出 `1e-6` 余量；球的纬圆和圆锥上下截面不能直接使用同一半径推断边长。其余子组按弓高与 `1e-5` 容差的比例控制难度。

圆柱螺旋带使用 `u=u0+theta*s, z=z0+H*s+delta*t`。圆柱带的真实展开宽度为 `delta*R*abs(theta)/hypot(R*theta,H)`；构造器先指定物理宽度，再求 delta。圆锥采用相同的线性角度/高度映射，并使用锥面度量设置横向宽度。`length`、`width` 是构造特征量；圆锥螺旋长度一列记录起始处局部长度尺度，不冒充整条螺旋线的积分长度。四条边另外按解析目标用 1024/2048 段折线积分测长，记录 `edge_min/max/measure_change`。测长是数值检查，不是严格弧长区间证书。

所有四边在展开 UV 中为线段或点。普通三角片使用矩形域的一条边压缩到 UV 顶点的双线性映射；极点片保留不同 UV 角点在同一三维极点相遇的情况。周期跨度最多一周，不包含任意大圆三角形、多周螺旋或自交圆环。

## 样条构造及独立证书

1. 在局部世界坐标中计算解析目标，先减去一个局部参考点，再进行 Chebyshev–Lobatto 节点上的张量 Bernstein 插值。两个方向分别选次数；长角跨度按不超过 1.5 rad 分片，内部连接为 C0。数据符合一般模块连续输入约束，不声称分片具有精确 C1 连续性。
2. 非有理输入直接拟合坐标。有理输入拟合 `W*F`，其中 `W=(1+0.2s)*(1-0.15t)>0`，再按齐次控制点恢复。它是有误差证书的有理逼近，不声称螺旋边界精确可表示为有限阶 NURBS。
3. 在 50–200 个控制点预算内尝试 5/7/9/11/13 次，线性方向使用 1 次（有理分子使用 2 次）；只根据来源误差证书选择构造精度。插入非均匀节点并固定退化边控制点，最后添加世界坐标平移。本批实际为 50–148 个控制点，最高使用 9 次。
4. 将最终控制网写入文件，重新读取后认证。常规来源上界要求 `<1e-6`；小弧组额外要求 `<0.1*h`，避免构造误差掩盖曲率。生成器不调用转换器、不依据转换结果选择数据。

`GeneralExtremeCertificate.cpp` 根据测试端解析公式和 UV 映射，使用向外舍入区间与 Bernstein 系数界证明全域同参数距离，包含四条边。它复用区间基础运算及 Bézier 提取，不复用生产分类、投影或验收决策。

为了处理大坐标，认证先将实际存储控制点减去一个参考点，同时显式计入这次减法的区间误差上界。参考母面原点的平移也用区间减法；不会用平移前的理想控制网代替落盘数据。正权重保证控制点平移误差的凸包界可以应用于整个曲面。

四边形 UV 映射要求角点 Jacobian 保持同号；普通三角片要求恰有一条参数边压缩、三个不同 UV 顶点不共线。其 Jacobian 是非零三角形行列式乘 `s`、`1-s`、`t` 或 `1-t`，因此只在指定边界为零，内部不折叠。认证预算耗尽不能视为通过。

`--checks` 使用本批普通三角片、极点片及原点/大坐标平面做定向检查，要求正确目标通过，点状 UV 图或偏移 0.001 的平面被拒绝。没有调用其他数据集。

## 执行、超时与结果口径

正式运行每例重新计算来源证书，然后只向 `ConvertGeneral` 传入样条和默认选项：每类 8 次初值尝试、16 个粗评点、距离容差 `1e-5`。母面真值、类型和 PCurve 均不注入分类或优化。

成功结果同时要求生产连续验收通过、独立的有限解析面片同参数上界 `<1e-5`、SDK PCurve 为一次曲线且共享端点完全一致，并检查 SDK 与规范参数边界的三维求值一致。输出可以是其他类型，只要它满足要求；小面片不强制恢复来源类型。独立同参数界是保守充分条件，未证明时记录 `output_audit_unresolved`，不直接称作几何超差。

运行器逐例启动 SDK 子进程，默认 30 秒外部上限包含加载、输入认证、转换和输出认证。通过 `BEGIN/END` 阶段标记区分 `input_timeout`、`convert_timeout`、`audit_timeout`。转换一旦完成即保存阶段计时，输出认证超时仍保留已知的转换结果；转换尚未结束时，阶段时间写为 NaN，禁止伪造为零。

结果区分 `same_type_success`、`alternate_type_success`、`conversion_failed`、`input_unqualified`、`output_audit_unresolved` 和阶段超时。任何未通过案例都导致运行器返回 1。保留全部 320 行和输入 SHA-256；不删除、替换失败输入或扩大生产容差。存在未完成计时的组，只报告已完成案例的最大分类加优化时间，不能声称得到整组最大值。

## 独立构建及复现

使用 Python 3.9 或更高版本，从仓库根目录执行以下命令，仅构建新测试目标：

```powershell
cmake -S tests/general-extreme320 -B build-general-extreme320 -G "Visual Studio 17 2022" -A x64
cmake --build build-general-extreme320 --config Release --target BSpline2AnalyticGeneralExtremeTests --parallel 2
python tests/SummarizeGeneralExtreme.py --data tests/data/general-extreme320
ctest --test-dir build-general-extreme320 -C Release -R '^BSpline2Analytic\.GeneralExtreme320$' --output-on-failure
python tests/SummarizeGeneralExtreme.py --data tests/data/general-extreme320 --results build-general-extreme320/results.csv --json build-general-extreme320/summary.json
```

若宿主向 Windows 构建器传递了大小写重复的 PATH/Path 环境项，可通过 Python `subprocess.run(..., env=dict(os.environ))` 启动 CMake。本轮采用此方式。

运行器可直接配置范围、超时和输出文件：

```powershell
python tests/RunGeneralExtreme.py --exe build-general-extreme320/Release/BSpline2AnalyticGeneralExtremeTests.exe --sdk ../../../SGGK --data tests/data/general-extreme320 --output build-general-extreme320/results.csv --timeout 30
```

`--first/--last` 为半开编号范围；`--resume` 要求同一构建与选项，并检查已有结果的输入哈希。重新运行不要覆盖需要保留的历史结果。

生成器入口为 `--generate <新目录> [first last]`，拒绝覆盖已有输入；需要 SDK DLL 路径，并从 SDK 工作目录运行。全部文件就绪后写出独立 `manifest.csv`。固定输入已经入库，常规测试不重新生成。额外的 1e10/1e12 坐标精度探测未纳入本批 320 个正例。
