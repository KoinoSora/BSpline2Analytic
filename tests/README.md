# 测试代码与固定输入

## 当前测试材料导航

完整的构造方式、设计思路、C++ 接口及单例／多例命令见 [测试材料说明](../docs/testing/materials/README.md)。材料提供全部固定输入与测试源文件索引，区分当前串联、模块诊断及已停用的历史入口。下文历史段落不能覆盖 `AGENTS.md` 的禁测约束。

## 等参优先、失败兜底的统一测试

本轮结果见 [1,160 例串联测试报告](../docs/testing/SEQUENTIAL_NON_EXTREME_2026-10-16.md)：3,480 次全部实测不超过 50 ms，947 例三次均小于 10 ms；985 个普通正例通过 920 个，183 例预算失败，回归未全部通过。

`SequentialConversion.h::ConvertIsoThenFallback` 对每个输入先运行等参分支（有理输入为精确识别，非有理输入为增量等参分类＋母面验证），失败后运行一般分类、当前有理边界验证，必要时优化母面参数并重新验证。所有候选和两分支共享 50 ms 总预算，等参子预算默认为 20 ms；触及子预算可进入兜底，总预算退出或实测超过 50 ms 均返回失败。两个分支保留各自证书合同，不能将兜底的内部采样解释为连续母面距离证明。

`RunSequentialNonExtreme.py` 串行运行专用 500 例、新版兜底 main 600 例和 extra 60 例，每例三次新缓存测试，不做逐例预热。极端和压力集不运行。每例 `total_ms` 为 API 调用前后的真实墙钟，包含分类、全部失败候选、优化、验收和返回开销，不含输入读取、来源证书、独立输出审计、CSV 或控制台输出。单例是否小于 10 ms 取三次中的最大值。

所有输入使用统一默认配置，不向算法传入母面真值、参数图或测试标签。extra 的已知母面拒绝、非线性边界、内部采样探针、非法选项及节点预算标签，以及专用数据的坍缩边，单列为诊断数据；非法选项/节点预算文件在本次是几何输入覆盖，不冒充原先的定母面/特殊配置测试。正例必须三次成功且输出审计通过；无解负例必须三次非超时拒绝，超时不能算正确拒绝。复现：

```powershell
cmake -S . -B build-sequential-20261016 -A x64
cmake --build build-sequential-20261016 --config Release --target BSpline2AnalyticSequentialTests -- /m:1 /nr:false
python tests/RunSequentialNonExtreme.py --exe build-sequential-20261016/Release/BSpline2AnalyticSequentialTests.exe --sdk ../../../SGGK --data tests/data --output build-sequential-20261016/full --runs 3
```

`cases.csv` 提供逐例正确性、三次最大耗时和 10/50 ms 判断；`runs.csv` 保留原始测量；`failures.csv`、`over-10ms.csv`、`timeouts.csv`、`over-50ms.csv` 分别列出未通过、10 ms 未达标、预算退出和真实超时。脚本返回 1 表示存在未通过/预算失败，返回 2 表示执行或记录完整性错误。

当前一般兜底数据为 [`data/general-rational-v2/`](data/general-rational-v2/)，共 985 例（600 主集、60 专项、320 极端、5 压力）。正例从母面 **NURBS 参数域直线边界**构造；20 个非线性边界探针继续作为单独的合同探针。旧 `dataset500/general-*.b2a`、`general-extra-*.b2a` 和 `general-extreme320/` 已废弃，仅保留历史证据。专用路径 `case-*.b2a` 不受影响。

985 例两方向次数均为 **3–5**，控制顶点总数 **≤200**，算法几何容差保持 **1e-5**；有意超差的负例/探针不改变验收阈值。生成器、清单发布和运行器均检查表示约束。

后续三角形参数域、微小分段和圆锥对应修正见 [99 例定向回归](../docs/testing/GENERAL_RATIONAL_REPAIRS_2026-10-16.md)：73 个目标失败恢复 66 个，26 个控制案例结果不变，剩余 7 个失败单列。默认允许源/目标三维直径界均不超过 `1e-8` 的小段豁免切向，但仍连续认证 `1e-5` 位置约束，并在输出 `tangentExemptions` 记录区间；角度上界不覆盖这些豁免段。设置 `microSegmentDiameter=0` 可禁用该策略。CSV 新增 `micro_merged_segments`、`tangent_exemption_count`、`recovery_3d_refinements`，逐例另存 `*-exemptions.csv`。

后续已获准并完成新版全量回归：**985 例 × 3 次**，默认单例共享 **50ms** 分类＋验证预算，实际最大 **18.7393ms**，无时间预算退出；850 个正例中 710 个认证、140 个未认证。见 [全量耗时与失败报告](../docs/testing/GENERAL_RATIONAL_985_50MS_2026-10-16.md)。单测及残差包络精确有理数审计也已执行通过。以下“未运行”记录描述此前生成/残差替换阶段，不再表示最新测试状态。

`RunGeneralRationalBoundaryFull.py --max-case-ms 50` 默认保留 5ms 清理余量、45ms 协作取消。CSV 记录真实耗时、`time_budget_exceeded`、`timeout_stage` 和 `within_time_limit`；汇总将时间取消与几何未认证分开。`over-time-limit.csv` 列出实际超过上限的案例，`time-budget-cases.csv` 列出预算中止案例。运行器对进程错误或实际超时返回非零；正例未认证须查看 `failures.csv`，不能以进程运行成功代替回归通过。`--max-case-ms 0` 只用于明确的不限时诊断；默认运行绝不取消 50ms 约束。

2026-10-16 本轮仅生成并独立认证数据、编译迁移后的入口，**没有运行兜底单测、转换、分类或性能测试**。来源认证通过不代表兜底识别通过，也不产生新的 10 ms 性能结论。构造方法、次数变化、来源证书和命令见 [V2 数据说明](../docs/testing/GENERAL_RATIONAL_DATASET_V2.md)。

后续已将一般兜底的切向叉乘后备替换为正比例缩放残差，见 [实现与编译记录](../docs/testing/GENERAL_TANGENT_RESIDUAL_2026-10-16.md)。本轮仅编译和静态检查，未运行单测、包络导出/数值审计或数据集。全量及选例 CSV 新增 `tangent_residual_checks`、`tangent_residual_accepted`、`tangent_residual_interval_fallbacks`；后续获准运行时用于区分残差构造次数、通过次数及区间残差调用次数（含直接调用）。`AuditGeneralBoundaryEnvelope.py` 已改为精确有理数残差系数和角度界审计，需要重新导出的残差格式；不读取旧点积/叉乘导出作为新实现证据。

只生成和认证数据：

```powershell
cmake -S . -B build-rational-v2 -A x64
cmake --build build-rational-v2 --config Release --target BSpline2AnalyticGeneralRationalDatasetTests
python tests/GenerateGeneralRationalDataset.py --exe build-rational-v2/Release/BSpline2AnalyticGeneralRationalDatasetTests.exe --sdk ../../../SGGK --legacy-main tests/data/dataset500 --legacy-extreme tests/data/general-extreme320 --output tests/data/general-rational-v2 --resume --audit
python tests/AuditGeneralRationalMigration.py --data tests/data/general-rational-v2 --legacy-main tests/data/dataset500 --legacy-extreme tests/data/general-extreme320
```

后续获准测试时，`RunGeneralRationalBoundaryFull.py`、`RunGeneralRationalBoundarySelected.py`、`CompareGeneralRationalBoundarySelected.py`、`BenchmarkGeneralBoundaryEnvelope.py` 的 `--data` 均传入新版根目录。运行器检查 V2 合同、完整清单和文件 SHA256；不再接收 `--extreme-data`。测试程序检查每个文件的 V2 标记，并重新认证实际输入。批量入口使用新的 `--rational-v2-only` 标记，必须重新构建，不能将旧可执行文件或旧 240 例失败清单混入新版批次。

全量运行器默认采用“分类＋当前验证器”，不向分类器提供来源母面或参数图；`--fixed-chart` 是单列的已知母面/参数图诊断，压力及指定专项也会明确标记入口。来源认证、输出审计和诊断均在算法计时外。`total_ms` 包括分类及验证；汇总从新版清单读取预期，不继承旧案例 254 的未认证结果。新恢复诊断要求新版基线及相同输入哈希，不再固定要求 240 例。

旧 `ValidateGeneralHybrid`、`ValidateConversion` 的禁测约束继续有效。`Build.ps1 -GeneralOnly` / `-GeneralPositiveOnly`、旧一般数据/极端数据/分类计时入口直接报废弃，旧 CTest 已禁用。历史 [985 例报告](../docs/testing/GENERAL_RATIONAL_FULL_2026-10-15.md) 和 [240 例诊断](../docs/testing/GENERAL_RATIONAL_MAIN240_DIAGNOSIS_2026-10-16.md) 仅供查阅，同编号不表示新旧输入相同。

| 文件/目录 | 职责 |
| --- | --- |
| `RecognitionTests.cpp`、`ApproximateConversionTests.cpp` | 有理识别与非有理转换回归。 |
| `PlateauTests.cpp`、`JacobianTests.cpp`、`ConeTests.cpp` | 平台扰动、导数及五类扰动、圆锥专项。 |
| `Dataset500Tests.cpp`、`Dataset500.*`、`PairedSurfaceBounds.*` | 固定数据加载、配额与来源证书校验、数据生成。 |
| `ReferenceCases.*`、`ApproximateCases.*` | 各测试共享的精确与非有理案例。 |
| `ApproximateBenchmark.cpp`、`SummarizeDataset500.py` | 性能记录与 CSV 汇总。 |
| `PairedValidationTests.cpp`、`SummarizePairedValidation.py` | 独立参数对应验证器专项、非有理球面及完整 Dataset500 分流对照；不替换正式转换入口。 |
| `PairedPreparationBenchmark.cpp`、`BenchmarkPairedPreparation.py` | 同一程序对比指定提交的独立验证器、当前首次验证与源面缓存命中；记录准备子项和首次缓存建立成本。 |
| `data/dataset500/` | `case-*.b2a` 当前为 `dataset500-degree5-v2`：两方向最高五次、控制点严格小于 200；一般路径数据另有独立清单。 |
| `AuditDataset500Migration.py` | 对比迁移前后的 500 个文件，检查身份、来源、解析参数和变更范围，并记录全部 SHA-256。 |
| `GeneralDataset.*`、`GeneralDatasetCertificate.cpp`、`GeneralDatasetTests.cpp` | 保留的旧一般数据实现；旧测试入口已废弃。 |
| `SummarizeGeneralDataset.py` | 仅用于查阅历史一般兜底 CSV。 |

原一般兜底数据、生成器和转换计时说明现已归档，不再按 [GENERAL600.md](../docs/testing/GENERAL600.md) 或历史分类报告中的命令重跑。原数据只允许迁移生成器读取规格。新版独立参考和数据工具见上文。

构建和复现命令从仓库根目录执行，见 [Build.ps1 用法](../README.md#构建与运行)及[数据集说明](../docs/testing/DATASET500.md)。测试通过 CMake 私有引用 `src/` 内部头文件，调用方只依赖 `include/`。
`ExtremeDataset.*` 和 `ExtremeTests.cpp` 提供独立的 200 例极端输入与 7 个模块诊断。`RunExtremeDataset.py` 隔离逐例执行并记录超时，`SummarizeExtremeDataset.py` 校验完整性、汇总分组耗时。`./Build.ps1 -ExtremeBoundaryOnly` 当前运行新的 `data/business-extreme200/` 业务尺寸矩阵，见 [BUSINESS_EXTREME200.md](../docs/testing/BUSINESS_EXTREME200.md)；历史 `data/extreme200/` 与[旧说明](../docs/testing/EXTREME200.md)保留。新数据生成用 `--generate-business`，旧数据生成用 `--generate`，新旧同编号不可混用。完整结果与每组最大分类加优化耗时见[业务尺寸极端测试报告](../docs/testing/BUSINESS_EXTREME200_REPORT_2026-10-14.md)。

`ClassificationGateTests.cpp` 检查短圆弧、数值共线、双候选保留、失败状态与通用回退，并读取原先返回 `InvalidInput` 的 8 个固定输入检查边界拟合。构建目标为 `BSpline2AnalyticClassificationGateTests`；定向运行：`ctest --test-dir build -C Release -R '^BSpline2Analytic\.ClassificationGate$' --output-on-failure`。

极端测试运行器支持 `--first 96 --last 128`，只运行角跨度门槛组（右端不含）。在 SDK 工作目录并配置其 DLL 路径后，执行 `python <repo>/tests/RunExtremeDataset.py --exe <repo>/build/Release/BSpline2AnalyticExtremeTests.exe --data <repo>/tests/data/extreme200 --output <result.csv> --first 96 --last 128`。`--resume` 仅接受同一选定范围的连续前缀。每个子进程额外保存 `case-NNN-fits.csv`，保留专用分支的逐等参线拟合诊断，包括回退前的拒绝原因；采样拟合误差不等于连续界。32 例的几何转换仍有未通过结果时，运行器返回 1，不能把定向单元测试通过等同于整组转换通过。结果见 [门槛修复报告](../docs/testing/CLASSIFICATION_GATE_REPORT_2026-10-14.md)。

## 等参专用路径计时

当前三套固定数据全量入口为 `BSpline2AnalyticMotherDistanceDatasetTests` 和
`RunIsoSpecializedFull.py`：Dataset500 500 例、业务尺寸极端集 200 例、历史极端集 200 例。
非有理输入接入 `IsoClassificationSearch` 的两子区间起步、按需加密及候选重试，再调用当前
`mother::ValidatePrepared`；有理输入调用 `Recognize`。不调用旧验证器、一般兜底或曲面优化器。
每例三次新缓存运行，不做逐例预热；所有失败候选验证累计进分类验证。来源认证和局部坐标输出审计在计时外。
超过五次的极端历史输入保留并报告不支持，不自动改写数据。结果与完整复现命令见
[900 例全量回归报告](../docs/testing/ISO_SPECIALIZED_FULL_2026-10-16.md)。
该脚本返回 1 表示存在未认证案例或诊断不符合预期，返回 2 表示执行/记录完整性错误；不能将回归已完成等同于全部通过。

第三个独立母面距离验证器使用目标 `BSpline2AnalyticMotherDistanceTests`。构建后执行 `ctest --test-dir <build> -C Release -R '^BSpline2Analytic.MotherDistance(Unit|SphereTorus)$' --output-on-failure --no-tests=error`，只运行代表性覆盖与 100 个非有理球面/圆环输入的三轮对照，输出 `mother-sphere-torus-Release.csv`。对照为当前 `ValidatePairedConversion`，不调用废弃验证器。`python tests/SummarizeMotherDistance.py <csv> <output>` 核验 600 条记录并生成五模块、逐案例和状态汇总。`algorithm_ms` 仅为分类计算加首候选验证，未运行优化器。合同差别、覆盖范围和结果见[专项报告](../docs/testing/MOTHER_DISTANCE_SPHERE_TORUS_2026-10-15.md)。

流式提取的定向性能测试使用 `BSpline2AnalyticMotherDistanceSelectedTests` 与 `BenchmarkMotherStreaming.py`，只选 16 个代表性及历史慢例，分别记录首次验证和源缓存复验。`AuditMotherExtraction.py` 用 80 位 Decimal 独立核对双精度 Bézier 控制点的误差包络；这些审计在计时外执行。改前、改后五模块耗时、测试命令及 10 ms 目标的剩余差距见[流式提取报告](../docs/testing/MOTHER_STREAMING_2026-10-16.md)。

无叉乘残差版本新增 `CompareMotherResidualSelected.py`，对保存的上一版联合预算可执行文件核验哈希，按案例交替执行修改前后版本，各记录五轮首次/缓存复验。新列 `tangent_residual_fits`、`tangent_residual_splits` 统计比例拟合和细分。最终单测 1585 项断言通过；16 例新版本 160 次验证与独立审计通过，五模块及逐例比较见[残差认证报告](../docs/testing/MOTHER_RESIDUAL_SELECTED_2026-10-16.md)。

2026-10-16 最新要求将母面验证器边界角度检查改为 U/V 参数切向联合预算，严格要求同一子区间的 `sin²(U夹角上界)+sin²(V夹角上界) < sin²(0.1°)`，不进行直接法向认证或失败回退，不检查相邻面 G1。新 CSV 分别输出 `tangent_u_bound_deg`、`tangent_v_bound_deg`、`joint_tangent_bound_deg`、`joint_sine_squared_bound`、`joint_tangent_failure`；联合界是各共同子区间联合界的最大值，不能用全边 U/V 各自最大值重新相加替代。第二验证器仍单独输出 `paired_tangent_bound`、`paired_normal_bound`，母面行留空；`boundary_angle_ms` 计入两个方向及联合预算认证。超限见证、无法认证和节点预算耗尽分别用 `JointTangentConditionExceeded`、`JointTangentConditionUnresolved`、`JointTangentBudgetExceeded` 说明，失败不等同于实际法向角超限。联合预算实现时仅编译；随后获准执行的 [16 例定向测试](../docs/testing/MOTHER_JOINT_SELECTED_2026-10-16.md) 已完成，160 次验证与独立审计通过。当前选例运行器为 `RunMotherJointSelected.py`，计时外直接用解析母面导数核对两个切向方向。历史 CSV 保留原合同和原数值，未重跑全 500 例。

2026-10-15 起，旧投影验证器 `ValidateConversion` 废弃，保留代码和历史数据但不再测试、预热或作对照。调用它的 CTest 项已禁用，`RunSpecializedDataset500.py` 已停止执行。正式 `ConvertApproximate` 尚未接入新验证器，其完整流程测试暂缓；以下历史报告中的旧命令不再作为执行入口。

当前使用 `./Build.ps1 -Dataset500Only`（500 例）或 `-SpherePolynomialOnly`（50 个非有理球面），输出分别为 `build/paired-dataset500-Release.csv` 和 `build/paired-sphere-Release.csv`。程序 `BSpline2AnalyticPairedValidationTests` 的 `--dataset` / `--sphere` 与 `--dataset-new-only` / `--sphere-new-only` 均仅使用新验证器。非有理输入验证实际分类首候选，不运行候选重试或优化；200 个有理输入另外调用精确识别 `Recognize`，不套用五模块计时。

CSV 的 `validation_ran` 标明是否实际进入新验证器。仅对此列为 1 的行汇总 `prepare_ms`、`surface_ms`、`subdivide_ms`、`boundary_ms`、`report_ms`，它们之和等于 `validation_total_ms`；边界内部细分包含在 `boundary_ms` 中。有理识别输出 `exact_status`、`exact_ms` 及分类计算/验证明细。未调用时旧兼容列中的零不是一次零耗时测量；`old_status` 恒为 `NotRun`。来源证书、预热、候选一致性审计与输出采样检查均不计入这五项。

使用 `python tests/SummarizePairedDataset500.py <500例CSV> <输出目录>` 核验完整性、排除未调用的零值并按曲面类型统计五模块。最新结果和口径见 [500 例报告](../docs/testing/PAIRED_DATASET500_NEW_ONLY_2026-10-15.md)。

位置误差界优化沿用比较工具：`python tests/BenchmarkPairedPreparation.py --baseline d345b30 --build-dir build-specialized-timing --output <comparison.csv> --runs 3`。新增位置内部计时 `compose_ms`、`product_ms`、`difference_ms`、`fast_bound_ms`、`tight_bound_ms`，均为 `surface_ms` 的子项；`fast_passes` 和 `tight_checks` 统计快速界成功及紧界调用的节点数。CSV 汇总兼容旧字段集合。实现范围、专项和同批比较见[位置误差界优化报告](../docs/testing/PAIRED_POSITION_2026-10-15.md)。

准备阶段优化的复现方式：先构建当前 Release 库，再执行 `python tests/BenchmarkPairedPreparation.py --baseline 3985a1b --build-dir build-specialized-timing --output <comparison.csv> --runs 3`。脚本提取旧实现到忽略的构建目录，只重命名私有命名空间和公开函数，并增加两个子计时器；两版共用未修改的区间基础运算和同一编译设置。40 个普通非有理球面的三种模式轮换顺序，各测三次；缓存首次填充另列，不能把缓存命中耗时当作首次验证耗时。见[准备阶段优化报告](../docs/testing/PAIRED_PREPARATION_2026-10-15.md)。

历史完整转换测试直接调用 `Recognize`（200 个有理输入）及 `ConvertApproximate`（300 个非有理输入）。现已暂停这一旧验证器流程；固定数据的 40 个退化案例保留退化预期。

低次数据集的迁移、证书与范围核验见[迁移报告](../docs/testing/DATASET500_LOW_DEGREE_2026-10-15.md)。数据检查可独立运行 `BSpline2AnalyticDataset500Tests --audit <dir>`；所有业务输入的 U/V 次数必须不超过 5，控制点必须小于 200。其他数据集和验证器高次运算单测保留各自合同。

历史 `RunSpecializedDataset500.py` 曾提供完整转换的逐例看门狗，现保留源码但入口直接退出。历史 CSV 中进程时限包括初始化、读盘和认证，不等同于算法时间；超时／异常的未知耗时不能作为零耗时或通过。

独立验证测试仍检查来源证书、误差和覆盖率；覆盖不足依旧返回失败。首候选验证不等同于带候选重试和优化的完整转换。

CSV 标记 `timing_schema=specialized-v2`，分别输出分类计算、初值验证、优化计算、优化后验证及其调用次数。分类总计和优化总计分别包含对应验证；`validate_ms` 为零，不能再把四项明细加到阶段总计上。用 `python tests/SummarizeSpecializedDataset500.py <csv> <目录>` 校验、生成 `summary.json` 与按类型/分支、控制点档分组的 `groups.csv`。旧报告使用旧计时定义，不能直接比较 `classify_ms` 和 `optimize_ms`。

验证内部五步计时标记为 `validation_timing_schema=exclusive-validation-v1`，`initial_*` 和 `refined_*` 分别累计初始候选及优化后验证。后缀为 `prepare`、`surface_bounds`、`subdivide`、`boundary`、`report`，每项输出 `_ms` 与 `_calls`，五项之和为对应 `*_validation_total_ms`。有理路径这些内部明细为零。`SummarizeSphereValidationTiming.py <一轮或多轮 CSV> --baseline <旧完整 CSV> --output <目录>` 仅接受完整的 50 个非有理球面，校验计时分解并比较状态、误差界及工作量；保留既有状态异常的非零退出码。见[三轮测量报告](../docs/testing/SPHERE_VALIDATION_TIMING_2026-10-14.md)。
