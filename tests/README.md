# 测试代码与固定输入

本目录保存测试实现、案例生成器和数据；测试说明与报告统一放在 [docs/testing/](../docs/testing/TEST_REPORT.md)。

独立混合验证器的单测为 `BSpline2Analytic.GeneralHybridUnit`。批量运行使用真实分类候选，不执行母面优化、不调用旧验证器；源数据证书与输出独立证书在计时外执行。示例（仓库根目录）：

```powershell
cmake -S . -B build-general-hybrid -A x64
cmake --build build-general-hybrid --config Release --target BSpline2AnalyticGeneralHybridTests --parallel 2
ctest --test-dir build-general-hybrid -C Release -R '^BSpline2Analytic.GeneralHybridUnit$' --output-on-failure
python tests/RunGeneralHybrid.py --exe build-general-hybrid/Release/BSpline2AnalyticGeneralHybridTests.exe --sdk ../../../SGGK --data tests/data/dataset500 --output build-general-hybrid/results.csv --runs 3 --timeout 30
```

默认覆盖 600 个主案例和 60 个专项；`--first/--last` 使用半开区间，0–499 为正例、500–599 为负例、600–659 为专项。CSV 保留每次测量；摘要按同例三次中位数统计分位数，另列三次均接受且均不超过 10 ms 的案例数。每例进程看门狗包含数据认证和输出审计；超时阶段、异常及未完成案例单独记录，不能从验收分母中剔除。运行器退出 0 表示完成测量且接受结果全部通过独立审计，不等同于全部正例通过或性能达标。结果合同与局限见[测试报告](../docs/testing/GENERAL_HYBRID_VALIDATION_2026-10-15.md)。

| 文件/目录 | 职责 |
| --- | --- |
| `RecognitionTests.cpp`、`ApproximateConversionTests.cpp` | 有理识别与非有理转换回归。 |
| `PlateauTests.cpp`、`JacobianTests.cpp`、`ConeTests.cpp` | 平台扰动、导数及五类扰动、圆锥专项。 |
| `Dataset500Tests.cpp`、`Dataset500.*`、`PairedSurfaceBounds.*` | 固定数据加载、配额与来源证书校验、数据生成。 |
| `ReferenceCases.*`、`ApproximateCases.*` | 各测试共享的精确与非有理案例。 |
| `ApproximateBenchmark.cpp`、`SummarizeDataset500.py` | 性能记录与 CSV 汇总。 |
| `PairedValidationTests.cpp`、`SummarizePairedValidation.py` | 独立参数对应验证器专项、非有理球面及完整 Dataset500 分流对照；不替换正式转换入口。 |
| `PairedPreparationBenchmark.cpp`、`BenchmarkPairedPreparation.py` | 同一程序对比指定提交的独立验证器、当前首次验证与源面缓存命中；记录准备子项和首次缓存建立成本。 |
| `data/dataset500/` | 500 个固定 `.b2a` 输入及 `manifest.csv`；整理目录不重新生成数据。 |
| `GeneralDataset.*`、`GeneralDatasetCertificate.cpp`、`GeneralDatasetTests.cpp` | 新增 600 个带独立证书的主案例、60 个专项及仅调用一般兜底的测试入口。 |
| `SummarizeGeneralDataset.py` | 校验新 CSV 完整性并按来源类型、测试组汇总结果。 |

新文件 `general-*.b2a`、`general-extra-*.b2a` 与旧 `case-*.b2a` 共存于 `data/dataset500/`，各自使用独立清单；旧测试按 `case-` 前缀统计。构造、证书和运行方式见 [GENERAL600.md](../docs/testing/GENERAL600.md)。

仅复测一般兜底的 500 个正例：`./Build.ps1 -GeneralPositiveOnly`。结果为 `build/general-positive500-Release.csv`，使用 `python tests/SummarizeGeneralDataset.py --positives <csv>` 汇总。测试仍校验来源证书和独立输出证书，并检查边界提取后的端点、导数及细分方向；不执行负例转换。

一般兜底 CSV 中，`raw_candidates` 为五类母面生成的候选总数（去重及分批选择前），`duplicate_candidates` 为评分前剔除的重复候选数，`score_evaluations` 为开始评分的非重复候选数（包括评分失败的调用，不将同一候选粗评和补评重复计数）。三者满足 `raw_candidates = duplicate_candidates + score_evaluations`；汇总脚本同时支持不含这些列的历史 CSV。

新入口 `BSpline2AnalyticGeneralDatasetTests --compare-positive <data-dir> <csv> 0 500` 对每个正例依次运行 `64-full`、`8-full`、`8-staged` 三种配置，逐例轮换六种执行顺序，并分别执行独立输出核验。只共享只读输入；转换内部采样上下文每次新建。输出共 1500 行，使用 `python tests/SummarizeGeneralDataset.py --positives --configuration 8-staged <csv>` 提取、核验和汇总某一配置，也可加 `--merge <csv>` 保存该配置的 500 行。

新增计数：`seed_attempts` 是五类实际尝试组数，`coarse_points` 是配置的粗评点数（0 表示完整评分），`coarse_scores`/`full_scores` 是进入粗评/完整补评的候选数（含失败），`completed_scores` 是完成评分数，`deferred_candidates` 是首批选择后尚未补评的候选数，`backup_candidates` 是后续批次选出数。`interior_point_evaluations`/`boundary_point_evaluations` 统计评分中的点投影求距次数，不含验证和优化。正例入口还用已有正例穷尽候选池，检查粗评点独立性、缓存复用、分批排序、后备候选与尝试预算。

圆环初值现在使用四法线公共交线的直接代数求解，每组最多两条轴线；子午圆拟合额外使用初值生成分区中最多 12 个分散点，不使用独立评分点。CSV 的 `torus_groups` 记录构造组数，`torus_rank_deficient`、`torus_axis_degenerate`、`torus_no_real_axis` 记录整组拒绝原因；`torus_axes`/`torus_rejected_axes` 记录合法/被拒轴线，`torus_circle_failures`/`torus_invalid_radii` 记录圆拟合/半径失败，`torus_candidates`/`torus_dual_groups` 记录去重前的圆环候选数/双候选组数。汇总脚本检查每组最多两轴，以及每条合法轴线对应一次圆构造结果。正例入口额外检查 144 组已知圆环点的排列、法向反向、尺度与小角度构造，并用原失败圆环代表例检查完整评分和分阶段评分的一致性。

对照配置均使用当前实现，不等同于检出历史提交。保留的历史分类比较脚本仍严格要求与旧版候选完全一致；当前去重时机、奇异法向处理和圆环构造器均有变化，不承诺继续满足该历史等价性检查；应使用逐例转换和独立几何验收比较结果。

历史评分前去重的分类阶段配对计时（显式保持 64 次、完整评分配置）：先构建当前 Release 库，再运行 `python tests/BenchmarkGeneralClassification.py --baseline bd3c179`。脚本从指定 Git 提交提取旧的双参数 `ClassifyGeneral`，仅重命名函数，并与当前库链接到同一个测试程序；对 500 个正例交替执行新旧分类，各重复三轮，每次使用独立采样上下文。逐一比较保留候选的顺序、几何参数和评分是否完全一致。结果保存为 `build/general-classification-paired.csv`，临时基线源码和编译产物均在忽略的 `build/` 下；该计时不包含优化或验证。

构建和复现命令从仓库根目录执行，见 [Build.ps1 用法](../README.md#构建与运行)及[数据集说明](../docs/testing/DATASET500.md)。测试通过 CMake 私有引用 `src/` 内部头文件，调用方只依赖 `include/`。
`ExtremeDataset.*` 和 `ExtremeTests.cpp` 提供独立的 200 例极端输入与 7 个模块诊断。`RunExtremeDataset.py` 隔离逐例执行并记录超时，`SummarizeExtremeDataset.py` 校验完整性、汇总分组耗时。`./Build.ps1 -ExtremeBoundaryOnly` 当前运行新的 `data/business-extreme200/` 业务尺寸矩阵，见 [BUSINESS_EXTREME200.md](../docs/testing/BUSINESS_EXTREME200.md)；历史 `data/extreme200/` 与[旧说明](../docs/testing/EXTREME200.md)保留。新数据生成用 `--generate-business`，旧数据生成用 `--generate`，新旧同编号不可混用。完整结果与每组最大分类加优化耗时见[业务尺寸极端测试报告](../docs/testing/BUSINESS_EXTREME200_REPORT_2026-10-14.md)。

`ClassificationGateTests.cpp` 检查短圆弧、数值共线、双候选保留、失败状态与通用回退，并读取原先返回 `InvalidInput` 的 8 个固定输入检查边界拟合。构建目标为 `BSpline2AnalyticClassificationGateTests`；定向运行：`ctest --test-dir build -C Release -R '^BSpline2Analytic\.ClassificationGate$' --output-on-failure`。

极端测试运行器支持 `--first 96 --last 128`，只运行角跨度门槛组（右端不含）。在 SDK 工作目录并配置其 DLL 路径后，执行 `python <repo>/tests/RunExtremeDataset.py --exe <repo>/build/Release/BSpline2AnalyticExtremeTests.exe --data <repo>/tests/data/extreme200 --output <result.csv> --first 96 --last 128`。`--resume` 仅接受同一选定范围的连续前缀。每个子进程额外保存 `case-NNN-fits.csv`，保留专用分支的逐等参线拟合诊断，包括回退前的拒绝原因；采样拟合误差不等于连续界。32 例的几何转换仍有未通过结果时，运行器返回 1，不能把定向单元测试通过等同于整组转换通过。结果见 [门槛修复报告](../docs/testing/CLASSIFICATION_GATE_REPORT_2026-10-14.md)。

## 等参专用路径计时

位置误差界优化沿用比较工具：`python tests/BenchmarkPairedPreparation.py --baseline d345b30 --build-dir build-specialized-timing --output <comparison.csv> --runs 3`。新增位置内部计时 `compose_ms`、`product_ms`、`difference_ms`、`fast_bound_ms`、`tight_bound_ms`，均为 `surface_ms` 的子项；`fast_passes` 和 `tight_checks` 统计快速界成功及紧界调用的节点数。CSV 汇总兼容旧字段集合。实现范围、专项和同批比较见[位置误差界优化报告](../docs/testing/PAIRED_POSITION_2026-10-15.md)。

准备阶段优化的复现方式：先构建当前 Release 库，再执行 `python tests/BenchmarkPairedPreparation.py --baseline 3985a1b --build-dir build-specialized-timing --output <comparison.csv> --runs 3`。脚本提取旧实现到忽略的构建目录，只重命名私有命名空间和公开函数，并增加两个子计时器；两版共用未修改的区间基础运算和同一编译设置。40 个普通非有理球面的三种模式轮换顺序，各测三次；缓存首次填充另列，不能把缓存命中耗时当作首次验证耗时。见[准备阶段优化报告](../docs/testing/PAIRED_PREPARATION_2026-10-15.md)。

`Build.ps1 -Dataset500Only` 当前直接调用 `Recognize`（200 个有理输入）及 `ConvertApproximate`（300 个非有理输入），不进入一般兜底。固定 500 例全部重载并核验来源证书，40 个退化案例保留退化预期。

CSV 标记 `timing_schema=specialized-v2`，分别输出分类计算、初值验证、优化计算、优化后验证及其调用次数。分类总计和优化总计分别包含对应验证；`validate_ms` 为零，不能再把四项明细加到阶段总计上。用 `python tests/SummarizeSpecializedDataset500.py <csv> <目录>` 校验、生成 `summary.json` 与按类型/分支、控制点档分组的 `groups.csv`。旧报告使用旧计时定义，不能直接比较 `classify_ms` 和 `optimize_ms`。

验证内部五步计时标记为 `validation_timing_schema=exclusive-validation-v1`，`initial_*` 和 `refined_*` 分别累计初始候选及优化后验证。后缀为 `prepare`、`surface_bounds`、`subdivide`、`boundary`、`report`，每项输出 `_ms` 与 `_calls`，五项之和为对应 `*_validation_total_ms`。有理路径这些内部明细为零。`SummarizeSphereValidationTiming.py <一轮或多轮 CSV> --baseline <旧完整 CSV> --output <目录>` 仅接受完整的 50 个非有理球面，校验计时分解并比较状态、误差界及工作量；保留既有状态异常的非零退出码。见[三轮测量报告](../docs/testing/SPHERE_VALIDATION_TIMING_2026-10-14.md)。
