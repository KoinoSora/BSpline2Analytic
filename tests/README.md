# 测试代码与固定输入

本目录保存测试实现、案例生成器和数据；测试说明与报告统一放在 [docs/testing/](../docs/testing/TEST_REPORT.md)。

| 文件/目录 | 职责 |
| --- | --- |
| `RecognitionTests.cpp`、`ApproximateConversionTests.cpp` | 有理识别与非有理转换回归。 |
| `PlateauTests.cpp`、`JacobianTests.cpp`、`ConeTests.cpp` | 平台扰动、导数及五类扰动、圆锥专项。 |
| `Dataset500Tests.cpp`、`Dataset500.*`、`PairedSurfaceBounds.*` | 固定数据加载、配额与来源证书校验、数据生成。 |
| `ReferenceCases.*`、`ApproximateCases.*` | 各测试共享的精确与非有理案例。 |
| `ApproximateBenchmark.cpp`、`SummarizeDataset500.py` | 性能记录与 CSV 汇总。 |
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
`ExtremeDataset.*` 和 `ExtremeTests.cpp` 提供独立的 200 例极端输入与 7 个模块诊断。`RunExtremeDataset.py` 隔离逐例执行并记录超时，`SummarizeExtremeDataset.py` 校验完整性、汇总分组耗时。只运行此组：`./Build.ps1 -ExtremeBoundaryOnly`；构造与验收口径见 [EXTREME200.md](../docs/testing/EXTREME200.md)。

`ClassificationGateTests.cpp` 检查短圆弧、数值共线、双候选保留、失败状态与通用回退，并读取原先返回 `InvalidInput` 的 8 个固定输入检查边界拟合。构建目标为 `BSpline2AnalyticClassificationGateTests`；定向运行：`ctest --test-dir build -C Release -R '^BSpline2Analytic\.ClassificationGate$' --output-on-failure`。

极端测试运行器支持 `--first 96 --last 128`，只运行角跨度门槛组（右端不含）。在 SDK 工作目录并配置其 DLL 路径后，执行 `python <repo>/tests/RunExtremeDataset.py --exe <repo>/build/Release/BSpline2AnalyticExtremeTests.exe --data <repo>/tests/data/extreme200 --output <result.csv> --first 96 --last 128`。`--resume` 仅接受同一选定范围的连续前缀。每个子进程额外保存 `case-NNN-fits.csv`，保留专用分支的逐等参线拟合诊断，包括回退前的拒绝原因；采样拟合误差不等于连续界。32 例的几何转换仍有未通过结果时，运行器返回 1，不能把定向单元测试通过等同于整组转换通过。结果见 [门槛修复报告](../docs/testing/CLASSIFICATION_GATE_REPORT_2026-10-14.md)。
