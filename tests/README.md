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

构建和复现命令从仓库根目录执行，见 [Build.ps1 用法](../README.md#构建与运行)及[数据集说明](../docs/testing/DATASET500.md)。测试通过 CMake 私有引用 `src/` 内部头文件，调用方只依赖 `include/`。
