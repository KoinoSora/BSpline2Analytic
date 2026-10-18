/**
 * @file Dataset500.h
 * @brief 固定等参数据及共享文件格式接口。
 * 输入、真值和数据资格信息分开保存，生产转换器只接收输入曲面。
 */
#pragma once
#include "PairedSurfaceBounds.h"
#include <filesystem>
namespace bspline2analytic::tests
{
struct DatasetCase
{
    unsigned id = 0, pair = 0, type = 0, bucket = 0, variant = 0, seed = 0;
    bool refitted = false;
    unsigned collapsedMask = 0;
    AnalyticParameters truth;
    std::array<double, 4> analyticDomain{};
    sggk::BSplineSurfacePtr source, input;
    double storedBound = 0;
};
sggk::SurfaceType DatasetType(unsigned type);
std::string CaseFilename(unsigned id);
void GenerateDataset500(const std::filesystem::path &directory);
// 此限制仅适用于 Dataset500；共享读取器也服务其他高次测试输入。
void CheckDataset500Contract(const DatasetCase &, unsigned expectedId);
void AuditDataset500(const std::filesystem::path &directory);
void RebuildLowDegreeDataset500(const std::filesystem::path &inputDirectory,
                                const std::filesystem::path &outputDirectory);
DatasetCase ReadDatasetCase(const std::filesystem::path &file);
unsigned CollapsedMask(const sggk::BSplineSurface &s);
sggk::BSplineSurfacePtr RebuildDatasetSource(unsigned type, unsigned bucket, unsigned variant);
// 共享测试数据工具，本身不加载或运行历史 500 例套件。
void SaveDatasetCase(const DatasetCase &, const std::filesystem::path &);
sggk::BSplineSurfacePtr NormalizeDatasetSurface(const sggk::BSplineSurface &);
sggk::BSplineSurfacePtr FitDatasetPolynomial(const sggk::BSplineSurface &, unsigned, unsigned);
void PadDatasetSurface(sggk::BSplineSurface &, unsigned, unsigned, unsigned);
} // namespace bspline2analytic::tests
