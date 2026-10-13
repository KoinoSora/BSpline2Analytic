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
DatasetCase ReadDatasetCase(const std::filesystem::path &file);
unsigned CollapsedMask(const sggk::BSplineSurface &s);
sggk::BSplineSurfacePtr RebuildDatasetSource(unsigned type, unsigned bucket, unsigned variant);
// Shared fixture utilities. These do not load or run the legacy 500-case suite.
void SaveDatasetCase(const DatasetCase &, const std::filesystem::path &);
sggk::BSplineSurfacePtr NormalizeDatasetSurface(const sggk::BSplineSurface &);
sggk::BSplineSurfacePtr FitDatasetPolynomial(const sggk::BSplineSurface &, unsigned, unsigned);
void PadDatasetSurface(sggk::BSplineSurface &, unsigned, unsigned, unsigned);
} // namespace bspline2analytic::tests
