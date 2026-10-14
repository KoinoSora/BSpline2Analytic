#pragma once
#include "Dataset500.h"
#include "Numerics.h"
#include <string>
namespace bspline2analytic::tests
{
struct ExtremeSpec
{
    unsigned pair = 0, type = 0, orientation = 0;
    std::string group, label;
    AnalyticParameters truth;
    std::array<double, 4> domain{};
    double rho = 0, theta = 0, knotRatio = 1;
    bool collapsed = false;
};
struct ExtremeQualification
{
    bool passed = false, dimensionsPassed = true;
    PairBound pair;
    unsigned controls = 0;
    std::array<double, 4> lengths{};
    double lengthDisagreement = 0, minWeight = 1, maxWeight = 1, ulp = 0;
    std::array<double,2> sourceAngles{};
    std::string message;
};
std::vector<ExtremeSpec> ExtremeSpecs();
std::string ExtremeFilename(unsigned id);
std::vector<ExtremeSpec> BusinessExtremeSpecs();
bool IsBusinessExtreme(const std::filesystem::path &);
bool BusinessDimensions(const DatasetCase &);
std::array<double,2> BusinessBoundaryAngles(const sggk::BSplineSurface &, const sggk::BSplineSurface &);
void PadBusinessExtreme(sggk::BSplineSurface &, unsigned bucket);
void GenerateExtremeDataset(const std::filesystem::path &, bool business = false);
ExtremeQualification QualifyExtreme(const DatasetCase &, bool business = false);
// Independent test formulas for mother distance/normal (no production projection).
double ExtremeMotherDistance(const AnalyticParameters &, detail::V);
detail::V ExtremeMotherNormal(const AnalyticParameters &, detail::V);
} // namespace bspline2analytic::tests
