#pragma once
#include "GeneralConversion.h"
#include "IntervalBounds.h"
#include <filesystem>

namespace bspline2analytic::tests
{
// Files contain input geometry and an independent oracle; only input reaches
// ConvertGeneral. Extra cases have their own expectations, never a no-solution label.
enum class GeneralFixtureGroup { Feasible, Infeasible, NonlinearPCurve, Periodic,
    SpherePriority, Tolerance, InvalidOptions, ValidationBudget };
struct GeneralFixture
{
    unsigned id = 0, type = 0, bucket = 0, variant = 0, seed = 0;
    bool rational = false, extra = false;
    unsigned collapsedMask = 0;
    GeneralFixtureGroup group = GeneralFixtureGroup::Feasible;
    AnalyticParameters truth;
    std::array<sggk::Point2D, 4> corners;
    // Optional curved-UV test: v += bend*s*(1-s). The 500 positives use zero.
    double bend = 0;
    sggk::BSplineSurfacePtr input;
    double storedUpper = 0, exclusionRadius = 0, storedRankMargin = 0, pairTolerance = 1e-6;
    // V2: corners belong to a specified rational mother chart, not analytic UV.
    bool rationalDomain = false;
    std::array<double,4> chartDomain{};
    double storedTangentDeg = INFINITY;
};
struct FixtureCertificate
{
    bool passed = false;
    double upper = INFINITY, rankMargin = 0;
    unsigned pieces = 0, witnesses = 0;
    std::string message;
};
sggk::SurfaceType GeneralDatasetType(unsigned);
const char *FixtureGroupName(GeneralFixtureGroup);
std::string GeneralFixtureFilename(unsigned, bool extra = false);
GeneralFixture ReadGeneralFixture(const std::filesystem::path &);
void WriteGeneralFixture(const GeneralFixture &, const std::filesystem::path &);
void GenerateGeneralDataset(const std::filesystem::path &, unsigned first = 0, unsigned last = 660);
// Canonical analytic formula and chart are test-owned, not production projection code.
detail::V FixtureTarget(const GeneralFixture &, double, double);
FixtureCertificate CertifyAnalyticPair(const GeneralFixture &, double tolerance,
                                      unsigned maxPieces = 20000);
// Excludes EVERY nonzero polynomial of total degree <= 4, hence all five mothers.
FixtureCertificate CertifyNoQuartic(const sggk::BSplineSurface &, double radius);
unsigned FixtureCollapsedMask(const sggk::BSplineSurface &);
}
