#pragma once
#include "GeneralConversion.h"
#include "ApproximateInternal.h"
#include <functional>
namespace bspline2analytic::detail
{
V ProjectToMother(const AnalyticParameters &, V);
struct GeneralSeed
{
    AnalyticParameters parameters;
    GeneralCandidateScore score;
};
struct GeneralTrim
{
    AnalyticParameters parameters;
    std::array<sggk::Point2D, 4> corners; // 00,10,01,11
    std::array<LinearPCurve, 4> edges;
    // Ordered projected UV samples fix the periodic branch for each source edge.
    std::array<std::vector<std::pair<double, std::array<double, 2>>>, 4> walks;
};
bool TryPairedTrimBound(const sggk::BSplineSurface &, const std::vector<BezierCell> &,
                       const GeneralTrim &, const ValidationOptions &, GeneralValidationReport &);
bool RecoverLinearPCurves(SurfaceSampleContext &, const GeneralConversionOptions &,
                         const AnalyticParameters &, GeneralTrim &, GeneralValidationReport &,
                         const std::function<bool(const GeneralTrim &,GeneralValidationReport &)> & = {});
std::vector<AnalyticParameters> RefinedSphereCharts(SurfaceSampleContext &, const GeneralConversionOptions &,
    const AnalyticParameters &, const std::vector<AnalyticParameters> &, unsigned &);
std::array<double, 2> TrimScale(const GeneralTrim &);
bool InsideTrim(const GeneralTrim &, std::array<double, 2>, double tolerance = 0);
std::vector<std::array<sggk::Point2D,3>> TrimTriangles(const GeneralTrim &);
std::array<double, 2> NearestTrimUV(const GeneralTrim &, std::array<double, 2>);
std::array<double, 2> EdgeUV(const LinearPCurve &, double);
double EdgeReference(const GeneralTrim &, unsigned, double, unsigned);
std::vector<GeneralSeed> ClassifyGeneral(SurfaceSampleContext &, const GeneralConversionOptions &);
GeneralValidationReport ValidateGeneral(const sggk::BSplineSurface &, const AnalyticParameters &,
                                        const GeneralConversionOptions &, SurfaceSampleContext &,
                                        GeneralTrim * = nullptr);
GeneralConversionResult ConvertGeneralWithSamples(const sggk::BSplineSurface &,
                                                  const GeneralConversionOptions &,
                                                  SurfaceSampleContext &);
RefinementProblem MotherRefinementProblem(const SampleCache &, const ConversionOptions &);
}
