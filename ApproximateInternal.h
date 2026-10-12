#pragma once
#include "AnalyticReconstruction.h"
#include "ApproximateConversion.h"
#include "IntervalBounds.h"
namespace bspline2analytic::detail
{
struct Frame3
{
    V o, x, y, z;
    V local(V p) const
    {
        p = p - o;
        return {Dot(p, x), Dot(p, y), Dot(p, z)};
    }
    V global(V p) const
    {
        return o + x * p.x + y * p.y + z * p.z;
    }
    V vector(V p) const
    {
        return x * p.x + y * p.y + z * p.z;
    }
};
Frame3 FrameOf(const AnalyticParameters &);
void MakeNative(AnalyticPatch &);
V CanonicalPoint(const AnalyticParameters &, double, double);
V CanonicalNormal(const AnalyticParameters &, double, double);
std::array<double, 2> MotherParameters(const AnalyticParameters &, V local);
V ProjectInteriorPoint(const AnalyticPatch &, V world);
struct EdgeGeometry
{
    unsigned fixed = 0;
    double value = 0, start = 0, end = 1;
    bool circle = false;
    V o, d, ex, ey, n0, nx, ny;
    double radius = 0, g0 = 0, g1 = 1;
};
EdgeGeometry MakeEdge(const AnalyticPatch &, unsigned);
double EdgeCoordinate(const EdgeGeometry &, V local, double reference);
std::array<EdgeGeometry, 4> BoundaryMaps(const sggk::BSplineSurface &, const AnalyticPatch &);
struct EdgePoint
{
    V point, tangent, normal;
    double parameter = 0;
};
EdgePoint MapBoundaryPoint(const AnalyticPatch &, const EdgeGeometry &, V world, double fraction);
AnalyticPatch RecoverApproxPatch(const sggk::BSplineSurface &, const AnalyticParameters &);
struct Sample
{
    double u = 0, v = 0;
    V point, tangent, normal;
    int edge = -1;
    double fraction = 0;
    sggk::DerivType sideU = sggk::DerivType::NONE, sideV = sggk::DerivType::NONE;
};
struct SampleCache
{
    std::vector<Sample> points;
};
SampleCache SampleSurface(const sggk::BSplineSurface &, const SamplingOptions &);
unsigned AddWitnessSamples(const sggk::BSplineSurface &, const ConversionValidationReport &,
                           const SamplingOptions &, SampleCache &);
std::vector<AnalyticPatch> ClassifyApproximate(const sggk::BSplineSurface &,
                                               const ConversionOptions &);
struct RefinementResult
{
    bool improved = false;
    RefinementStopReason stopReason = RefinementStopReason::NoSamples;
};
RefinementResult RefineCandidate(const sggk::BSplineSurface &, AnalyticPatch &, const SampleCache &,
                                 const ConversionOptions &, ConversionStatistics &);
bool SolveLeastSquares(std::vector<std::vector<double>>, std::vector<double>,
                       std::vector<double> &);
void CheckApproximateInput(const sggk::BSplineSurface &, const ValidationOptions &);
double Angle(V, V);
} // namespace bspline2analytic::detail
