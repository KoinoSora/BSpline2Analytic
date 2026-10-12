#pragma once
#include "AnalyticReconstruction.h"
#include "ApproximateConversion.h"
#include "IntervalBounds.h"
#include <set>
#include <tuple>
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
double EdgeCoordinateRaw(const EdgeGeometry &, V local);
double EdgeCoordinate(const EdgeGeometry &, V local, double reference);
std::array<EdgeGeometry, 4> BoundaryMaps(const sggk::BSplineSurface &, const AnalyticPatch &);
struct EdgePoint
{
    V point, tangent, normal;
    double parameter = 0;
};
EdgePoint MapBoundaryPoint(const AnalyticPatch &, const EdgeGeometry &, V world, double fraction);
EdgePoint MapBoundaryPointAndTangent(const Frame3 &, const EdgeGeometry &, double coordinate);
V ProjectInteriorPointLocal(const AnalyticPatch &, V local, const std::array<EdgeGeometry, 4> &,
                            double coneTangent);
void MaterializePatch(AnalyticPatch &);
AnalyticPatch RecoverApproxPatch(const sggk::BSplineSurface &, const AnalyticParameters &);
struct Sample
{
    double u = 0, v = 0;
    V point, tangent, normal;
    int edge = -1;
    double fraction = 0;
    sggk::DerivType sideU = sggk::DerivType::NONE, sideV = sggk::DerivType::NONE;
};
using SampleKey = std::tuple<double, double, int, int, int>;
struct SourceSampleCache
{
    std::vector<Sample> points;
    std::array<std::vector<size_t>, 4> boundaryOrder;
    std::set<SampleKey> keys;
};
struct SampleCache : SourceSampleCache
{
    // Immutable base is shared across candidates; points owns feedback only.
    std::shared_ptr<const SourceSampleCache> base;
    unsigned version = 0;
    SampleCache() = default;
    explicit SampleCache(std::shared_ptr<const SourceSampleCache>);
    size_t size() const;
    const Sample &operator[](size_t) const;
    void RebuildIndices();
};
SampleCache SampleSurface(const sggk::BSplineSurface &, const SamplingOptions &);
unsigned AddWitnessSamples(const sggk::BSplineSurface &, const ConversionValidationReport &,
                           const SamplingOptions &, SampleCache &);
std::vector<AnalyticPatch> ClassifyApproximate(const sggk::BSplineSurface &,
                                               const ConversionOptions &);
struct ResidualBuffer
{
    // Distance and boundary tangent only. Normal errors are acceptance-only.
    std::vector<double> leastSquares, constraints;
    double cost = 0, worst = 0;
};
struct EvaluationContext
{
    Frame3 frame;
    std::array<EdgeGeometry, 4> edges;
    double coneTangent = 0;
};
struct EvaluationWorkspace
{
    EvaluationContext context;
    std::vector<V> local, mappedPoints, mappedTangents;
    std::vector<double> phase;
    ResidualBuffer current, best, trial, plus, minus;
    std::vector<std::vector<double>> jacobian, constraintJacobian;
    unsigned evaluations = 0, resizes = 0;
};
void EvaluateResiduals(const AnalyticPatch &, const SampleCache &, const ValidationOptions &,
                       EvaluationWorkspace &, ResidualBuffer &);
struct ActiveParameters
{
    std::vector<unsigned> indices;
};
ActiveParameters ParametersFor(sggk::SurfaceType, ParameterLevel);
// Linearizes at p, independently of mutable trial geometry in the workspace.
// Returns false if neither side of a required difference can be evaluated.
bool BuildJacobian(const AnalyticPatch &, const SampleCache &, const ConversionOptions &,
                   const ActiveParameters &, double scale, EvaluationWorkspace &,
                   RefinementDiagnostics &);
AnalyticPatch Increment(const AnalyticPatch &, const std::vector<double> &, double scale);
double SampleScale(const AnalyticPatch &, const SampleCache &);
struct RefinementResult
{
    bool improved = false;
    RefinementStopReason stopReason = RefinementStopReason::NoSamples;
};
RefinementResult RefineCandidate(const sggk::BSplineSurface &, AnalyticPatch &, const SampleCache &,
                                 const ConversionOptions &, ConversionStatistics &);
RefinementResult RefineCandidateLevel(AnalyticPatch &, const SampleCache &,
                                      const ConversionOptions &, ConversionStatistics &,
                                      ParameterLevel, unsigned iterationBudget, bool forceAttempt,
                                      EvaluationWorkspace &);
RefinementResult RefinePlaneOffset(AnalyticPatch &, const SampleCache &, const ConversionOptions &,
                                   ConversionStatistics &, EvaluationWorkspace &);
bool RefineAndValidateCandidate(const sggk::BSplineSurface &, AnalyticPatch &,
                                const ConversionOptions &, ConversionResult &,
                                std::shared_ptr<const SourceSampleCache> &,
                                unsigned candidateIndex);
bool SolveLeastSquares(std::vector<std::vector<double>>, std::vector<double>,
                       std::vector<double> &);
void CheckApproximateInput(const sggk::BSplineSurface &, const ValidationOptions &);
double Angle(V, V);
} // namespace bspline2analytic::detail
