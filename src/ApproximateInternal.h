/**
 * @file ApproximateInternal.h
 * @brief 近似分类、恢复、采样与优化共用的内部契约。
 * 将候选生成、残差求解和独立认证分离，缓存所有权与反馈样本版本由各类型明确管理。
 */
#pragma once
#include "AnalyticReconstruction.h"
#include "ApproximateConversion.h"
#include "IntervalBounds.h"
#include "SurfaceSampleContext.h"
#include <functional>
#include <set>
#include <tuple>
namespace bspline2analytic::detail
{
ConversionStatus ApproximateFailureStatus(Status) noexcept;
struct CurveFits
{
    BoundaryInfo line, circle;
    bool hasLine = false, hasCircle = false;
};
struct IsoCurveFitOptions
{
    unsigned subdivisionsPerSpan = 2;
    double fastTolerance = 2.5e-6;
    bool forceGeometric = false;
};
CurveFits FitApproximateIsoCurve(const sggk::Curve3DPtr &, BoundaryId, double tolerance,
                                unsigned maxSamples, IsoCurveFitDiagnostic * = nullptr,
                                const std::function<V(double)> &sample = {},
                                const IsoCurveFitOptions & = {});
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
struct ClassificationSamples
{
    std::shared_ptr<const int> stamp = std::make_shared<const int>(0);
    std::array<std::vector<V>, 4> boundaries;
    std::array<std::vector<double>, 4> boundaryParameters;
    std::vector<V> scorePoints;
    V centerPoint, centerDu, centerDv;
    std::array<double, 4> sourceDomain;
};
struct PatchRecoveryTrace
{
    std::array<std::vector<std::array<double, 2>>, 4> paths;
};
struct RecoveredCandidate
{
    AnalyticPatch patch;
    PatchRecoveryTrace trace;
    AnalyticParameters traceParameters;
    bool hasTrace = false;
    std::shared_ptr<const int> sourceStamp;
};
bool SameMotherParameters(const AnalyticParameters &, const AnalyticParameters &);
PatchRecoveryTrace TraceApproximateBoundaries(const ClassificationSamples &,
                                             const AnalyticParameters &);
ClassificationSamples BuildClassificationSamples(const sggk::BSplineSurface &,
                                                 ClassificationDiagnostics * = nullptr,
                                                 SurfaceSampleContext * = nullptr);
AnalyticPatch RecoverApproxPatch(const ClassificationSamples &, const AnalyticParameters &,
                                 PatchRecoveryTrace * = nullptr);
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
    // 不可变基础采样由各候选共享；points 只持有反馈新增点。
    std::shared_ptr<const SourceSampleCache> base;
    unsigned version = 0;
    SampleCache() = default;
    explicit SampleCache(std::shared_ptr<const SourceSampleCache>);
    size_t size() const;
    const Sample &operator[](size_t) const;
    void RebuildIndices();
};
SampleCache SampleSurface(const sggk::BSplineSurface &, const SamplingOptions &,
                          SurfaceSampleContext * = nullptr);
unsigned AddWitnessSamples(const sggk::BSplineSurface &, const ConversionValidationReport &,
                           const SamplingOptions &, SampleCache &, SurfaceSampleContext * = nullptr);
std::vector<AnalyticPatch> ClassifyApproximate(
    const sggk::BSplineSurface &, const ConversionOptions &, ClassificationDiagnostics * = nullptr,
    bool directionalCones = true, SurfaceSampleContext * = nullptr,
    std::vector<RecoveredCandidate> * = nullptr, ClassificationSamples * = nullptr);
struct ApproximateSeed
{
    AnalyticParameters parameters;
    double score = 0;
    bool fromBoundary = false, directionalCone = false;
};
// 可续行的候选搜索：调用方先验证返回的单个候选，失败后再请求下一个。
// 搜索器内部不调用验证器或优化器。
class IsoClassificationSearch
{
  public:
    IsoClassificationSearch(SurfaceSampleContext &, const ConversionOptions &,
                            ClassificationDiagnostics &);
    ~IsoClassificationSearch();
    IsoClassificationSearch(const IsoClassificationSearch &) = delete;
    IsoClassificationSearch &operator=(const IsoClassificationSearch &) = delete;
    std::optional<RecoveredCandidate> Next(ClassificationSamples &);
  private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
// 保留原来的最佳候选和边界候选，作为新增圆锥初值的后备。
std::vector<AnalyticPatch> RecoverRankedCandidates(const ClassificationSamples &,
                                                   const std::vector<ApproximateSeed> &,
                                                   double distanceTolerance,
                                                   ClassificationDiagnostics * = nullptr,
                                                   std::vector<RecoveredCandidate> * = nullptr);
struct ResidualBuffer
{
    // 优化残差仅含距离与边界切向；法向误差只参与验收。
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
// 几何适配器提供残差、导数和参数更新；迭代推进、步接受与停滞判定复用历史求解框架。
struct RefinementProblem
{
    std::function<void(const AnalyticPatch &, EvaluationWorkspace &, ResidualBuffer &)> evaluate;
    std::function<bool(const AnalyticPatch &, const ActiveParameters &, double,
                       EvaluationWorkspace &, RefinementDiagnostics &)> jacobian;
    std::function<AnalyticPatch(const AnalyticPatch &, const std::vector<double> &, double)> increment;
    std::function<void(AnalyticPatch &)> materialize;
};
ActiveParameters ParametersFor(sggk::SurfaceType, ParameterLevel);
// 在 p 处线性化，不依赖工作区中的可变试探几何；所需差分两侧均不能求值时返回 false。
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
                                      EvaluationWorkspace &, const RefinementProblem * = nullptr);
RefinementResult RefinePlaneOffset(AnalyticPatch &, const SampleCache &, const ConversionOptions &,
                                   ConversionStatistics &, EvaluationWorkspace &);
bool RefineAndValidateCandidate(const sggk::BSplineSurface &, AnalyticPatch &,
                                const ConversionOptions &, ConversionResult &,
                                std::shared_ptr<const SourceSampleCache> &,
                                unsigned candidateIndex, SurfaceSampleContext * = nullptr);
bool SolveLeastSquares(std::vector<std::vector<double>>, std::vector<double>,
                       std::vector<double> &);
void CheckApproximateInput(const sggk::BSplineSurface &, const ValidationOptions &,
                           bool checkNormalTolerance = true);
ConversionResult ConvertApproximateWithSamples(const sggk::BSplineSurface &,
                                               const ConversionOptions &, SurfaceSampleContext &);
double Angle(V, V);
} // namespace bspline2analytic::detail
