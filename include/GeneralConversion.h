/**
 * @file GeneralConversion.h
 * @brief 解析 UV 线性裁剪的历史一般转换接口与共享分类配置。
 * GeneralAnalyticPatch 包含母面及有限裁剪；它与当前兜底的有理参数域边界输出不是同一种表示。
 */
#pragma once
#include "ApproximateConversion.h"
#include <memory>
#include <cstdint>
#include <Geometry/2D/Curve2D.h>

namespace bspline2analytic
{
struct GeneralConversionOptions
{
    ConversionOptions fitting;
    // 每类曲面的固定采样组尝试预算，构造失败也消耗一次尝试。
    unsigned seedAttemptsPerType = 8;
    // 每批优先保留的同类候选数量；当前批全部失败后，仍可继续使用候选池中其他已去重母面。
    unsigned candidatesPerType = 3;
    double similarScoreAbsoluteRatio = .1;
    double similarScoreRelativeTolerance = .05;
    unsigned boundarySamplesPerSpan = 8;
    // 为 0 时直接完整评分；否则先用至多此数量的独立内部评分点粗评，再补评有竞争力的候选。
    unsigned coarseScorePointCount = 16;
};
enum class PCurveFailure
{
    None, NonlinearPCurve, PCurveOutOfRange, NonMonotonePCurve,
    InconsistentPeriodicBranch, UnresolvedDegeneratePCurve,
    InvalidBoundaryLoop, BoundaryGeometryExceeded, ValidationUnresolved, ChartRecoveryUnresolved
};
const char *PCurveFailureName(PCurveFailure);
struct LinearPCurve
{
    // 展开周期后的母面规范 UV，保留源边顺序：Umin/Umax 沿源 v 变化，Vmin/Vmax 沿源 u 变化。
    sggk::Point2D start, end;
    sggk::Point3D projectedStart, projectedEnd;
    bool collapsed3D = false;
    // SDK 母面自身 UV 约定下的有界一次参数曲线。
    sggk::Curve2DPtr sdkCurve;
    sggk::Point2D CalcPoint(double fraction) const;
};
// 解析母面上由四条线性 PCurve 围成的有限简单参数四边形，允许一边坍缩而成为三角形；不持有源样条。
// motherSurface 自身不含裁剪信息。
// 与 AnalyticPatch 一致，规范 UV 与 SDK UV 分别存储，不能混用。
class GeneralAnalyticPatch
{
public:
    const AnalyticParameters &Parameters() const noexcept { return parameters; }
    const sggk::SurfacePtr &MotherSurface() const noexcept { return motherSurface; }
    const std::array<double, 4> &SourceDomain() const noexcept { return domain; }
    const std::array<LinearPCurve, 4> &PCurves() const noexcept { return pcurves; }
    const std::array<sggk::Point2D, 4> &Corners() const noexcept { return corners; }
    // 输入为规范解析 UV，不是源样条的参数坐标。
    sggk::Point3D CalcPoint(double u, double v) const;
    sggk::Point3D CalcBoundaryPoint(BoundaryId edge, double fraction) const;
    GeneralAnalyticPatch(const AnalyticParameters &, const std::array<double, 4> &,
                        const std::array<sggk::Point2D, 4> &,
                        const std::array<LinearPCurve, 4> &);
private:
    AnalyticParameters parameters;
    sggk::SurfacePtr motherSurface;
    std::array<double, 4> domain;
    std::array<sggk::Point2D, 4> corners;
    std::array<LinearPCurve, 4> pcurves;
};
struct GeneralCandidateScore
{
    sggk::SurfaceType type = sggk::SurfaceType::Unknown;
    unsigned degreesOfFreedom = 0;
    double interiorMax = 0, interiorRms = 0, boundaryMax = 0, normalMax = 0;
    unsigned seedGroup = 0; // 对应 classification.seedGroups 中的索引。
};
enum class GeneralSeedFailure
{
    None, InsufficientPoints, RankDeficient, AmbiguousAxis, InvalidParameters,
    NonFinite, ConstructionFailed
};
const char *GeneralSeedFailureName(GeneralSeedFailure);
struct GeneralSeedGroupDiagnostic
{
    sggk::SurfaceType type = sggk::SurfaceType::Unknown;
    unsigned attempt = 0, sampleCount = 0;
    bool multiPoint = false;
    GeneralSeedFailure failure = GeneralSeedFailure::None;
    unsigned generated = 0, duplicates = 0, scored = 0, fullyScored = 0, selected = 0;
    // 相对奇异值与代数残差用于描述初值构造的数值状况，不代表已认证的几何精度；较小但仍通过秩判定的数值不直接等于失败。
    double relativeRank = 1, axisSeparation = 0, constructionResidual = 0;
};
struct GeneralTorusSeedDiagnostics
{
    // 实际尝试构造的四法线组数量，退化组仍消耗一次尝试。
    unsigned sampleGroups = 0;
    unsigned rankDeficientGroups = 0, degenerateAxisGroups = 0, noRealAxisGroups = 0;
    // 通过代数残差检查的有限轴线数量，发生在子午圆拟合之前。
    unsigned axisCount = 0, rejectedAxes = 0;
    unsigned circleFitFailures = 0, invalidRadii = 0;
    // 几何去重与独立评分前的常规环形圆环面候选数量。
    unsigned candidateCount = 0, dualCandidateGroups = 0;
};
struct GeneralClassificationDiagnostics
{
    // 去重与批次筛选之前生成的母面候选总数。
    unsigned rawCandidateCount = 0;
    unsigned duplicateCandidateCount = 0;
    // 进入评分的去重候选数，包括评分失败的候选。
    unsigned scoreEvaluationCount = 0;
    unsigned seedAttemptCount = 0;
    unsigned coarseScoreCount = 0, fullScoreCount = 0;
    unsigned completedScoreCount = 0;
    unsigned deferredCandidateCount = 0; // 首批完成后仍未完成完整评分的候选数。
    unsigned backupCandidateCount = 0, batchCount = 0;
    std::uint64_t interiorPointEvaluations = 0, boundaryPointEvaluations = 0;
    GeneralTorusSeedDiagnostics torus;
    std::vector<GeneralSeedGroupDiagnostic> seedGroups;
};
struct GeneralValidationReport
{
    CheckState geometry = CheckState::NotChecked;
    CheckState boundaryHausdorff = CheckState::NotChecked;
    CheckState pcurveLinearity = CheckState::NotChecked;
    PCurveFailure pcurveFailure = PCurveFailure::None;
    BoundaryId failedBoundary = BoundaryId::Umin;
    double maxSampledPCurveDeviation = 0;
    double distanceUpperBound = std::numeric_limits<double>::infinity();
    double boundaryHausdorffUpperBound = std::numeric_limits<double>::infinity();
    unsigned piecesVisited = 0;
    unsigned chartAttempts = 0;
    unsigned chartIterations = 0;
    bool hasWitness = false;
    std::array<double, 2> witness{};
    std::string message;
    bool Passed() const noexcept
    {
        return geometry == CheckState::Passed && boundaryHausdorff == CheckState::Passed &&
               pcurveLinearity == CheckState::Passed;
    }
};
struct GeneralConversionResult
{
    ConversionStatus status = ConversionStatus::NoCandidate;
    // 仅在独立连续验证通过后填充。
    std::optional<GeneralAnalyticPatch> patch;
    GeneralValidationReport validation;
    std::vector<GeneralCandidateScore> candidates;
    struct Attempt
    {
        sggk::SurfaceType type=sggk::SurfaceType::Unknown;
        PCurveFailure pcurveFailure=PCurveFailure::None;
        BoundaryId failedBoundary=BoundaryId::Umin;
        bool succeeded=false;
        std::string message;
        unsigned seedGroup = 0, optimizationIterations = 0;
    };
    std::vector<Attempt> attempts;
    ConversionStatistics statistics;
    ConversionTiming timing;
    unsigned chartAttempts = 0;
    unsigned chartIterations = 0;
    std::string message;
    GeneralClassificationDiagnostics classification;
};
GeneralConversionResult ConvertGeneral(const sggk::BSplineSurface &,
                                       const GeneralConversionOptions & = {});
}
