/**
 * @file ApproximateConversion.h
 * @brief 近似转换公共类型、优化配置与历史投影验证接口。
 * 这些数据类型仍被当前实现复用；ConvertApproximate/ValidateConversion 的完整历史验证链已禁测，当前串联接口见
 * SequentialConversion.h。
 */
#pragma once
#include "BSpline2Analytic.h"
#include <vector>

namespace bspline2analytic
{
/// 检查状态：未检查、通过、有超限结论、无法判定；Unresolved 不能当成几何不可行证明。
enum class CheckState
{
    NotChecked,
    Passed,
    Exceeded,
    Unresolved
};
enum class ConversionStatus
{
    Success,
    InvalidInput,
    Unsupported,
    NoCandidate,
    Degenerate,
    OptimizationFailed,
    VerificationFailed,
    Unresolved
};
/// 共享验证选项：长度跟随输入单位，角度显式使用度；各验证器实际采用的字段见其接口合同。
struct ValidationOptions
{
    double distanceTolerance = 1e-5;
    double boundaryTangentToleranceDeg = 0.1;
    double boundaryNormalToleranceDeg = 0.1;
    unsigned maxSubdivisionDepth = 18;
    unsigned maxValidationPieces = 40000;
};
struct SamplingOptions
{
    unsigned interiorIsoCurves = 3;
    unsigned interiorPointsPerSpan = 2;
    unsigned boundaryPointsPerSpan = 5;
    unsigned maxSamples = 20000;
};
enum class JacobianMode
{
    AnalyticWithFallback,
    FiniteDifference
};
/// 采样优化的迭代和停滞配置；放宽候选筛选或改变迭代预算不会放宽独立验证的几何容差。
struct RefinementOptions
{
    JacobianMode jacobianMode = JacobianMode::AnalyticWithFallback;
    unsigned maxIterations = 30;
    unsigned maxRounds = 5;
    double candidateToleranceFactor = 20;
    unsigned maxWarmupIterations = 5;
    unsigned stagnationWindow = 3;
    double worstImprovementAbsTol = 1e-5;
    double worstImprovementRelTol = 1e-4;
    double costImprovementRelTol = 1e-4;
    unsigned maxConsecutiveFailures = 3;
    double minTrustRadius = 1e-9;
    double smallModelChangeTolerance = 1e-3;
    bool recordIterationTrace = false;
    bool enableProgressiveVariables = true;
    // 启用渐进变量时，仅在 SizeOrOffset 层使用该平面偏移求解器。
    bool enablePlaneOffsetSolver = true;
    unsigned maxScalarIterations = 48;
    double scalarDistanceToleranceRatio = 1e-4;
};
struct ConversionOptions
{
    ValidationOptions validation;
    SamplingOptions sampling;
    RefinementOptions refinement;
};
struct AnalyticPatch
{
    // 验证以此几何参数为准。
    // 对于成功返回的结果，修改 parameters/domain 时必须同步 surface/sdkDomain，避免两套表示不一致。
    AnalyticParameters parameters;
    sggk::SurfacePtr surface;
    // 规范参数 u 为方位角；v 为高度（圆锥采用轴向高度）、球面极角（z=-R*cos(v)）或圆环面管截面角。
    // 平面采用笛卡尔 u、v。
    std::array<double, 4> domain{0, 1, 0, 1}; // 参数域顺序：umin、umax、vmin、vmax。
    ParameterCorrespondence mapping;
    std::optional<sggk::UVRange> sdkDomain;
    int normalSign = 1;
};
struct ConstraintReport
{
    CheckState state = CheckState::NotChecked;
    // 已检查采样点处的误差上界，不是整个连续域的最大值，也不是均方根误差。
    double maxSampledUpperBound = 0;
    double upperBound = std::numeric_limits<double>::infinity();
    BoundaryId boundary = BoundaryId::Umin;
    std::array<double, 4> parameterRegion{};
    bool hasWitness = false;
};
// 验证器内部互斥的墙钟计时阶段；五阶段之和为 totalMs。
// 嵌套和递归工作只归属当前活动阶段，避免重复计时。
struct ValidationTiming
{
    StageTiming preparation, surfaceBounds, subdivision, boundaryChecks, reporting;
    double totalMs = 0;
};
struct ConversionValidationReport
{
    ConstraintReport geometry, boundaryTangent, boundaryNormal;
    CheckState correspondence = CheckState::NotChecked;
    unsigned piecesVisited = 0;
    std::string message;
    ValidationTiming timing;
    bool Passed() const noexcept
    {
        return geometry.state == CheckState::Passed &&
               boundaryTangent.state == CheckState::Passed &&
               boundaryNormal.state == CheckState::Passed && correspondence == CheckState::Passed;
    }
};
enum class RefinementPhase
{
    Warmup,
    Constrained
};
enum class ParameterLevel
{
    SizeOrOffset,
    MotherSurface,
    FullPatch
};
enum class RefinementEntryReason
{
    Initial,
    WitnessSamples,
    MoreVariables
};
enum class RefinementStopReason
{
    NoSamples,
    MarginReached,
    SampleFeasibleStagnation,
    Stagnation,
    NoAcceptedStep,
    TrustTooSmall,
    SmallStep,
    IterationLimit,
    SolverFailure,
    PlaneOffsetSolved,
    ScalarIterationLimit
};
struct RefinementIteration
{
    unsigned iteration = 0;
    RefinementPhase phase = RefinementPhase::Warmup;
    double worst = 0, cost = 0, bestWorst = 0;
    double target = 0; // 仅约束优化阶段使用此目标值。
    double trustBefore = 0, trustAfter = 0, stepNorm = 0, modelChange = 0;
    bool accepted = false, stepSolved = false;
};
struct RefinementDiagnostics
{
    ParameterLevel level = ParameterLevel::FullPatch;
    RefinementEntryReason entryReason = RefinementEntryReason::Initial;
    unsigned activeDimension = 0, sampleVersion = 0;
    unsigned planeOffsetSolves = 0, scalarIterations = 0;
    unsigned residualEvaluations = 0, workspaceResizes = 0;
    double planeOffset = 0;
    unsigned candidateIndex = 0,
             roundIndex = 0; // 由转换器赋值，编号从 1 开始。
    unsigned warmupIterations = 0, constrainedIterations = 0, acceptedSteps = 0;
    unsigned jacobianEvaluations = 0, lineSearchEvaluations = 0;
    unsigned analyticSampleBlocks = 0, fallbackSampleBlocks = 0;
    unsigned nonsmoothConstraintRows = 0, fallbackResidualEvaluations = 0;
    double jacobianMs = 0;
    double initialWorst = 0, bestWorst = 0, bestCost = 0;
    double finalTrust = 0, lastStepNorm = 0;
    RefinementStopReason stopReason = RefinementStopReason::NoSamples;
    std::vector<RefinementIteration> trace; // 仅在显式要求记录迭代轨迹时填充。
};
enum class CircleFitState
{
    NotAttempted, Accepted, NumericallyLinear, RankDeficient,
    NonFinite, InvalidRadius, UnresolvedExtent, ResidualExceeded
};
const char *CircleFitStateName(CircleFitState);
struct IsoCurveFitDiagnostic
{
    unsigned fixedDirection = 0; // 0 表示固定源曲面的 U 参数，1 表示固定 V 参数。
    double fixedParameter = 0;
    bool boundary = false;
    BoundaryId edge = BoundaryId::Umin;
    unsigned sampleCount = 0;
    double candidateTolerance = 0, numericFloor = 0, transverseRatio = 0;
    double lineError = std::numeric_limits<double>::infinity();
    double circleError = std::numeric_limits<double>::infinity();
    double radius = 0, sweepAngle = 0, sagitta = 0;
    bool hasLine = false, hasCircle = false;
    unsigned subdivisionsPerSpan = 0, geometricIterations = 0;
    bool skippedGeometricFit = false;
    double samplingMs = 0, algebraicMs = 0, geometricMs = 0;
    CircleFitState circleState = CircleFitState::NotAttempted;
};
struct ClassificationDiagnostics
{
    unsigned sampleCacheBuilds = 0;
    // 分类缓存发起的源曲面求值次数，不包括等参曲线拟合的曲线求值。
    unsigned cachedPointEvaluations = 0, cachedDerivativeEvaluations = 0;
    unsigned rawCandidateCount = 0, directionalConeCandidates = 0;
    unsigned patchRecoveryAttempts = 0;
    // 增量等参分类器各互斥阶段的耗时；验证由调用方计时，不计入这里的统计。
    double extractionSamplingMs = 0, algebraicFitMs = 0, geometricFitMs = 0;
    double candidateScoringMs = 0, patchRecoveryMs = 0;
    unsigned curveExtractions = 0, curvePointEvaluations = 0, sourcePointEvaluations = 0;
    unsigned pointCacheHits = 0, fitPasses = 0, geometricIterations = 0;
    unsigned fastCircleFits = 0, samplingRefinements = 0, classificationStages = 0;
    unsigned distanceEvaluations = 0, sharedMotherScores = 0, duplicateCandidates = 0;
    // 采样拟合诊断，不构成连续几何误差证书。
    std::vector<IsoCurveFitDiagnostic> isoCurveFits;
};
struct ConversionStatistics
{
    ClassificationDiagnostics classification;
    unsigned baseSamplingBuilds = 0, witnessSamplesAdded = 0;
    // 统计通用 LM 与约束步的尝试次数；标量求解工作另见 refinementRuns。
    unsigned candidateCount = 0, optimizationIterations = 0;
    unsigned refinementRounds = 0, sampleCount = 0, validationPieces = 0;
    std::vector<RefinementDiagnostics> refinementRuns;
    // 累计所有候选与所有轮次的验证工作，不仅记录最后一次验证。
    ValidationTiming initialValidationTiming, refinementValidationTiming;
};
struct ConversionResult
{
    ConversionStatus status = ConversionStatus::InvalidInput;
    // 失败时可能保留候选用于诊断，只有 status 为 Success 才可接受输出。
    std::optional<AnalyticPatch> patch;
    ConversionValidationReport validation;
    ConversionStatistics statistics;
    ConversionTiming timing;
    std::string message;
    ClassificationFailureReason classificationFailure = ClassificationFailureReason::None;
};
/// 历史完整近似转换入口，仍使用旧投影验证链；当前禁测，限时转换请用 ConvertIsoThenFallback。
ConversionResult ConvertApproximate(const sggk::BSplineSurface &input,
                                    const ConversionOptions &options = {});
/// 已废弃的投影验证接口，仅保留兼容与历史实现，不用于新回归或性能比较。
ConversionValidationReport ValidateConversion(const sggk::BSplineSurface &input,
                                              const AnalyticPatch &output,
                                              const ValidationOptions &options = {});
const char *ConversionStatusName(ConversionStatus status);
const char *CheckStateName(CheckState state);
const char *RefinementStopReasonName(RefinementStopReason reason);
} // namespace bspline2analytic
