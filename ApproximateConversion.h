#pragma once
#include "BSpline2Analytic.h"

namespace bspline2analytic
{
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
    // Used at SizeOrOffset level when progressive variables are enabled.
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
    // Authoritative geometry for validation. Do not edit parameters/domain
    // independently of surface/sdkDomain in a returned successful result.
    AnalyticParameters parameters;
    sggk::SurfacePtr surface;
    // Canonical u=azimuth; v=height (cone: axial height), polar angle
    // (sphere: z=-R*cos(v)), or tube angle (torus). Plane uses Cartesian u,v.
    std::array<double, 4> domain{0, 1, 0, 1}; // umin,umax,vmin,vmax
    ParameterCorrespondence mapping;
    std::optional<sggk::UVRange> sdkDomain;
    int normalSign = 1;
};
struct ConstraintReport
{
    CheckState state = CheckState::NotChecked;
    // Upper bound at inspected point samples; not a measured maximum or RMS.
    double maxSampledUpperBound = 0;
    double upperBound = std::numeric_limits<double>::infinity();
    BoundaryId boundary = BoundaryId::Umin;
    std::array<double, 4> parameterRegion{};
    bool hasWitness = false;
};
struct ConversionValidationReport
{
    ConstraintReport geometry, boundaryTangent, boundaryNormal;
    CheckState correspondence = CheckState::NotChecked;
    unsigned piecesVisited = 0;
    std::string message;
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
    double target = 0; // Only applicable in the constrained phase.
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
             roundIndex = 0; // One-based, assigned by the converter.
    unsigned warmupIterations = 0, constrainedIterations = 0, acceptedSteps = 0;
    unsigned jacobianEvaluations = 0, lineSearchEvaluations = 0;
    unsigned analyticSampleBlocks = 0, fallbackSampleBlocks = 0;
    unsigned nonsmoothConstraintRows = 0, fallbackResidualEvaluations = 0;
    double jacobianMs = 0;
    double initialWorst = 0, bestWorst = 0, bestCost = 0;
    double finalTrust = 0, lastStepNorm = 0;
    RefinementStopReason stopReason = RefinementStopReason::NoSamples;
    std::vector<RefinementIteration> trace; // Populated only when explicitly requested.
};
struct ConversionStatistics
{
    unsigned baseSamplingBuilds = 0, witnessSamplesAdded = 0;
    // Counts generic LM/constrained attempts; scalar work is in refinementRuns.
    unsigned candidateCount = 0, optimizationIterations = 0;
    unsigned refinementRounds = 0, sampleCount = 0, validationPieces = 0;
    std::vector<RefinementDiagnostics> refinementRuns;
};
struct ConversionResult
{
    ConversionStatus status = ConversionStatus::InvalidInput;
    // A failed attempt may retain a diagnostic candidate. Accept only Success.
    std::optional<AnalyticPatch> patch;
    ConversionValidationReport validation;
    ConversionStatistics statistics;
    ConversionTiming timing;
    std::string message;
};
ConversionResult ConvertApproximate(const sggk::BSplineSurface &input,
                                    const ConversionOptions &options = {});
ConversionValidationReport ValidateConversion(const sggk::BSplineSurface &input,
                                              const AnalyticPatch &output,
                                              const ValidationOptions &options = {});
const char *ConversionStatusName(ConversionStatus status);
const char *CheckStateName(CheckState state);
const char *RefinementStopReasonName(RefinementStopReason reason);
} // namespace bspline2analytic
