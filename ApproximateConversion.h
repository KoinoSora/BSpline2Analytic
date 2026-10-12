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
struct StageTiming
{
    double elapsedMs = 0;
    unsigned calls = 0;
};
struct ConversionTiming
{
    StageTiming classification, optimization, validation;
    double otherMs = 0, totalMs = 0;
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
struct RefinementOptions
{
    unsigned maxIterations = 30;
    unsigned maxRounds = 5;
    double candidateToleranceFactor = 20;
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
struct ConversionStatistics
{
    unsigned candidateCount = 0, optimizationIterations = 0;
    unsigned refinementRounds = 0, sampleCount = 0, validationPieces = 0;
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
} // namespace bspline2analytic
