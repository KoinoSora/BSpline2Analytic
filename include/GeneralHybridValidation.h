#pragma once
#include "GeneralConversion.h"

namespace bspline2analytic
{
// Standalone experimental validator. Interior acceptance is sampled evidence,
// NOT a continuous surface certificate. No optimizer or legacy validator runs.
enum class HybridValidationStatus { Passed, NotCertified, InvalidInput, Unsupported, BudgetExceeded };
struct GeneralHybridOptions
{
    GeneralConversionOptions general;
    unsigned maxBoundaryPieces = 2048;
    unsigned maxBoundaryDepth = 12;
};
struct GeneralHybridReport
{
    static constexpr bool certifiesInteriorContinuum = false;
    HybridValidationStatus status = HybridValidationStatus::NotCertified;
    CheckState interiorSamples = CheckState::NotChecked;
    CheckState boundaryPosition = CheckState::NotChecked;
    CheckState boundaryTangent = CheckState::NotChecked;
    // On successful validation, true marks an exactly point-collapsed SOURCE
    // boundary. Its tangent is undefined and exempt, not certified as zero error.
    std::array<bool,4> tangentNotApplicable{};
    // Deliberately not named distanceUpperBound: no interior continuum claim.
    double interiorSampleMaxDistance = 0;
    double interiorSampleMaxFiniteDistance = 0;
    double boundaryDistanceUpperBound = INFINITY;
    double boundaryTangentUpperBoundDeg = INFINITY;
    unsigned interiorSampleCount = 0, boundaryNodes = 0, boundarySegments = 0;
    unsigned mappingRefinements = 0, collapsedEdges = 0, chartAttempts = 0, chartIterations = 0;
    unsigned failedBoundary = 0;
    double samplingMs = 0, recoveryMs = 0, preparationMs = 0, boundaryMs = 0, totalMs = 0;
    std::optional<GeneralAnalyticPatch> patch;
    std::string message;
    bool Passed() const noexcept
    {
        return status == HybridValidationStatus::Passed && patch.has_value() &&
               interiorSamples == CheckState::Passed && boundaryPosition == CheckState::Passed &&
               boundaryTangent == CheckState::Passed;
    }
};
// Retired on 2026-10-15. Preserve code; do not run tests or use as a benchmark.
// Current validator: ValidateGeneralRationalBoundaries. Shared status stays usable.
GeneralHybridReport ValidateGeneralHybrid(const sggk::BSplineSurface &, const AnalyticParameters &,
                                         const GeneralHybridOptions & = {});
const char *HybridValidationStatusName(HybridValidationStatus) noexcept;
}
