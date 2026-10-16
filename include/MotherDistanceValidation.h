#pragma once
#include "PairedConversionValidation.h"

namespace bspline2analytic
{
// Independent experimental validator for polynomial sources. Full-domain distance
// uses the exact-recognition settings and a floating-point Bernstein guard.
// Four paired boundary positions and oriented tangents use intervals.
// This is NOT an interior paired-distance or trimmed Hausdorff certificate.
struct MotherDistanceValidationOptions
{
    RecognitionOptions position;
    // position.distanceTolerance also controls boundary position and mapping.
    // boundary.validation.distanceTolerance and boundaryNormalToleranceDeg are ignored;
    // tangent/budget fields apply. Surface normal orientation is not certified.
    PairedValidationOptions boundary;
};
struct MotherDistanceValidationReport
{
    PairedValidationStatus status = PairedValidationStatus::NotCertified;
    ConstraintReport motherDistance, boundaryPosition, boundaryTangent;
    CheckState correspondence = CheckState::NotChecked;
    ValidationTiming timing;
    std::string message;
    PairedParameterMap uMap, vMap;
    unsigned surfaceNodes = 0, surfaceDepth = 0, boundaryPositionNodes = 0;
    unsigned boundaryNodes = 0, mappingAttempts = 0, sampleCount = 0;
    unsigned sourceExtractions = 0, sourceCacheHits = 0;
    double sourceExtractionMs = 0, sourceBoundaryMs = 0, targetPreparationMs = 0;
    double mappingMs = 0, sampleMs = 0, residualBuildMs = 0, residualBoundMs = 0;
    double boundaryPositionMs = 0, boundaryAngleMs = 0; // Angle timing is tangent-only.
    double sourceEnclosureRadius = 0;
    // Streaming extraction diagnostics. All source times are preparation subsets.
    unsigned sourceCells = 0, extractedCells = 0, cachedCellHits = 0;
    unsigned extractionPlanBuilds = 0, tightenedCells = 0, boundaryFallbacks = 0;
    double extractionPlanMs = 0, sourceContractionMs = 0, sourceEnvelopeMs = 0;
    double maxExtractionError = 0;
    bool Passed() const noexcept
    {
        return status == PairedValidationStatus::Passed &&
               motherDistance.state == CheckState::Passed &&
               boundaryPosition.state == CheckState::Passed &&
               boundaryTangent.state == CheckState::Passed && correspondence == CheckState::Passed;
    }
};
MotherDistanceValidationReport ValidateMotherDistanceConversion(
    const sggk::BSplineSurface &, const AnalyticPatch &,
    const MotherDistanceValidationOptions & = {});
} // namespace bspline2analytic
