#pragma once
#include "PairedConversionValidation.h"

namespace bspline2analytic
{
// Independent experimental validator for polynomial sources. Full-domain distance
// uses the exact-recognition settings and a floating-point Bernstein guard.
// Four paired boundary positions and both parameter tangents use intervals.
// This is NOT an interior paired-distance or trimmed Hausdorff certificate.
struct MotherDistanceValidationOptions
{
    RecognitionOptions position;
    // position.distanceTolerance also controls boundary position and mapping.
    // The legacy distance/tangent/normal tolerances in boundary.validation are
    // ignored. Its subdivision and mapping budgets still apply.
    PairedValidationOptions boundary;
    // Require sin^2(angle(Su,Qu)) + sin^2(angle(Sv,Qv)) < sin^2(tolerance)
    // on every boundary interval, with positive direction dot products.
    // Target parameter directions are analytically orthogonal. No normal field
    // is constructed, and failure never falls back to direct normal validation.
    double jointTangentToleranceDeg = 0.1;
};
struct MotherDistanceValidationReport
{
    PairedValidationStatus status = PairedValidationStatus::NotCertified;
    ConstraintReport motherDistance, boundaryPosition;
    ConstraintReport boundaryTangentU, boundaryTangentV, boundaryJointTangent;
    // Joint bound is max over shared leaf intervals, NOT the sum of separately
    // maximized U/V bounds. boundaryJointTangent.upperBound is in degrees and
    // bounds asin(sqrt(the joint sine-squared budget)), not a measured normal.
    double jointTangentSineSquaredUpperBound = INFINITY;
    // Empty on success/not checked. Stable failure codes distinguish a witnessed
    // joint-condition violation from an inconclusive or budget-limited bound.
    std::string jointTangentFailure;
    CheckState correspondence = CheckState::NotChecked;
    ValidationTiming timing;
    std::string message;
    PairedParameterMap uMap, vMap;
    unsigned surfaceNodes = 0, surfaceDepth = 0, boundaryPositionNodes = 0;
    unsigned boundaryNodes = 0, mappingAttempts = 0, sampleCount = 0;
    unsigned sourceExtractions = 0, sourceCacheHits = 0;
    // Positive-scale tangent residual fits (two per joint node) and joint splits.
    unsigned tangentResidualFits = 0, tangentResidualSplits = 0;
    double sourceExtractionMs = 0, sourceBoundaryMs = 0, targetPreparationMs = 0;
    double mappingMs = 0, sampleMs = 0, residualBuildMs = 0, residualBoundMs = 0;
    double boundaryPositionMs = 0, boundaryAngleMs = 0; // Both tangents + joint budget.
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
               boundaryTangentU.state == CheckState::Passed &&
               boundaryTangentV.state == CheckState::Passed &&
               boundaryJointTangent.state == CheckState::Passed &&
               correspondence == CheckState::Passed;
    }
};
MotherDistanceValidationReport ValidateMotherDistanceConversion(
    const sggk::BSplineSurface &, const AnalyticPatch &,
    const MotherDistanceValidationOptions & = {});
} // namespace bspline2analytic
