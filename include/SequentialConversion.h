#pragma once
#include "MotherDistanceValidation.h"
#include "GeneralRationalBoundaryValidation.h"

namespace bspline2analytic
{
struct SequentialConversionOptions
{
    // Positive, at most 50 ms. Both branches share this one wall-clock budget.
    double maxTotalMs = 50;
    // Reserve time for fallback when iso search is inconclusive. Same start
    // instant as the total deadline; exceeding this is a branch failure only.
    double maxIsoMs = 20;
    RecognitionOptions exact;
    ConversionOptions iso;
    MotherDistanceValidationOptions mother;
    GeneralRationalBoundaryOptions fallback;
};
struct SequentialConversionResult
{
    bool accepted = false, isoAttempted = false, fallbackAttempted = false;
    bool timeBudgetExceeded = false;
    bool isoTimeBudgetExceeded = false;
    std::string path, message, isoMessage, timeoutStage;
    double totalMs = 0, isoMs = 0, fallbackMs = 0;
    unsigned isoCandidates = 0, fallbackCandidates = 0, optimizationIterations = 0;
    RecognitionResult exact;
    std::optional<AnalyticPatch> isoPatch;
    MotherDistanceValidationReport mother;
    GeneralRationalBoundaryReport fallback;
    bool Passed() const noexcept { return accepted && !timeBudgetExceeded; }
};
// Current contracts: exact rational recognition / polynomial mother-distance
// and joint U/V boundary certificate; fallback sampled interior + continuous
// rational boundary position/tangent certificate. No retired validator is used.
// The fallback gets only input geometry, never fixture truth or a supplied chart.
SequentialConversionResult ConvertIsoThenFallback(const sggk::BSplineSurface &,
                                                  const SequentialConversionOptions & = {});
}
