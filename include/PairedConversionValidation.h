#pragma once
#include "ApproximateConversion.h"

namespace bspline2analytic
{
// Experimental, standalone certificate. Never invokes the projection validator
// or an optimizer. Non-success does NOT prove geometric infeasibility.
enum class PairedValidationStatus
{
    Passed,
    NotCertified,
    Unsupported,
    InvalidInput,
    BudgetExceeded,
    NumericalFailure
};
struct PairedValidationOptions
{
    ValidationOptions validation;
    bool enablePiecewiseMapping = true;
    unsigned maxMappingSegments = 512;
    unsigned maxMappingDepth = 12;
    unsigned maxCommonCells = 16384;
};
struct PairedParameterMap
{
    // Target coordinates refer to the oriented, normalized rational chart,
    // NOT to linearly interpolated analytic angles.
    std::vector<double> source, target;
};
struct PairedValidationReport
{
    PairedValidationStatus status = PairedValidationStatus::NotCertified;
    // Exceeded here refers only to the reported parameter correspondence.
    ConversionValidationReport certificate;
    PairedParameterMap uMap, vMap;
    unsigned mappingAttempts = 0, mappingPointEvaluations = 0;
    unsigned commonCells = 0, surfaceNodes = 0, boundaryNodes = 0;
    unsigned sourceExtractions = 0, targetPatches = 0;
    // Subsets of certificate.timing.preparation; do not add to totalMs again.
    double nurbsMs = 0, mappingMs = 0;
    bool Passed() const noexcept
    {
        return status == PairedValidationStatus::Passed && certificate.Passed();
    }
};
PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &, const AnalyticPatch &,
                                                const PairedValidationOptions & = {});
const char *PairedValidationStatusName(PairedValidationStatus) noexcept;
} // namespace bspline2analytic
