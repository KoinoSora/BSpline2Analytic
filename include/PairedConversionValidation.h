#pragma once
#include "ApproximateConversion.h"
#include <memory>

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
    unsigned sourceCacheHits = 0, uExtractionMatrices = 0, vExtractionMatrices = 0;
    // Subsets of certificate.timing.preparation; do not add to totalMs again.
    double nurbsMs = 0, mappingMs = 0;
    double sourceExtractionMs = 0, sourceCacheMs = 0, localTransformMs = 0;
    unsigned fastBoundPasses = 0, tightBoundChecks = 0;
    // Subsets of certificate.timing.surfaceBounds, exclusive of subdivision.
    double surfaceComposeMs = 0, surfaceProductMs = 0, surfaceDifferenceMs = 0;
    double surfaceFastBoundMs = 0, surfaceTightBoundMs = 0;
    // Schema 2: source/target preparation are subsets of preparation; edge
    // composition/algebra/bounds/subdivision are subsets of boundaryChecks.
    unsigned preparationSchema = 2;
    unsigned sourceBoundaryBuilds = 0, boundarySegments = 0;
    unsigned recoveryTraceHits = 0, recoveryTraceBuilds = 0;
    double sourceCenterMs = 0, sourceBoundaryMs = 0, targetTransformMs = 0;
    double targetBoundaryMs = 0, edgeComposeMs = 0, edgeAlgebraMs = 0;
    double edgeBoundsMs = 0, edgeSubdivisionMs = 0;
    bool Passed() const noexcept
    {
        return status == PairedValidationStatus::Passed && certificate.Passed();
    }
};
// Optional cache for sequential validations of one source against changing
// candidates. Exact source-content comparison invalidates stale entries.
// Candidate geometry, domains, maps and certificates are never cached.
// Use a separate context per concurrent caller; it is not thread-safe.
class PairedValidationContext
{
  public:
    PairedValidationContext();
    ~PairedValidationContext();
    PairedValidationContext(PairedValidationContext &&) noexcept;
    PairedValidationContext &operator=(PairedValidationContext &&) noexcept;
    PairedValidationContext(const PairedValidationContext &) = delete;
    PairedValidationContext &operator=(const PairedValidationContext &) = delete;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    friend PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &,
                                                           const AnalyticPatch &,
                                                           PairedValidationContext &,
                                                           const PairedValidationOptions &);
};
PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &, const AnalyticPatch &,
                                                const PairedValidationOptions & = {});
PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &, const AnalyticPatch &,
                                                PairedValidationContext &,
                                                const PairedValidationOptions & = {});
const char *PairedValidationStatusName(PairedValidationStatus) noexcept;
} // namespace bspline2analytic
