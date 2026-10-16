#pragma once
#include "GeneralHybridValidation.h"
#include <Geometry/3D/Curve/BSplineCurve3D.h>
#include <vector>

namespace bspline2analytic
{
struct GeneralRationalBoundaryOptions
{
    GeneralConversionOptions general;
    // Optional fixed analytic box defining the rational chart of this mother.
    // No analytic PCurve is constructed. Omitted: infer a box from boundary walks.
    std::optional<std::array<double,4>> chartDomain;
    unsigned maxChartAttempts = 12, maxBoundaryNodes = 4096, maxBoundaryDepth = 14;
    // Wall-clock budget; 0 explicitly disables it for diagnostics. Includes all
    // validation stages. Classified callers share one deadline across candidates.
    double maxTotalMs = 50;
    // Merge only pieces whose SOURCE and TARGET spatial diameter bounds are
    // both <= this value. Position still needs distanceTolerance. Their tangent
    // is explicitly exempted, never reported as an angle certificate. Allowed
    // range [0,1e-8]; 0 disables this policy.
    double microSegmentDiameter = 1e-8;
};
struct GeneralBoundaryTangentExemption
{
    unsigned edge = 0;
    // Source uses the original varying surface parameter; target uses [0,1].
    double sourceFirst = 0, sourceLast = 0, targetFirst = 0, targetLast = 0;
    double sourceDiameterUpper = 0, targetDiameterUpper = 0, distanceUpper = 0;
};
struct GeneralBoundaryPatch
{
    AnalyticParameters mother;
    // Source order Umin,Umax,Vmin,Vmax, each parameterized on [0,1].
    // Positive traversal of the wire: Vmin,Umax,-Vmax,-Umin.
    std::array<sggk::BSplineCurve3DPtr,4> boundaries;
    std::array<sggk::Point3D,4> vertices; // 00,10,01,11; shared exactly by curves.
    std::array<bool,4> tangentNotApplicable{};
    std::vector<GeneralBoundaryTangentExemption> tangentExemptions;
    // Region selector and construction metadata, NOT analytic linear PCurves.
    std::array<double,4> chartDomain{},sourceDomain{};
    std::array<sggk::Point2D,4> rationalCorners;
    sggk::Point3D interiorReference;
};
struct GeneralRationalBoundaryReport
{
    static constexpr bool certifiesInteriorContinuum = false;
    HybridValidationStatus status = HybridValidationStatus::NotCertified;
    std::optional<GeneralBoundaryPatch> patch;
    unsigned samples = 0, chartAttempts = 0, boundaryNodes = 0, boundarySegments = 0;
    unsigned mappingRefinements = 0, targetSegments = 0, maxTargetDegree = 0;
    unsigned microMergedSegments = 0, recovery3DRefinements = 0;
    unsigned fusedErrorBuilds = 0, positionOnlyNodes = 0;
    unsigned tangentFastChecks = 0, tangentFastAccepted = 0, tangentFallbacks = 0;
    // Positive-scale residual attempts (including interval recomputations).
    unsigned tangentResidualChecks = 0, tangentResidualAccepted = 0;
    // Includes direct interval attempts when the double envelope is disabled.
    unsigned tangentResidualIntervalFallbacks = 0;
    unsigned envelopeChecks = 0, envelopeAccepted = 0, envelopeFallbacks = 0;
    unsigned envelopePositionWitnesses = 0;
    unsigned envelopeTangentWitnesses = 0;
    double interiorSampleMaxDistance = 0, interiorSampleMaxFiniteDistance = 0;
    // Tangent bound covers noncollapsed pieces OUTSIDE patch.tangentExemptions.
    double boundaryDistanceUpperBound = INFINITY, boundaryTangentUpperBoundDeg = INFINITY;
    // Certified output-curve distance to the ideal rational trace on the mother.
    double outputMotherDistanceUpperBound = INFINITY;
    double samplingMs = 0, recoveryMs = 0, preparationMs = 0, boundaryMs = 0, totalMs = 0;
    bool timeBudgetExceeded = false;
    std::string timeoutStage;
    std::string message;
    bool Passed() const noexcept {return status==HybridValidationStatus::Passed&&patch.has_value();}
};
// Sampled interior distance; continuous exported-boundary position and tangent
// outside explicitly documented micro-segment exemptions.
// No surface normal coefficients are built and no normal angles are checked.
GeneralRationalBoundaryReport ValidateGeneralRationalBoundaries(const sggk::BSplineSurface &,
    const AnalyticParameters &,const GeneralRationalBoundaryOptions & = {});
}
