#pragma once
#include "GeneralHybridValidation.h"
#include <Geometry/3D/Curve/BSplineCurve3D.h>

namespace bspline2analytic
{
struct GeneralRationalBoundaryOptions
{
    GeneralConversionOptions general;
    // Optional fixed analytic box defining the rational chart of this mother.
    // No analytic PCurve is constructed. Omitted: infer a box from boundary walks.
    std::optional<std::array<double,4>> chartDomain;
    unsigned maxChartAttempts = 12, maxBoundaryNodes = 4096, maxBoundaryDepth = 14;
};
struct GeneralBoundaryPatch
{
    AnalyticParameters mother;
    // Source order Umin,Umax,Vmin,Vmax, each parameterized on [0,1].
    // Positive traversal of the wire: Vmin,Umax,-Vmax,-Umin.
    std::array<sggk::BSplineCurve3DPtr,4> boundaries;
    std::array<sggk::Point3D,4> vertices; // 00,10,01,11; shared exactly by curves.
    std::array<bool,4> tangentNotApplicable{};
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
    unsigned fusedErrorBuilds = 0, positionOnlyNodes = 0;
    unsigned tangentFastChecks = 0, tangentFastAccepted = 0, tangentFallbacks = 0;
    double interiorSampleMaxDistance = 0, interiorSampleMaxFiniteDistance = 0;
    double boundaryDistanceUpperBound = INFINITY, boundaryTangentUpperBoundDeg = INFINITY;
    // Certified output-curve distance to the ideal rational trace on the mother.
    double outputMotherDistanceUpperBound = INFINITY;
    double samplingMs = 0, recoveryMs = 0, preparationMs = 0, boundaryMs = 0, totalMs = 0;
    std::string message;
    bool Passed() const noexcept {return status==HybridValidationStatus::Passed&&patch.has_value();}
};
GeneralRationalBoundaryReport ValidateGeneralRationalBoundaries(const sggk::BSplineSurface &,
    const AnalyticParameters &,const GeneralRationalBoundaryOptions & = {});
}
