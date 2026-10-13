#pragma once
#include "ApproximateConversion.h"
#include <memory>
#include <Geometry/2D/Curve2D.h>

namespace bspline2analytic
{
struct GeneralConversionOptions
{
    ConversionOptions fitting;
    unsigned seedAttemptsPerType = 64;
    unsigned candidatesPerType = 3;
    double similarScoreAbsoluteRatio = .1;
    double similarScoreRelativeTolerance = .05;
    unsigned boundarySamplesPerSpan = 8;
};
enum class PCurveFailure
{
    None, NonlinearPCurve, PCurveOutOfRange, NonMonotonePCurve,
    InconsistentPeriodicBranch, UnresolvedDegeneratePCurve,
    InvalidBoundaryLoop, BoundaryGeometryExceeded, ValidationUnresolved, ChartRecoveryUnresolved
};
const char *PCurveFailureName(PCurveFailure);
struct LinearPCurve
{
    // Canonical mother UV, with periods unfolded. Source edge order is retained:
    // Umin/Umax vary with source v; Vmin/Vmax vary with source u.
    sggk::Point2D start, end;
    sggk::Point3D projectedStart, projectedEnd;
    bool collapsed3D = false;
    // A degree-one bounded curve in the SDK motherSurface's own UV convention.
    sggk::Curve2DPtr sdkCurve;
    sggk::Point2D CalcPoint(double fraction) const;
};
// Finite simple UV quadrilateral (a collapsed side may give a triangle) on the
// analytic mother, bounded by four linear PCurves. It owns no source spline.
// motherSurface alone does not include these trims. Canonical UV and SDK UV
// are deliberately exposed separately, as in the existing AnalyticPatch.
class GeneralAnalyticPatch
{
public:
    const AnalyticParameters &Parameters() const noexcept { return parameters; }
    const sggk::SurfacePtr &MotherSurface() const noexcept { return motherSurface; }
    const std::array<double, 4> &SourceDomain() const noexcept { return domain; }
    const std::array<LinearPCurve, 4> &PCurves() const noexcept { return pcurves; }
    const std::array<sggk::Point2D, 4> &Corners() const noexcept { return corners; }
    // Canonical analytic UV, not the source spline parameters.
    sggk::Point3D CalcPoint(double u, double v) const;
    sggk::Point3D CalcBoundaryPoint(BoundaryId edge, double fraction) const;
    GeneralAnalyticPatch(const AnalyticParameters &, const std::array<double, 4> &,
                        const std::array<sggk::Point2D, 4> &,
                        const std::array<LinearPCurve, 4> &);
private:
    AnalyticParameters parameters;
    sggk::SurfacePtr motherSurface;
    std::array<double, 4> domain;
    std::array<sggk::Point2D, 4> corners;
    std::array<LinearPCurve, 4> pcurves;
};
struct GeneralCandidateScore
{
    sggk::SurfaceType type = sggk::SurfaceType::Unknown;
    unsigned degreesOfFreedom = 0;
    double interiorMax = 0, interiorRms = 0, boundaryMax = 0, normalMax = 0;
};
struct GeneralValidationReport
{
    CheckState geometry = CheckState::NotChecked;
    CheckState boundaryHausdorff = CheckState::NotChecked;
    CheckState pcurveLinearity = CheckState::NotChecked;
    PCurveFailure pcurveFailure = PCurveFailure::None;
    BoundaryId failedBoundary = BoundaryId::Umin;
    double maxSampledPCurveDeviation = 0;
    double distanceUpperBound = std::numeric_limits<double>::infinity();
    double boundaryHausdorffUpperBound = std::numeric_limits<double>::infinity();
    unsigned piecesVisited = 0;
    unsigned chartAttempts = 0;
    unsigned chartIterations = 0;
    bool hasWitness = false;
    std::array<double, 2> witness{};
    std::string message;
    bool Passed() const noexcept
    {
        return geometry == CheckState::Passed && boundaryHausdorff == CheckState::Passed &&
               pcurveLinearity == CheckState::Passed;
    }
};
struct GeneralConversionResult
{
    ConversionStatus status = ConversionStatus::NoCandidate;
    // Populated only after independent continuous validation succeeds.
    std::optional<GeneralAnalyticPatch> patch;
    GeneralValidationReport validation;
    std::vector<GeneralCandidateScore> candidates;
    struct Attempt
    {
        sggk::SurfaceType type=sggk::SurfaceType::Unknown;
        PCurveFailure pcurveFailure=PCurveFailure::None;
        BoundaryId failedBoundary=BoundaryId::Umin;
        bool succeeded=false;
        std::string message;
    };
    std::vector<Attempt> attempts;
    ConversionStatistics statistics;
    ConversionTiming timing;
    unsigned chartAttempts = 0;
    unsigned chartIterations = 0;
    std::string message;
};
GeneralConversionResult ConvertGeneral(const sggk::BSplineSurface &,
                                       const GeneralConversionOptions & = {});
}
