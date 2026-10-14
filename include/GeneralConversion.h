#pragma once
#include "ApproximateConversion.h"
#include <memory>
#include <cstdint>
#include <Geometry/2D/Curve2D.h>

namespace bspline2analytic
{
struct GeneralConversionOptions
{
    ConversionOptions fitting;
    // Fixed sample-group budget, including unsuccessful constructions.
    unsigned seedAttemptsPerType = 8;
    // Preferred candidates per family in each batch; remaining unique mothers
    // stay available if all candidates in this batch fail conversion.
    unsigned candidatesPerType = 3;
    double similarScoreAbsoluteRatio = .1;
    double similarScoreRelativeTolerance = .05;
    unsigned boundarySamplesPerSpan = 8;
    // Zero selects eager full scoring; otherwise use up to this many held-out
    // interior points before completing competitive candidates.
    unsigned coarseScorePointCount = 16;
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
struct GeneralTorusSeedDiagnostics
{
    // Four-normal groups actually constructed; a degenerate group still uses its attempt.
    unsigned sampleGroups = 0;
    unsigned rankDeficientGroups = 0, degenerateAxisGroups = 0, noRealAxisGroups = 0;
    // Finite axes passing the algebraic residual check, before meridian circle fitting.
    unsigned axisCount = 0, rejectedAxes = 0;
    unsigned circleFitFailures = 0, invalidRadii = 0;
    // Ring-torus candidates before geometric deduplication and held-out scoring.
    unsigned candidateCount = 0, dualCandidateGroups = 0;
};
struct GeneralClassificationDiagnostics
{
    // Generated mother candidates, before deduplication and batch selection.
    unsigned rawCandidateCount = 0;
    unsigned duplicateCandidateCount = 0;
    // Unique candidates entering scoring, including scoring failures.
    unsigned scoreEvaluationCount = 0;
    unsigned seedAttemptCount = 0;
    unsigned coarseScoreCount = 0, fullScoreCount = 0;
    unsigned completedScoreCount = 0;
    unsigned deferredCandidateCount = 0; // Incomplete after the first batch.
    unsigned backupCandidateCount = 0, batchCount = 0;
    std::uint64_t interiorPointEvaluations = 0, boundaryPointEvaluations = 0;
    GeneralTorusSeedDiagnostics torus;
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
    GeneralClassificationDiagnostics classification;
};
GeneralConversionResult ConvertGeneral(const sggk::BSplineSurface &,
                                       const GeneralConversionOptions & = {});
}
