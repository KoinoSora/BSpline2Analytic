#include "MotherDistanceInternal.h"

namespace bspline2analytic::detail::mother
{
namespace
{
bool Samples(paired::AnalysisContext &context, const AnalyticParameters &p,
             const RecognitionOptions &options, MotherDistanceValidationReport &out)
{
    Measure time(out.sampleMs);
    bool passed = true;
    auto check = [&](const SurfaceSample &sample) {
        V point = sample.point;
        ++out.sampleCount;
        double d = Distance(p, Point(point));
        if (!std::isfinite(d)) throw Failure(Status::VerificationFailed, "Nonfinite mother sample");
        double guard = 512 * std::numeric_limits<double>::epsilon() *
            (Norm(point) + Norm(Get(p.origin)) + std::abs(p.radius) + std::abs(p.majorRadius));
        out.motherDistance.maxSampledUpperBound = std::max(out.motherDistance.maxSampledUpperBound, d + guard);
        if (d - guard >= options.distanceTolerance)
        {
            passed = false;
            out.motherDistance.parameterRegion = {sample.u, sample.u, sample.v, sample.v};
        }
    };
    for (const auto &sample : context.samples.Greville(false)) check(sample);
    const auto &s = context.samples.source;
    for (bool along : {true, false})
        for (double t : SampleParameters(along ? s.KnotsU() : s.KnotsV(),
                 along ? s.MinParamU() : s.MinParamV(), along ? s.MaxParamU() : s.MaxParamV(),
                 options.samplesPerSpan))
            for (bool last : {false, true})
                check(context.samples.At(along ? t : (last ? s.MaxParamU() : s.MinParamU()),
                     along ? (last ? s.MaxParamV() : s.MinParamV()) : t,
                     sggk::DerivType::NONE, sggk::DerivType::NONE, false));
    if (!passed)
    {
        out.motherDistance.state = CheckState::Exceeded;
        out.motherDistance.hasWitness = true;
    }
    return passed;
}
bool Run(paired::AnalysisContext &context, const AnalyticPatch &patch,
         const RecoveredCandidate *recovered, const MotherDistanceValidationOptions &options,
         MotherDistanceValidationReport &out, PairedValidationReport &aux, ValidationTimer &timer)
{
    const auto &s = context.samples.source;
    auto boundary = options.boundary;
    boundary.validation.distanceTolerance = options.position.distanceTolerance;
    ValidateOptions(options.position);
    CheckApproximateInput(s, boundary.validation);
    paired::CheckPatch(patch);
    if (!IntervalEnvironmentSupported())
        throw Failure(Status::VerificationFailed, "Unsupported interval environment");
    if (!options.position.requireBound)
        throw Failure(Status::InvalidInput, "Mother validator requires a continuous bound");
    if (s.DegreeU() > 5 || s.DegreeV() > 5)
        throw Failure(Status::Unsupported, "Mother validator supports polynomial degree <=5");
    if (!boundary.maxMappingSegments || boundary.maxMappingSegments > 4096 ||
        boundary.maxMappingDepth > 20 || !boundary.maxCommonCells || boundary.maxCommonCells > 1000000)
        throw Failure(Status::InvalidInput, "Invalid correspondence budget");
    if (s.KnotsU().size() - 1 > boundary.maxCommonCells / (s.KnotsV().size() - 1))
        throw paired::Budget("Source cell budget exhausted");
    if (context.source.Matches(s)) ++aux.sourceCacheHits;
    else
    {
        auto prepared = paired::PrepareSource(s, aux);
        ++aux.sourceExtractions;
        context.source.Store(s, std::move(prepared));
    }
    const auto &source = context.source.prepared;
    if (source.cells.empty()) throw Failure(Status::InvalidInput, "Empty source domain");
    {
        ValidationTimer::Scope position(timer, out.timing.surfaceBounds);
        if (!Samples(context, patch.parameters, options.position, out)) return false;
    }
    if (!ValidatePosition(source, patch.parameters, options.position, out, timer)) return false;
    paired::Products products;
    paired::Chart chart;
    paired::Boundaries target;
    {
        Measure time(out.targetPreparationMs);
        chart = paired::BuildChart(patch);
        paired::PlaceChart(chart, patch.parameters, source.origin);
        target = paired::PrepareTargetBoundaries(chart, products);
    }
    aux.uMap = {{s.MinParamU(), s.MaxParamU()}, {0, 1}};
    aux.vMap = {{s.MinParamV(), s.MaxParamV()}, {0, 1}};
    for (unsigned attempt = 0; attempt < (boundary.enablePiecewiseMapping ? 2u : 1u); ++attempt)
    {
        // A failed refinement can partially change a map. Never attach a bound
        // from the previous correspondence to the maps returned after failure.
        out.boundaryPosition = {}; out.boundaryTangent = {}; out.boundaryNormal = {};
        ++out.mappingAttempts;
        if (attempt) paired::RefineBoundaryMaps(context, patch, chart, boundary, aux, recovered);
        paired::CompleteChartMaps(aux, chart);
        bool positionOK;
        paired::BoundaryResult angular;
        {
            ValidationTimer::Scope edges(timer, out.timing.boundaryChecks);
            positionOK = ValidateEdgePositions(source, chart, boundary, aux, products, out);
            if (positionOK)
            {
                Measure time(out.boundaryAngleMs);
                angular = paired::ValidateBoundaries(source.boundaries, target, aux.uMap,
                    aux.vMap, boundary, aux, products);
            }
        }
        {
            ValidationTimer::Scope reporting(timer, out.timing.reporting);
            if (positionOK && angular.passed)
            {
                out.boundaryTangent.upperBound = (Atan(I(angular.tangentRatio)) * I(180) / IPi()).hi;
                out.boundaryNormal.upperBound = (Atan(I(angular.normalRatio)) * I(180) / IPi()).hi;
                if (out.boundaryTangent.upperBound < boundary.validation.boundaryTangentToleranceDeg &&
                    out.boundaryNormal.upperBound < boundary.validation.boundaryNormalToleranceDeg)
                {
                    out.boundaryTangent.state = out.boundaryNormal.state = out.correspondence = CheckState::Passed;
                    return true;
                }
            }
            out.boundaryTangent.state = angular.tangentExceeded ? CheckState::Exceeded : CheckState::Unresolved;
            out.boundaryNormal.state = angular.normalExceeded ? CheckState::Exceeded : CheckState::Unresolved;
            out.boundaryTangent.upperBound = out.boundaryNormal.upperBound = INFINITY;
        }
    }
    return false;
}
MotherDistanceValidationReport Validate(paired::AnalysisContext &context, const AnalyticPatch &patch,
    const RecoveredCandidate *recovered, const MotherDistanceValidationOptions &options)
{
    MotherDistanceValidationReport out;
    PairedValidationReport aux;
    ValidationTimer timer(out.timing);
    try
    {
        bool passed = Run(context, patch, recovered, options, out, aux, timer);
        out.status = passed ? PairedValidationStatus::Passed : PairedValidationStatus::NotCertified;
        out.message = passed ? "Floating-point guarded Bernstein mother-distance bound, with interval "
            "paired boundary position/tangent/normal checks; not an interior paired or Hausdorff certificate"
            : "Mother distance or boundary correspondence not certified; inspect individual check states";
    }
    catch (const paired::Budget &e)
    {
        out.status = PairedValidationStatus::BudgetExceeded; out.message = e.what();
    }
    catch (const Failure &e)
    {
        out.status = e.status == Status::InvalidInput ? PairedValidationStatus::InvalidInput :
            e.status == Status::Unsupported || e.status == Status::Degenerate ? PairedValidationStatus::Unsupported :
            PairedValidationStatus::NotCertified;
        out.message = e.what();
    }
    catch (const std::exception &e)
    {
        out.status = PairedValidationStatus::NumericalFailure; out.message = e.what();
    }
    timer.Switch(out.timing.reporting);
    out.sourceExtractions = aux.sourceExtractions; out.sourceCacheHits = aux.sourceCacheHits;
    out.sourceExtractionMs = aux.sourceExtractionMs; out.sourceBoundaryMs = aux.sourceBoundaryMs;
    out.mappingMs = aux.mappingMs; out.boundaryNodes = aux.boundaryNodes;
    out.uMap = std::move(aux.uMap); out.vMap = std::move(aux.vMap);
    if (out.status != PairedValidationStatus::Passed) out.correspondence = CheckState::Unresolved;
    // Only completed individual checks retain finite bounds after an early exit.
    for (auto *check : {&out.motherDistance, &out.boundaryPosition, &out.boundaryTangent, &out.boundaryNormal})
        if (check->state != CheckState::Passed) check->upperBound = INFINITY;
    timer.Finish();
    return out;
}
} // namespace
MotherDistanceValidationReport ValidatePrepared(paired::AnalysisContext &context,
    const RecoveredCandidate &candidate, const MotherDistanceValidationOptions &options)
{
    return Validate(context, candidate.patch, &candidate, options);
}
} // namespace bspline2analytic::detail::mother

namespace bspline2analytic
{
MotherDistanceValidationReport ValidateMotherDistanceConversion(const sggk::BSplineSurface &s,
    const AnalyticPatch &p, const MotherDistanceValidationOptions &options)
{
    detail::paired::AnalysisContext context(s);
    return detail::mother::Validate(context, p, nullptr, options);
}
} // namespace bspline2analytic
