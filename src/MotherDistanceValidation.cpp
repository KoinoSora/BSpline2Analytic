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
        if (!std::isfinite(d))
            throw Failure(Status::VerificationFailed, "Nonfinite mother sample");
        double guard =
            512 * std::numeric_limits<double>::epsilon() *
            (Norm(point) + Norm(Get(p.origin)) + std::abs(p.radius) + std::abs(p.majorRadius));
        out.motherDistance.maxSampledUpperBound =
            std::max(out.motherDistance.maxSampledUpperBound, d + guard);
        if (d - guard >= options.distanceTolerance)
        {
            passed = false;
            out.motherDistance.parameterRegion = {sample.u, sample.u, sample.v, sample.v};
        }
    };
    for (const auto &sample : context.samples.Greville(false))
        check(sample);
    const auto &s = context.samples.source;
    for (bool along : {true, false})
        for (double t : SampleParameters(
                 along ? s.KnotsU() : s.KnotsV(), along ? s.MinParamU() : s.MinParamV(),
                 along ? s.MaxParamU() : s.MaxParamV(), options.samplesPerSpan))
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
    if (!(options.jointTangentToleranceDeg > 0 && options.jointTangentToleranceDeg < 90))
        throw Failure(Status::InvalidInput, "Joint tangent tolerance must be finite and in (0,90) degrees");
    // CheckApproximateInput validates this legacy field. Only the independent
    // joint budget controls the third validator's actual angle acceptance.
    boundary.validation.boundaryTangentToleranceDeg = options.jointTangentToleranceDeg;
    ValidateOptions(options.position);
    CheckApproximateInput(s, boundary.validation, false);
    paired::CheckPatch(patch, false);
    if (!IntervalEnvironmentSupported())
        throw Failure(Status::VerificationFailed, "Unsupported interval environment");
    if (!options.position.requireBound)
        throw Failure(Status::InvalidInput, "Mother validator requires a continuous bound");
    if (s.DegreeU() > 5 || s.DegreeV() > 5)
        throw Failure(Status::Unsupported, "Mother validator supports polynomial degree <=5");
    if (!boundary.maxMappingSegments || boundary.maxMappingSegments > 4096 ||
        boundary.maxMappingDepth > 20 || !boundary.maxCommonCells ||
        boundary.maxCommonCells > 1000000)
        throw Failure(Status::InvalidInput, "Invalid correspondence budget");
    if (s.KnotsU().size() - 1 > boundary.maxCommonCells / (s.KnotsV().size() - 1))
        throw paired::Budget("Source cell budget exhausted");
    if (!context.motherSource)
        context.motherSource = std::make_shared<SourceCache>();
    if (context.motherSource->Matches(s))
        ++out.sourceCacheHits;
    else
    {
        Measure time(out.sourceExtractionMs);
        context.motherSource->Initialize(s, out);
        ++out.sourceExtractions;
    }
    auto &source = *context.motherSource;
    out.sourceCells = static_cast<unsigned>(source.cells.size());
    if (source.cells.empty())
        throw Failure(Status::InvalidInput, "Empty source domain");
    {
        ValidationTimer::Scope position(timer, out.timing.surfaceBounds);
        if (!Samples(context, patch.parameters, options.position, out))
            return false;
    }
    if (!ValidatePosition(source, patch.parameters, options.position, out, timer, context.samples.deadline))
        return false;
    {
        Measure time(out.sourceBoundaryMs);
        source.PrepareBoundary(out);
    }
    paired::Products products;
    paired::Chart chart;
    BoundaryTangents target;
    {
        Measure time(out.targetPreparationMs);
        chart = paired::BuildChart(patch, false);
        paired::PlaceChart(chart, patch.parameters, source.origin);
        target = PrepareTargetTangents(chart, products);
    }
    aux.uMap = {{s.MinParamU(), s.MaxParamU()}, {0, 1}};
    aux.vMap = {{s.MinParamV(), s.MaxParamV()}, {0, 1}};
    bool boundaryTightened = false;
    for (unsigned attempt = 0; attempt < (boundary.enablePiecewiseMapping ? 2u : 1u); ++attempt)
    {
        // A failed refinement can partially change a map. Never attach a bound
        // from the previous correspondence to the maps returned after failure.
        out.boundaryPosition = {};
        out.boundaryTangentU = {};
        out.boundaryTangentV = {};
        out.boundaryJointTangent = {};
        out.jointTangentSineSquaredUpperBound = INFINITY;
        out.jointTangentFailure.clear();
        out.message.clear();
        ++out.mappingAttempts;
        if (attempt)
            paired::RefineBoundaryMaps(context, patch, chart, boundary, aux, recovered);
        paired::CompleteChartMaps(aux, chart);
        bool positionOK;
        bool angularOK = false;
        for (unsigned retry = 0; retry < 2; ++retry)
        {
            {
                ValidationTimer::Scope edges(timer, out.timing.boundaryChecks);
                positionOK =
                    ValidateEdgePositions(*source.boundary, chart, boundary, aux, products, out, context.samples.deadline);
                if (positionOK)
                {
                    Measure time(out.boundaryAngleMs);
                    angularOK = ValidateJointTangents(source.boundaryTangents, target, options,
                                                      aux, products, out, context.samples.deadline);
                }
            }
            if ((positionOK && angularOK) || boundaryTightened ||
                out.boundaryPosition.state == CheckState::Exceeded ||
                out.boundaryJointTangent.state == CheckState::Exceeded)
                break;
            // Only unresolved interval checks warrant a tighter local source.
            // Keep the same correspondence and share the existing node budget.
            Measure time(out.sourceBoundaryMs);
            source.PrepareBoundary(out, true);
            boundaryTightened = true;
            ++out.boundaryFallbacks;
            out.boundaryPosition = {};
            angularOK = false;
            out.boundaryTangentU = {};
            out.boundaryTangentV = {};
            out.boundaryJointTangent = {};
            out.jointTangentSineSquaredUpperBound = INFINITY;
            out.jointTangentFailure.clear();
            out.message.clear();
        }
        {
            ValidationTimer::Scope reporting(timer, out.timing.reporting);
            if (positionOK && angularOK)
            {
                out.correspondence = CheckState::Passed;
                return true;
            }
        }
    }
    return false;
}
MotherDistanceValidationReport Validate(paired::AnalysisContext &context,
                                        const AnalyticPatch &patch,
                                        const RecoveredCandidate *recovered,
                                        const MotherDistanceValidationOptions &options)
{
    MotherDistanceValidationReport out;
    PairedValidationReport aux;
    ValidationTimer timer(out.timing);
    try
    {
        bool passed = Run(context, patch, recovered, options, out, aux, timer);
        out.status = passed ? PairedValidationStatus::Passed : PairedValidationStatus::NotCertified;
        if (passed)
            out.message = "Floating-point guarded Bernstein mother-distance bound, with interval "
                "paired boundary positions and joint U/V tangent budget; no direct normal "
                "certification or fallback, no adjacent-face G1 check, not an interior paired "
                "or Hausdorff certificate";
        else if (out.message.empty())
            out.message = "Mother distance or boundary correspondence not certified; inspect "
                          "individual check states";
    }
    catch (const paired::Budget &e)
    {
        out.status = PairedValidationStatus::BudgetExceeded;
        out.message = e.what();
    }
    catch (const Failure &e)
    {
        out.status = e.status == Status::InvalidInput ? PairedValidationStatus::InvalidInput
                     : e.status == Status::Unsupported || e.status == Status::Degenerate
                         ? PairedValidationStatus::Unsupported
                         : PairedValidationStatus::NotCertified;
        out.message = e.what();
    }
    catch (const DeadlineExceeded &) { throw; }
    catch (const std::exception &e)
    {
        out.status = PairedValidationStatus::NumericalFailure;
        out.message = e.what();
    }
    timer.Switch(out.timing.reporting);
    out.mappingMs = aux.mappingMs;
    out.boundaryNodes = aux.boundaryNodes;
    out.uMap = std::move(aux.uMap);
    out.vMap = std::move(aux.vMap);
    if (out.status != PairedValidationStatus::Passed)
        out.correspondence = CheckState::Unresolved;
    // Only completed individual checks retain finite bounds after an early exit.
    for (auto *check : {&out.motherDistance, &out.boundaryPosition, &out.boundaryTangentU,
                        &out.boundaryTangentV, &out.boundaryJointTangent})
        if (check->state != CheckState::Passed)
            check->upperBound = INFINITY;
    if (out.boundaryJointTangent.state != CheckState::Passed)
        out.jointTangentSineSquaredUpperBound = INFINITY;
    if (!out.jointTangentFailure.empty() &&
        out.message.find(out.jointTangentFailure) == std::string::npos)
        out.message = out.jointTangentFailure + ": " + out.message +
                      "; joint condition was not certified, no direct normal fallback";
    timer.Finish();
    return out;
}
} // namespace
MotherDistanceValidationReport ValidatePrepared(paired::AnalysisContext &context,
                                                const RecoveredCandidate &candidate,
                                                const MotherDistanceValidationOptions &options)
{
    return Validate(context, candidate.patch, &candidate, options);
}
} // namespace bspline2analytic::detail::mother

namespace bspline2analytic
{
MotherDistanceValidationReport ValidateMotherDistanceConversion(
    const sggk::BSplineSurface &s, const AnalyticPatch &p,
    const MotherDistanceValidationOptions &options)
{
    detail::paired::AnalysisContext context(s);
    return detail::mother::Validate(context, p, nullptr, options);
}
} // namespace bspline2analytic
