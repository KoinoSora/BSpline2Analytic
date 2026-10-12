#include "ApproximateInternal.h"
#include "StageTimer.h"
namespace bspline2analytic
{
const char *ConversionStatusName(ConversionStatus s)
{
    switch (s)
    {
    case ConversionStatus::Success:
        return "Success";
    case ConversionStatus::InvalidInput:
        return "InvalidInput";
    case ConversionStatus::Unsupported:
        return "Unsupported";
    case ConversionStatus::NoCandidate:
        return "NoCandidate";
    case ConversionStatus::Degenerate:
        return "Degenerate";
    case ConversionStatus::OptimizationFailed:
        return "OptimizationFailed";
    case ConversionStatus::VerificationFailed:
        return "VerificationFailed";
    case ConversionStatus::Unresolved:
        return "Unresolved";
    }
    return "Unknown";
}
const char *CheckStateName(CheckState s)
{
    switch (s)
    {
    case CheckState::NotChecked:
        return "NotChecked";
    case CheckState::Passed:
        return "Passed";
    case CheckState::Exceeded:
        return "Exceeded";
    case CheckState::Unresolved:
        return "Unresolved";
    }
    return "Unknown";
}
const char *RefinementStopReasonName(RefinementStopReason reason)
{
    switch (reason)
    {
    case RefinementStopReason::NoSamples:
        return "NoSamples";
    case RefinementStopReason::MarginReached:
        return "MarginReached";
    case RefinementStopReason::SampleFeasibleStagnation:
        return "SampleFeasibleStagnation";
    case RefinementStopReason::Stagnation:
        return "Stagnation";
    case RefinementStopReason::NoAcceptedStep:
        return "NoAcceptedStep";
    case RefinementStopReason::TrustTooSmall:
        return "TrustTooSmall";
    case RefinementStopReason::SmallStep:
        return "SmallStep";
    case RefinementStopReason::IterationLimit:
        return "IterationLimit";
    case RefinementStopReason::PlaneOffsetSolved:
        return "PlaneOffsetSolved";
    case RefinementStopReason::ScalarIterationLimit:
        return "ScalarIterationLimit";
    case RefinementStopReason::SolverFailure:
        return "SolverFailure";
    }
    return "Unknown";
}
namespace detail
{
bool RefineAndValidateCandidate(const sggk::BSplineSurface &s, AnalyticPatch &candidate,
                                const ConversionOptions &o, ConversionResult &result,
                                std::shared_ptr<const SourceSampleCache> &base,
                                unsigned candidateIndex)
{
    auto validate = [&]() {
        StageTimer timer(result.timing.validation);
        result.validation = ValidateConversion(s, candidate, o.validation);
        result.statistics.validationPieces += result.validation.piecesVisited;
    };
    validate();
    if (result.validation.Passed() || !o.refinement.maxRounds)
        return false;
    SampleCache samples;
    {
        StageTimer timer(result.timing.optimization);
        if (!base)
        {
            auto initial = SampleSurface(s, o.sampling);
            auto source = std::make_shared<SourceSampleCache>();
            source->points = std::move(initial.points);
            source->boundaryOrder = std::move(initial.boundaryOrder);
            source->keys = std::move(initial.keys);
            base = std::move(source);
            ++result.statistics.baseSamplingBuilds;
        }
        samples = SampleCache(base);
    }
    auto addWitnesses = [&]() {
        StageTimer timer(result.timing.optimization);
        unsigned count = AddWitnessSamples(s, result.validation, o.sampling, samples);
        result.statistics.witnessSamplesAdded += count;
        result.statistics.sampleCount =
            std::max(result.statistics.sampleCount, static_cast<unsigned>(samples.size()));
        return count;
    };
    addWitnesses();
    EvaluationWorkspace workspace;
    ParameterLevel level = o.refinement.enableProgressiveVariables ? ParameterLevel::SizeOrOffset
                                                                   : ParameterLevel::FullPatch;
    RefinementEntryReason reason = RefinementEntryReason::Initial;
    bool progressed = false, forceAttempt = false;
    for (unsigned round = 0; round < o.refinement.maxRounds; ++round)
    {
        ++result.statistics.refinementRounds;
        unsigned remaining = o.refinement.maxIterations;
        bool retryWithWitnesses = false;
        for (;;)
        {
            unsigned before = result.statistics.optimizationIterations;
            {
                StageTimer timer(result.timing.optimization);
                RefinementResult refinement;
                if (level == ParameterLevel::SizeOrOffset &&
                    candidate.parameters.type == sggk::SurfaceType::Plane &&
                    o.refinement.enablePlaneOffsetSolver)
                    refinement =
                        RefinePlaneOffset(candidate, samples, o, result.statistics, workspace);
                else
                {
                    unsigned cap = level == ParameterLevel::SizeOrOffset    ? 6
                                   : level == ParameterLevel::MotherSurface ? 10
                                                                            : remaining;
                    refinement =
                        RefineCandidateLevel(candidate, samples, o, result.statistics, level,
                                             std::min(remaining, cap), forceAttempt, workspace);
                }
                progressed |= refinement.improved;
                auto &d = result.statistics.refinementRuns.back();
                d.candidateIndex = candidateIndex;
                d.roundIndex = round + 1;
                d.entryReason = reason;
            }
            remaining -= result.statistics.optimizationIterations - before;
            validate();
            if (result.validation.Passed())
                return progressed;
            if (!o.refinement.maxIterations || o.refinement.minTrustRadius >= .1)
                return progressed;
            unsigned added = addWitnesses();
            if (!remaining)
            {
                // A genuinely new sample set may start a new feedback round;
                // changing variable levels alone never resets this round's budget.
                retryWithWitnesses = added != 0;
                reason = RefinementEntryReason::WitnessSamples;
                forceAttempt = true;
                break;
            }
            const auto &d = result.statistics.refinementRuns.back();
            bool sampledFeasible = d.bestWorst < 1;
            bool residualCheckFailed =
                result.validation.geometry.state != CheckState::Passed ||
                result.validation.boundaryTangent.state != CheckState::Passed;
            if (sampledFeasible && residualCheckFailed && added &&
                round + 1 < o.refinement.maxRounds)
            {
                // A new witness can invalidate sampled feasibility. Retry this
                // level, with a new feedback round and its shared iteration budget.
                retryWithWitnesses = true;
                reason = RefinementEntryReason::WitnessSamples;
                forceAttempt = true;
                break;
            }
            if (level != ParameterLevel::FullPatch)
            {
                level = level == ParameterLevel::SizeOrOffset ? ParameterLevel::MotherSurface
                                                              : ParameterLevel::FullPatch;
                reason = RefinementEntryReason::MoreVariables;
                // A failed independent check must not be bypassed by the old
                // sampled margin exit. Normals remain acceptance-only.
                forceAttempt = true;
                continue;
            }
            retryWithWitnesses = added != 0;
            reason = RefinementEntryReason::WitnessSamples;
            forceAttempt = true;
            break;
        }
        if (!retryWithWitnesses)
            break;
    }
    return progressed;
}
} // namespace detail
ConversionResult ConvertApproximate(const sggk::BSplineSurface &s, const ConversionOptions &o)
{
    using namespace detail;
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    ConversionResult result;
    try
    {
        CheckApproximateInput(s, o.validation);
        if (o.sampling.interiorIsoCurves > 5 || !o.sampling.interiorPointsPerSpan ||
            o.sampling.interiorPointsPerSpan > 8 || o.sampling.boundaryPointsPerSpan < 3 ||
            o.sampling.boundaryPointsPerSpan > 32 || !o.sampling.maxSamples)
            throw Failure(Status::InvalidInput, "Invalid sampling options");
        if (!(o.refinement.candidateToleranceFactor >= 1 &&
              std::isfinite(o.refinement.candidateToleranceFactor)) ||
            o.refinement.maxIterations > 500 || o.refinement.maxRounds > 20)
            throw Failure(Status::InvalidInput, "Invalid refinement options");
        const auto &refine = o.refinement;
        auto nonnegative = [](double x) { return std::isfinite(x) && x >= 0; };
        if ((refine.jacobianMode != JacobianMode::AnalyticWithFallback &&
             refine.jacobianMode != JacobianMode::FiniteDifference) ||
            refine.maxWarmupIterations > 500 || !refine.stagnationWindow ||
            refine.stagnationWindow > 500 || !refine.maxConsecutiveFailures ||
            refine.maxConsecutiveFailures > 500 || !nonnegative(refine.worstImprovementAbsTol) ||
            !nonnegative(refine.worstImprovementRelTol) ||
            !nonnegative(refine.costImprovementRelTol) ||
            !(refine.minTrustRadius > 0 && refine.minTrustRadius <= .1) ||
            !std::isfinite(refine.minTrustRadius) ||
            !nonnegative(refine.smallModelChangeTolerance) || !refine.maxScalarIterations ||
            refine.maxScalarIterations > 500 ||
            !(refine.scalarDistanceToleranceRatio > 0 && refine.scalarDistanceToleranceRatio < 1) ||
            !std::isfinite(refine.scalarDistanceToleranceRatio))
            throw Failure(Status::InvalidInput, "Invalid refinement stopping options");
        std::vector<AnalyticPatch> candidates;
        {
            StageTimer timer(result.timing.classification);
            candidates = ClassifyApproximate(s, o);
        }
        result.statistics.candidateCount = static_cast<unsigned>(candidates.size());
        result.status = ConversionStatus::NoCandidate;
        result.message = "No analytic candidate could be initialized from line/circle families";
        unsigned candidateIndex = 0;
        std::shared_ptr<const SourceSampleCache> baseSamples;
        for (auto candidate : candidates)
        {
            ++candidateIndex;
            bool progressed =
                RefineAndValidateCandidate(s, candidate, o, result, baseSamples, candidateIndex);
            const auto &report = result.validation;
            result.patch = candidate;
            result.validation = report;
            result.message = report.message;
            if (report.Passed())
            {
                result.status = ConversionStatus::Success;
                break;
            }
            bool exceeded = report.geometry.state == CheckState::Exceeded ||
                            report.boundaryTangent.state == CheckState::Exceeded ||
                            report.boundaryNormal.state == CheckState::Exceeded;
            result.status = exceeded ? ConversionStatus::VerificationFailed
                                     : (result.timing.optimization.calls && !progressed
                                            ? ConversionStatus::OptimizationFailed
                                            : ConversionStatus::Unresolved);
        }
    }
    catch (const Failure &e)
    {
        result.message = e.what();
        result.status =
            e.status == Status::Unsupported
                ? ConversionStatus::Unsupported
                : (e.status == Status::Degenerate
                       ? ConversionStatus::Degenerate
                       : (e.status == Status::VerificationFailed ? ConversionStatus::Unresolved
                                                                 : ConversionStatus::InvalidInput));
    }
    catch (const std::exception &e)
    {
        result.status = ConversionStatus::InvalidInput;
        result.message = e.what();
    }
    result.timing.totalMs = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    result.timing.otherMs = result.timing.totalMs - result.timing.classification.elapsedMs -
                            result.timing.optimization.elapsedMs -
                            result.timing.validation.elapsedMs;
    return result;
}
} // namespace bspline2analytic
