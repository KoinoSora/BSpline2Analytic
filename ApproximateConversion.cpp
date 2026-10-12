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
    case RefinementStopReason::SolverFailure:
        return "SolverFailure";
    }
    return "Unknown";
}
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
        if (refine.maxWarmupIterations > 500 || !refine.stagnationWindow ||
            refine.stagnationWindow > 500 || !refine.maxConsecutiveFailures ||
            refine.maxConsecutiveFailures > 500 || !nonnegative(refine.worstImprovementAbsTol) ||
            !nonnegative(refine.worstImprovementRelTol) ||
            !nonnegative(refine.costImprovementRelTol) ||
            !(refine.minTrustRadius > 0 && refine.minTrustRadius <= .1) ||
            !std::isfinite(refine.minTrustRadius) || !nonnegative(refine.smallModelChangeTolerance))
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
        for (auto candidate : candidates)
        {
            ++candidateIndex;
            auto validate = [&]() {
                StageTimer timer(result.timing.validation);
                auto r = ValidateConversion(s, candidate, o.validation);
                result.statistics.validationPieces += r.piecesVisited;
                return r;
            };
            auto report = validate();
            SampleCache samples;
            bool progressed = false;
            for (unsigned round = 0; !report.Passed() && round < o.refinement.maxRounds; ++round)
            {
                {
                    StageTimer timer(result.timing.optimization);
                    if (samples.points.empty())
                        samples = SampleSurface(s, o.sampling);
                    unsigned added = AddWitnessSamples(s, report, o.sampling, samples);
                    // A stalled solve must not restart with identical data and
                    // parameters.
                    if (round > 0 && !added)
                        break;
                    ++result.statistics.refinementRounds;
                    result.statistics.sampleCount =
                        std::max(result.statistics.sampleCount,
                                 static_cast<unsigned>(samples.points.size()));
                    auto refinement = RefineCandidate(s, candidate, samples, o, result.statistics);
                    auto &diagnostic = result.statistics.refinementRuns.back();
                    diagnostic.candidateIndex = candidateIndex;
                    diagnostic.roundIndex = round + 1;
                    progressed |= refinement.improved;
                }
                report = validate();
            }
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
