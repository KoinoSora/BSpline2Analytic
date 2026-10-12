#include "ApproximateInternal.h"
namespace bspline2analytic::detail
{
RefinementResult RefinePlaneOffset(AnalyticPatch &p, const SampleCache &samples,
                                   const ConversionOptions &o, ConversionStatistics &statistics,
                                   EvaluationWorkspace &workspace)
{
    RefinementDiagnostics d;
    d.level = ParameterLevel::SizeOrOffset;
    d.activeDimension = 1;
    d.sampleVersion = samples.version;
    d.finalTrust = .1;
    unsigned evaluations = workspace.evaluations, resizes = workspace.resizes;
    auto finish = [&]() {
        d.residualEvaluations = workspace.evaluations - evaluations;
        d.workspaceResizes = workspace.resizes - resizes;
        statistics.refinementRuns.push_back(d);
        return RefinementResult{d.bestWorst < d.initialWorst, d.stopReason};
    };
    if (!samples.size())
        return finish();
    if (p.parameters.type != sggk::SurfaceType::Plane)
        throw Failure(Status::InvalidInput, "Plane offset solver requires a plane");
    auto &current = workspace.current;
    EvaluateResiduals(p, samples, o.validation, workspace, current);
    d.initialWorst = d.bestWorst = current.worst;
    d.bestCost = current.cost;
    if (!o.refinement.maxIterations)
    {
        d.stopReason = RefinementStopReason::IterationLimit;
        return finish();
    }
    if (.1 <= o.refinement.minTrustRadius)
    {
        d.stopReason = RefinementStopReason::TrustTooSmall;
        return finish();
    }
    ++d.planeOffsetSolves;
    const auto frame = workspace.context.frame;
    std::vector<double> height(samples.size()), lateral(samples.size());
    double hmin = INFINITY, hmax = -INFINITY, lower = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        V error = samples[i].point - workspace.mappedPoints[i];
        double h = Dot(error, frame.z);
        height[i] = h;
        // Avoid subtracting nearly equal squared lengths.
        lateral[i] = Norm(error - frame.z * h);
        hmin = std::min(hmin, h);
        hmax = std::max(hmax, h);
        lower = std::max(lower, lateral[i]);
    }
    double range = .1 * SampleScale(p, samples);
    double offset = std::clamp(hmin + (hmax - hmin) * .5, -range, range);
    auto distance = [&](double c) {
        double worst = 0;
        for (size_t i = 0; i < height.size(); ++i)
            worst = std::max(worst, std::hypot(lateral[i], height[i] - c));
        return worst;
    };
    double upper = distance(offset);
    lower = std::max(lower, (hmax - hmin) * .5);
    auto feasible = [&](double bound, double &c) {
        double lo = -range, hi = range;
        for (size_t i = 0; i < height.size(); ++i)
        {
            if (bound < lateral[i])
                return false;
            double width = std::sqrt(bound - lateral[i]) * std::sqrt(bound + lateral[i]);
            lo = std::max(lo, height[i] - width);
            hi = std::min(hi, height[i] + width);
            if (lo > hi)
                return false;
        }
        c = lo + (hi - lo) * .5;
        return true;
    };
    const double precision =
        o.validation.distanceTolerance * o.refinement.scalarDistanceToleranceRatio;
    // A sampled feasible midpoint goes directly to independent validation.
    // Otherwise solve the bounded convex scalar minimax problem.
    while (upper >= o.validation.distanceTolerance && upper - lower > precision &&
           d.scalarIterations < o.refinement.maxScalarIterations)
    {
        ++d.scalarIterations;
        double middle = lower + (upper - lower) * .5, c = 0;
        if (middle == lower || middle == upper)
            break;
        if (feasible(middle, c))
        {
            upper = middle;
            offset = c;
        }
        else
            lower = middle;
    }
    d.stopReason = upper >= o.validation.distanceTolerance && upper - lower > precision &&
                           d.scalarIterations == o.refinement.maxScalarIterations
                       ? RefinementStopReason::ScalarIterationLimit
                       : RefinementStopReason::PlaneOffsetSolved;
    AnalyticPatch trial = p;
    trial.parameters.origin = Point(frame.o + frame.z * offset);
    EvaluateResiduals(trial, samples, o.validation, workspace, workspace.trial);
    const auto &next = workspace.trial;
    if (next.worst < current.worst || (next.worst == current.worst && next.cost < current.cost))
    {
        p = std::move(trial);
        d.planeOffset = offset;
        d.bestWorst = next.worst;
        d.bestCost = next.cost;
        MaterializePatch(p);
    }
    return finish();
}
} // namespace bspline2analytic::detail
