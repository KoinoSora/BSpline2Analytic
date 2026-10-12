#include "ApproximateCases.h"
#include "ApproximateInternal.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
static void Require(bool ok, const char *message)
{
    if (!ok)
        throw std::runtime_error(message);
}
static void Check(const ConversionResult &r)
{
    Require(r.status == ConversionStatus::Success && r.validation.Passed(),
            "Plateau conversion did not pass independent validation");
    Require(r.patch && r.patch->parameters.type == sggk::SurfaceType::Plane,
            "Plateau did not recover a plane");
    Require(r.validation.geometry.upperBound < 1e-5 &&
                r.validation.boundaryTangent.upperBound < .1 &&
                r.validation.boundaryNormal.upperBound < .1,
            "Plateau tolerance exceeded");
    Require(r.statistics.optimizationIterations > 0 && r.timing.validation.calls >= 2,
            "Plateau did not exercise optimization and revalidation");
    Require(r.statistics.optimizationIterations < 30, "Plateau still exhausted 30 iterations");
    Require(r.statistics.refinementRuns.size() == 1 && r.statistics.refinementRounds == 1,
            "Plateau unexpectedly restarted the optimizer");
    const auto &d = r.statistics.refinementRuns.front();
    Require(d.warmupIterations + d.constrainedIterations == r.statistics.optimizationIterations &&
                d.warmupIterations <= 5 && d.acceptedSteps <= r.statistics.optimizationIterations,
            "Refinement diagnostics do not account for attempts");
    Require(d.bestWorst < 1 && d.bestWorst < d.initialWorst &&
                d.stopReason != RefinementStopReason::IterationLimit,
            "Plateau did not stop with a feasible sampled candidate");
    const auto &t = r.timing;
    Require(std::abs(t.classification.elapsedMs + t.optimization.elapsedMs +
                     t.validation.elapsedMs + t.otherMs - t.totalMs) < 1e-6,
            "Stage timing accounting mismatch");
}
static int Run(const char *path)
{
    auto references = tests::MakeCases();
    auto surface = tests::PlateauFixture(*references[0].spline);
    auto exact = Recognize(*references[0].spline);
    auto feasible = RecoverApproxPatch(*surface, exact.parameters);
    feasible.parameters.origin =
        Point(Get(feasible.parameters.origin) + Get(exact.parameters.axis) * 8e-6);
    MakeNative(feasible);
    Require(ValidateConversion(*surface, feasible).Passed(), "Known feasible plane rejected");
    Check(ConvertApproximate(*surface));
    std::ofstream csv(path);
    Require(bool(csv), "Cannot open plateau CSV");
    csv << "run,classify_ms,optimize_ms,validate_ms,total_ms,iterations,rounds,"
           "samples,distance_bound,tangent_deg_bound,normal_deg_bound,warmup_"
           "iterations,constrained_iterations,accepted_steps,jacobians,line_"
           "search_evaluations,best_worst,final_trust,stop_reason\n";
    csv << std::setprecision(17);
    std::array<double, 4> sums{};
    for (unsigned run = 0; run < 5; ++run)
    {
        auto r = ConvertApproximate(*surface);
        Check(r);
        const auto &t = r.timing;
        sums[0] += t.classification.elapsedMs;
        sums[1] += t.optimization.elapsedMs;
        sums[2] += t.validation.elapsedMs;
        sums[3] += t.totalMs;
        csv << run << ',' << t.classification.elapsedMs << ',' << t.optimization.elapsedMs << ','
            << t.validation.elapsedMs << ',' << t.totalMs << ','
            << r.statistics.optimizationIterations << ',' << r.statistics.refinementRounds << ','
            << r.statistics.sampleCount << ',' << r.validation.geometry.upperBound << ','
            << r.validation.boundaryTangent.upperBound << ','
            << r.validation.boundaryNormal.upperBound;
        const auto &d = r.statistics.refinementRuns.front();
        csv << ',' << d.warmupIterations << ',' << d.constrainedIterations << ',' << d.acceptedSteps
            << ',' << d.jacobianEvaluations << ',' << d.lineSearchEvaluations << ',' << d.bestWorst
            << ',' << d.finalTrust << ',' << RefinementStopReasonName(d.stopReason) << '\n';
        std::cout << "plateau run=" << run << " iterations=" << r.statistics.optimizationIterations
                  << " distance=" << std::setprecision(10) << r.validation.geometry.upperBound
                  << " tangent=" << r.validation.boundaryTangent.upperBound
                  << " normal=" << r.validation.boundaryNormal.upperBound
                  << " warmup=" << d.warmupIterations << " constrained=" << d.constrainedIterations
                  << " reason=" << RefinementStopReasonName(d.stopReason) << '\n';
    }
    Require(bool(csv), "Cannot write plateau CSV");
    // A separate diagnostic run does not enter the timed averages above.
    ConversionOptions traced;
    traced.refinement.recordIterationTrace = true;
    auto detail = ConvertApproximate(*surface, traced);
    Check(detail);
    std::ofstream trace(std::string(path) + ".trace.csv");
    Require(bool(trace), "Cannot open iteration trace");
    trace << "iteration,phase,worst,cost,best_worst,target,trust_before,trust_"
             "after,step_norm,model_change,step_solved,accepted\n";
    trace << std::setprecision(17);
    const auto &d = detail.statistics.refinementRuns.front();
    Require(d.trace.size() == detail.statistics.optimizationIterations,
            "Incomplete iteration trace");
    for (const auto &row : d.trace)
        trace << row.iteration << ','
              << (row.phase == RefinementPhase::Warmup ? "warmup" : "constrained") << ','
              << row.worst << ',' << row.cost << ',' << row.bestWorst << ',' << row.target << ','
              << row.trustBefore << ',' << row.trustAfter << ',' << row.stepNorm << ','
              << row.modelChange << ',' << row.stepSolved << ',' << row.accepted << '\n';
    Require(bool(trace), "Cannot write iteration trace");
    // Same plateau, deliberately restricted solves: a stopping reason must never
    // imply success.
    ConversionOptions disabled;
    disabled.refinement.maxIterations = 0;
    auto unoptimized = ConvertApproximate(*surface, disabled);
    Require(unoptimized.status != ConversionStatus::Success && !unoptimized.validation.Passed() &&
                unoptimized.statistics.optimizationIterations == 0 &&
                unoptimized.timing.validation.calls >= 2 &&
                unoptimized.statistics.refinementRounds == 1,
            "Zero-iteration plateau was accepted or repeatedly optimized without "
            "new samples");
    ConversionOptions tinyRange;
    tinyRange.refinement.minTrustRadius = .1;
    auto limited = ConvertApproximate(*surface, tinyRange);
    Require(limited.status != ConversionStatus::Success && !limited.validation.Passed() &&
                limited.statistics.optimizationIterations == 0 &&
                limited.statistics.refinementRuns.front().stopReason ==
                    RefinementStopReason::TrustTooSmall,
            "Trust-limited plateau was accepted or ignored the step range stop");
    // The same plateau at its known feasible midpoint has an unattainable 0.72 target.
    // This exercises constrained-phase stagnation and Jacobian reuse after rejected steps.
    auto cache = SampleSurface(*surface, {});
    ConversionOptions constrained;
    constrained.refinement.maxWarmupIterations = 0;
    auto midpoint = feasible;
    ConversionStatistics stationary;
    RefineCandidate(*surface, midpoint, cache, constrained, stationary);
    Require(stationary.optimizationIterations == 3 &&
                stationary.refinementRuns.front().warmupIterations == 0 &&
                stationary.refinementRuns.front().jacobianEvaluations == 1 &&
                stationary.refinementRuns.front().acceptedSteps == 0 &&
                stationary.refinementRuns.front().stopReason ==
                    RefinementStopReason::SampleFeasibleStagnation &&
                ValidateConversion(*surface, midpoint).Passed(),
            "Feasible plateau midpoint did not stop after three rejected constrained steps");
    // Unlike the zero-attempt check above, this reaches the shared post-failure range guard.
    constrained.refinement.minTrustRadius = .05;
    midpoint = feasible;
    ConversionStatistics shrinking;
    RefineCandidate(*surface, midpoint, cache, constrained, shrinking);
    Require(shrinking.optimizationIterations == 1 &&
                shrinking.refinementRuns.front().stopReason ==
                    RefinementStopReason::TrustTooSmall &&
                ValidateConversion(*surface, midpoint).Passed(),
            "A failed constrained step bypassed the trust-radius floor");
    std::cout << "plateau means(ms): classification=" << sums[0] / 5
              << " optimization=" << sums[1] / 5 << " validation=" << sums[2] / 5
              << " total=" << sums[3] / 5
              << "\nPlateau checks passed (one warm-up, five measured runs).\n";
    return 0;
}
int main(int argc, char **argv)
{
    sggk::init();
    int result = 1;
    try
    {
        result = Run(argc > 1 ? argv[1] : "plateau.csv");
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
    sggk::fini();
    return result;
}
