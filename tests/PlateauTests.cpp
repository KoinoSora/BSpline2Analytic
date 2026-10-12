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
static void Check(const ConversionResult &r, bool fast = true)
{
    Require(r.status == ConversionStatus::Success && r.validation.Passed(),
            "Plateau conversion did not pass independent validation");
    Require(r.patch && r.patch->parameters.type == sggk::SurfaceType::Plane,
            "Plateau did not recover a plane");
    Require(r.validation.geometry.upperBound < 1e-5 &&
                r.validation.boundaryTangent.upperBound < .1 &&
                r.validation.boundaryNormal.upperBound < .1,
            "Plateau tolerance exceeded");
    Require((r.statistics.optimizationIterations > 0 ||
             r.statistics.refinementRuns.front().planeOffsetSolves > 0) &&
                r.timing.validation.calls >= 2,
            "Plateau did not exercise optimization and revalidation");
    Require(r.statistics.optimizationIterations < 30, "Plateau still exhausted 30 iterations");
    Require(r.statistics.refinementRuns.size() == 1 && r.statistics.refinementRounds == 1,
            "Plateau unexpectedly restarted the optimizer");
    const auto &d = r.statistics.refinementRuns.front();
    Require(d.warmupIterations + d.constrainedIterations == r.statistics.optimizationIterations &&
                d.warmupIterations <= 5 && d.acceptedSteps <= r.statistics.optimizationIterations,
            "Refinement diagnostics do not account for attempts");
    if (fast)
        Require(d.planeOffsetSolves == 1 && d.scalarIterations == 0 && d.jacobianEvaluations == 0 &&
                    d.residualEvaluations == 2 && r.statistics.optimizationIterations == 0 &&
                    d.level == ParameterLevel::SizeOrOffset && d.activeDimension == 1 &&
                    std::abs(d.planeOffset - 8e-6) < 1e-10 && r.statistics.baseSamplingBuilds == 1,
                "Plateau did not use the direct midpoint solve");
    Require(d.bestWorst < 1 && d.bestWorst < d.initialWorst &&
                d.stopReason != RefinementStopReason::IterationLimit,
            "Plateau did not stop with a feasible sampled candidate");
    const auto &t = r.timing;
    Require(std::abs(t.classification.elapsedMs + t.optimization.elapsedMs +
                     t.validation.elapsedMs + t.otherMs - t.totalMs) < 1e-6,
            "Stage timing accounting mismatch");
}
// Reference path retains the pre-cache projection and boundary mapping helpers.
static void CheckEvaluation(const AnalyticPatch &p, const SampleCache &samples)
{
    EvaluationWorkspace w;
    ResidualBuffer r;
    EvaluateResiduals(p, samples, {}, w, r);
    unsigned resizes = w.resizes;
    EvaluateResiduals(p, samples, {}, w, r);
    Require(w.resizes == resizes, "Repeated evaluation resized stable buffers");
    auto frame = FrameOf(p.parameters);
    std::array<EdgeGeometry, 4> edges;
    for (unsigned k = 0; k < 4; ++k)
    {
        auto e = MakeEdge(p, k);
        double last = e.start;
        for (size_t j = 0; j < samples.boundaryOrder[k].size(); ++j)
        {
            last = EdgeCoordinate(e, frame.local(samples[samples.boundaryOrder[k][j]].point), last);
            if (!j)
                e.g0 = last;
        }
        e.g1 = last;
        edges[k] = e;
    }
    size_t row = 0;
    double cost = 0, worst = 0;
    auto compare = [&](V value, double metric) {
        V actual{r.leastSquares[3 * row], r.leastSquares[3 * row + 1], r.leastSquares[3 * row + 2]};
        Require(Norm(actual - value) < 5e-8 && std::abs(r.constraints[row] - metric) < 5e-8,
                "Cached evaluator changed normalized residuals");
        cost += Dot(value, value);
        worst = std::max(worst, metric);
        ++row;
    };
    for (size_t i = 0; i < samples.size(); ++i)
    {
        const auto &x = samples[i];
        V q;
        EdgePoint edgePoint;
        if (x.edge < 0)
            q = ProjectInteriorPoint(p, x.point);
        else
        {
            edgePoint = MapBoundaryPoint(p, edges[x.edge], x.point, x.fraction);
            q = edgePoint.point;
        }
        Require(Norm(q - w.mappedPoints[i]) < 1e-12, "Cached correspondence point changed");
        V d = (x.point - q) / 1e-5;
        compare(d, Norm(d));
        if (x.edge >= 0)
        {
            Require(Norm(edgePoint.tangent - w.mappedTangents[i]) < 1e-12,
                    "Cached correspondence tangent changed");
            double a = .1 * Pi / 180;
            compare((x.tangent - edgePoint.tangent) / (2 * std::sin(a / 2)),
                    Angle(x.tangent, edgePoint.tangent) / a);
        }
    }
    Require(std::abs(cost - r.cost) < 1e-7 * std::max(1.0, cost) &&
                std::abs(worst - r.worst) < 5e-8,
            "Cached aggregate errors changed");
}
static void CheckCachesAndFallbacks(const sggk::BSplineSurface &surface,
                                    const AnalyticPatch &initial)
{
    auto cache = SampleSurface(surface, {});
    CheckEvaluation(initial, cache);
    std::vector<double> delta(8, 0);
    delta[0] = 8e-6;
    delta[1] = 1e-4;
    delta[4] = 4e-5;
    CheckEvaluation(Increment(initial, delta, 1), cache);
    auto base = std::make_shared<SourceSampleCache>();
    base->points = cache.points;
    base->boundaryOrder = cache.boundaryOrder;
    base->keys = cache.keys;
    SampleCache a(base), b(base);
    Require(&a[0] == &b[0] && a.points.empty() && b.points.empty(), "Base samples were copied");
    ConversionValidationReport witness;
    witness.geometry.state = CheckState::Exceeded;
    double u = surface.MinParamU() + .371 * (surface.MaxParamU() - surface.MinParamU());
    double v = surface.MinParamV() + .423 * (surface.MaxParamV() - surface.MinParamV());
    witness.geometry.parameterRegion = {u, u, v, v};
    unsigned version = a.version;
    Require(AddWitnessSamples(surface, witness, {}, a) == 1 && a.version == version + 1 &&
                AddWitnessSamples(surface, witness, {}, a) == 0 && a.version == version + 1 &&
                b.size() == base->points.size(),
            "Feedback deduplication or cache isolation failed");
    CheckEvaluation(initial, a);
    std::shared_ptr<const SourceSampleCache> shared;
    ConversionResult result;
    for (unsigned test = 0; test < 2; ++test)
    {
        auto candidate = initial;
        if (!test)
        {
            std::vector<double> tilt(8, 0);
            tilt[1] = 1e-4;
            candidate = Increment(candidate, tilt, 1);
        }
        else
            candidate.domain[0] += 4e-5;
        MaterializePatch(candidate);
        size_t firstRun = result.statistics.refinementRuns.size();
        RefineAndValidateCandidate(surface, candidate, {}, result, shared, test + 1);
        Require(result.validation.Passed(), "Progressive plane fallback failed validation");
        bool upgraded = false, full = false;
        std::array<unsigned, 5> roundIterations{};
        for (size_t i = firstRun; i < result.statistics.refinementRuns.size(); ++i)
        {
            const auto &d = result.statistics.refinementRuns[i];
            upgraded |= d.entryReason == RefinementEntryReason::MoreVariables;
            full |= d.level == ParameterLevel::FullPatch;
            roundIterations.at(d.roundIndex - 1) += d.warmupIterations + d.constrainedIterations;
        }
        Require(upgraded && (!test || full),
                "Failed reduced solve did not open required variables");
        for (unsigned count : roundIterations)
            Require(count <= 30, "Variable levels reset the shared iteration budget");
        std::cout << "plateau fallback=" << test
                  << " distance=" << result.validation.geometry.upperBound
                  << " runs=" << result.statistics.refinementRuns.size() - firstRun << '\n';
    }
    Require(result.statistics.baseSamplingBuilds == 1, "Candidates resampled the shared source");
    auto clipped = initial;
    clipped.domain[0] += 4e-5;
    ConversionOptions bounded;
    bounded.refinement.maxScalarIterations = 1;
    bounded.refinement.scalarDistanceToleranceRatio = 1e-12;
    ConversionStatistics stats;
    EvaluationWorkspace w;
    RefinePlaneOffset(clipped, cache, bounded, stats, w);
    Require(stats.refinementRuns.front().scalarIterations == 1 &&
                stats.refinementRuns.front().stopReason ==
                    RefinementStopReason::ScalarIterationLimit &&
                !ValidateConversion(surface, clipped).Passed(),
            "Bounded scalar solve exceeded its budget or bypassed validation");
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
           "search_evaluations,best_worst,final_trust,stop_reason,plane_solves,scalar_iterations,"
           "residual_evaluations,workspace_resizes,active_dimension,offset\n";
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
            << ',' << d.finalTrust << ',' << RefinementStopReasonName(d.stopReason) << ','
            << d.planeOffsetSolves << ',' << d.scalarIterations << ',' << d.residualEvaluations
            << ',' << d.workspaceResizes << ',' << d.activeDimension << ',' << d.planeOffset
            << '\n';
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
    traced.refinement.jacobianMode = JacobianMode::FiniteDifference;
    traced.refinement.recordIterationTrace = true;
    traced.refinement.enablePlaneOffsetSolver = false;
    traced.refinement.enableProgressiveVariables = false;
    auto detail = ConvertApproximate(*surface, traced);
    Check(detail, false);
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
    ConversionOptions reduced;
    reduced.refinement.enablePlaneOffsetSolver = false;
    auto reducedResult = ConvertApproximate(*surface, reduced);
    Check(reducedResult, false);
    Require(reducedResult.statistics.refinementRuns.front().activeDimension == 1 &&
                reducedResult.statistics.refinementRuns.front().jacobianEvaluations > 0,
            "Radius/offset variable selection did not reach the generic solver");
    CheckCachesAndFallbacks(*surface, RecoverApproxPatch(*surface, exact.parameters));
    // Same plateau, deliberately restricted solves: a stopping reason must never
    // imply success.
    ConversionOptions disabled;
    disabled.refinement.maxIterations = 0;
    auto unoptimized = ConvertApproximate(*surface, disabled);
    Require(unoptimized.status != ConversionStatus::Success && !unoptimized.validation.Passed() &&
                unoptimized.statistics.optimizationIterations == 0 &&
                unoptimized.timing.validation.calls >= 2 &&
                unoptimized.statistics.refinementRounds <= unoptimized.statistics.candidateCount &&
                std::all_of(unoptimized.statistics.refinementRuns.begin(),
                            unoptimized.statistics.refinementRuns.end(),
                            [](const RefinementDiagnostics &d) { return d.roundIndex == 1; }),
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
    constrained.refinement.jacobianMode = JacobianMode::FiniteDifference;
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
