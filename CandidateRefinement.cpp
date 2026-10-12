#include "ApproximateInternal.h"
#include <deque>
#include <numeric>
namespace bspline2analytic::detail
{
static V Rotate(V p, V w)
{
    double a = Norm(w);
    if (a < 1e-14)
        return p + Cross(w, p);
    V n = w / a;
    return p * std::cos(a) + Cross(n, p) * std::sin(a) + n * (Dot(n, p) * (1 - std::cos(a)));
}
static unsigned Dimension(sggk::SurfaceType t)
{
    if (t == sggk::SurfaceType::Plane)
        return 8;
    if (t == sggk::SurfaceType::Cylinder)
        return 9;
    if (t == sggk::SurfaceType::Torus)
        return 11;
    return 10;
}
static AnalyticPatch Increment(const AnalyticPatch &base, const std::vector<double> &d,
                               double scale)
{
    AnalyticPatch p = base;
    auto &a = p.parameters;
    auto f = FrameOf(a);
    size_t i = 0;
    V shift{};
    if (a.type == sggk::SurfaceType::Plane)
        shift = f.z * (d[i++] * scale);
    else if (a.type == sggk::SurfaceType::Cylinder)
    {
        double x = d[i++], y = d[i++];
        shift = (f.x * x + f.y * y) * scale;
    }
    else
    {
        double x = d[i++], y = d[i++], z = d[i++];
        shift = (f.x * x + f.y * y + f.z * z) * scale;
    }
    a.origin = Point(f.o + shift);
    double rx = d[i++], ry = d[i++];
    V rotation = f.x * rx + f.y * ry;
    if (a.type == sggk::SurfaceType::Plane)
        rotation = rotation + f.z * d[i++];
    a.axis = Direction(Rotate(f.z, rotation));
    a.xDirection = Direction(Rotate(f.x, rotation));
    if (a.type == sggk::SurfaceType::Cone)
        a.halfAngle = std::atan(std::tan(a.halfAngle) * std::exp(d[i++]));
    else if (a.type != sggk::SurfaceType::Plane)
        a.radius *= std::exp(d[i++]);
    if (a.type == sggk::SurfaceType::Torus)
        a.majorRadius += d[i++] * scale;
    for (unsigned j = 0; j < 4; ++j)
    {
        bool length = a.type == sggk::SurfaceType::Plane ||
                      (j >= 2 && (a.type == sggk::SurfaceType::Cylinder ||
                                  a.type == sggk::SurfaceType::Cone));
        p.domain[j] += d[i++] * (length ? scale : 1);
    }
    if (!(p.domain[0] < p.domain[1] && p.domain[2] < p.domain[3]) ||
        (a.type != sggk::SurfaceType::Plane && p.domain[1] - p.domain[0] > 2 * Pi + 1e-12) ||
        (a.type == sggk::SurfaceType::Sphere && (p.domain[2] <= 0 || p.domain[3] >= Pi)) ||
        (a.type == sggk::SurfaceType::Cone && p.domain[2] <= 0))
        throw Failure(Status::InvalidInput, "Trial step left the regular chart");
    return p;
}
struct Residuals
{
    // leastSquares: normalized vector residuals; constraints: nonnegative
    // distance/tangent error ratios. Normal errors are acceptance-only.
    std::vector<double> leastSquares, constraints;
    // cost = sum of squared vector components; worst = max(constraints).
    double cost = 0, worst = 0;
};
static Residuals Evaluate(const AnalyticPatch &p, const SampleCache &samples,
                          const ValidationOptions &o,
                          const std::array<std::vector<size_t>, 4> &order)
{
    auto frame = FrameOf(p.parameters);
    std::array<EdgeGeometry, 4> edges;
    for (unsigned k = 0; k < 4; ++k)
    {
        auto e = MakeEdge(p, k);
        double last = e.start;
        for (size_t j = 0; j < order[k].size(); ++j)
        {
            last = EdgeCoordinate(e, frame.local(samples.points[order[k][j]].point), last);
            if (j == 0)
                e.g0 = last;
        }
        e.g1 = last;
        if (std::abs(e.g1 - e.g0) < 1e-12)
            throw Failure(Status::Degenerate, "Collapsed optimization correspondence");
        edges[k] = e;
    }
    Residuals r;
    auto add = [&](V a, double metric) {
        r.leastSquares.push_back(a.x);
        r.leastSquares.push_back(a.y);
        r.leastSquares.push_back(a.z);
        r.constraints.push_back(metric);
        r.cost += Dot(a, a);
        r.worst = std::max(r.worst, metric);
    };
    double ta = o.boundaryTangentToleranceDeg * Pi / 180;
    for (const auto &s : samples.points)
    {
        if (s.edge < 0)
        {
            V d = (s.point - ProjectInteriorPoint(p, s.point)) / o.distanceTolerance;
            add(d, Norm(d));
        }
        else
        {
            auto q = MapBoundaryPoint(p, edges[s.edge], s.point, s.fraction);
            V d = (s.point - q.point) / o.distanceTolerance;
            add(d, Norm(d));
            add((s.tangent - q.tangent) / (2 * std::sin(ta / 2)), Angle(s.tangent, q.tangent) / ta);
            // Boundary normals remain checked by ValidateConversion, but
            // contribute neither to the solver cost nor its max-error metric.
        }
    }
    if (!std::isfinite(r.cost) || !std::isfinite(r.worst))
        throw Failure(Status::VerificationFailed, "Nonfinite trial residual");
    return r;
}
// Minimum-norm step subject to linearized max-error constraints and a box
// trust region. Hildreth dual coordinate updates keep multipliers nonnegative.
static bool ConstrainedStep(const std::vector<std::vector<double>> &jac,
                            const std::vector<double> &residual, double target, double trust,
                            std::vector<double> &step)
{
    if (jac.empty())
        return false;
    size_t n = jac[0].size();
    auto rows = jac;
    std::vector<double> b;
    for (double r : residual)
        b.push_back(target - r);
    for (size_t j = 0; j < n; ++j)
        for (double sign : {-1.0, 1.0})
        {
            std::vector<double> row(n, 0);
            row[j] = sign;
            rows.push_back(row);
            b.push_back(trust);
        }
    std::vector<double> multiplier(rows.size(), 0), norm(rows.size(), 0);
    step.assign(n, 0);
    for (size_t i = 0; i < rows.size(); ++i)
        for (double a : rows[i])
            norm[i] += a * a;
    for (unsigned sweep = 0; sweep < 100; ++sweep)
    {
        double worst = 0;
        for (size_t i = 0; i < rows.size(); ++i)
        {
            double violation = -b[i];
            for (size_t j = 0; j < n; ++j)
                violation += rows[i][j] * step[j];
            worst = std::max(worst, violation);
            if (norm[i] < 1e-24)
                continue;
            double next = std::max(0.0, multiplier[i] + violation / norm[i]),
                   change = next - multiplier[i];
            for (size_t j = 0; j < n; ++j)
                step[j] -= change * rows[i][j];
            multiplier[i] = next;
        }
        if (worst < 1e-5)
            break;
    }
    for (double x : step)
        if (!std::isfinite(x) || std::abs(x) > trust * 1.01)
            return false;
    return true;
}
static double StepNorm(const std::vector<double> &step)
{
    double norm = 0;
    for (double x : step)
        norm = std::max(norm, std::abs(x));
    return norm;
}
static double ModelChange(const Residuals &a, const Residuals &b)
{
    double change = 0;
    for (size_t i = 0; i < a.leastSquares.size(); i += 3)
        change = std::max(change, std::hypot(a.leastSquares[i] - b.leastSquares[i],
                                             a.leastSquares[i + 1] - b.leastSquares[i + 1],
                                             a.leastSquares[i + 2] - b.leastSquares[i + 2]));
    return change;
}
RefinementResult RefineCandidate(const sggk::BSplineSurface &, AnalyticPatch &p,
                                 const SampleCache &samples, const ConversionOptions &o,
                                 ConversionStatistics &statistics)
{
    RefinementDiagnostics diagnostic;
    if (samples.points.empty())
    {
        statistics.refinementRuns.push_back(diagnostic);
        return {};
    }
    const auto &options = o.refinement;
    std::array<std::vector<size_t>, 4> order;
    double scale = 0;
    V origin = Get(p.parameters.origin);
    for (size_t i = 0; i < samples.points.size(); ++i)
    {
        const auto &x = samples.points[i];
        scale = std::max(scale, Norm(x.point - origin));
        if (x.edge >= 0)
            order[x.edge].push_back(i);
    }
    for (auto &list : order)
        std::sort(list.begin(), list.end(), [&](size_t a, size_t b) {
            return samples.points[a].fraction < samples.points[b].fraction;
        });
    scale = std::max(scale, 1e-8);
    auto current = Evaluate(p, samples, o.validation, order);
    auto best = current;
    AnalyticPatch bestPatch = p;
    diagnostic.initialWorst = current.worst;
    diagnostic.stopReason = RefinementStopReason::IterationLimit;
    double lambda = 1e-4, trust = 0.1;
    const unsigned n = Dimension(p.parameters.type);
    RefinementPhase phase =
        options.maxWarmupIterations ? RefinementPhase::Warmup : RefinementPhase::Constrained;
    unsigned failures = 0, smallChanges = 0;
    double phaseBestCost = current.cost;
    // A window of m attempts needs m+1 records, including the phase's initial
    // state.
    std::deque<std::pair<double, double>> history{{best.worst, phaseBestCost}};
    std::vector<std::vector<double>> J, C;
    bool jacobianValid = false;
    auto enterConstrained = [&]() {
        p = bestPatch;
        current = best;
        phase = RefinementPhase::Constrained;
        failures = smallChanges = 0;
        phaseBestCost = current.cost;
        history = {{best.worst, phaseBestCost}};
        jacobianValid = false;
    };
    for (unsigned iteration = 0; iteration < options.maxIterations; ++iteration)
    {
        if (best.worst <= 0.75)
        {
            diagnostic.stopReason = RefinementStopReason::MarginReached;
            break;
        }
        if (trust <= options.minTrustRadius)
        {
            diagnostic.stopReason = RefinementStopReason::TrustTooSmall;
            break;
        }
        ++statistics.optimizationIterations;
        if (phase == RefinementPhase::Warmup)
            ++diagnostic.warmupIterations;
        else
            ++diagnostic.constrainedIterations;
        RefinementIteration record;
        record.iteration = iteration + 1;
        record.phase = phase;
        record.trustBefore = trust;
        if (!jacobianValid)
        {
            ++diagnostic.jacobianEvaluations;
            J.assign(current.leastSquares.size(), std::vector<double>(n));
            C.assign(current.constraints.size(), std::vector<double>(n));
            for (unsigned j = 0; j < n; ++j)
            {
                const double h = 2e-6;
                std::vector<double> delta(n, 0);
                delta[j] = h;
                Residuals plus, minus;
                bool havePlus = false, haveMinus = false;
                try
                {
                    plus = Evaluate(Increment(p, delta, scale), samples, o.validation, order);
                    havePlus = true;
                }
                catch (const Failure &)
                {
                }
                delta[j] = -h;
                try
                {
                    minus = Evaluate(Increment(p, delta, scale), samples, o.validation, order);
                    haveMinus = true;
                }
                catch (const Failure &)
                {
                }
                if (!havePlus && !haveMinus)
                    continue;
                const auto &a = havePlus ? plus : current;
                const auto &b = haveMinus ? minus : current;
                double denominator = havePlus && haveMinus ? 2 * h : h;
                for (size_t i = 0; i < J.size(); ++i)
                    J[i][j] = (a.leastSquares[i] - b.leastSquares[i]) / denominator;
                for (size_t i = 0; i < C.size(); ++i)
                    C[i][j] = (a.constraints[i] - b.constraints[i]) / denominator;
            }
            jacobianValid = true;
        }
        std::vector<double> step;
        const bool constrained = phase == RefinementPhase::Constrained;
        if (constrained)
        {
            record.target = std::max(0.70, current.worst * 0.9);
            record.stepSolved = ConstrainedStep(C, current.constraints, record.target, trust, step);
        }
        else
        {
            auto rows = J;
            std::vector<double> rhs;
            for (double r : current.leastSquares)
                rhs.push_back(-r);
            for (unsigned j = 0; j < n; ++j)
            {
                double norm = 0;
                for (const auto &row : J)
                    norm = std::hypot(norm, row[j]);
                std::vector<double> row(n, 0);
                row[j] = std::sqrt(lambda) * std::max(norm, 1.0);
                rows.push_back(row);
                rhs.push_back(0);
            }
            record.stepSolved = SolveLeastSquares(rows, rhs, step);
            if (record.stepSolved && StepNorm(step) > trust)
            {
                const double factor = trust / StepNorm(step);
                for (auto &x : step)
                    x *= factor;
            }
        }
        if (record.stepSolved)
            for (unsigned line = 0; line < 12; ++line)
            {
                record.stepNorm = StepNorm(step);
                if (record.stepNorm == 0)
                    break;
                try
                {
                    auto trial = Increment(p, step, scale);
                    ++diagnostic.lineSearchEvaluations;
                    auto next = Evaluate(trial, samples, o.validation, order);
                    record.modelChange = ModelChange(current, next);
                    // An unchanged sampled model cannot improve either acceptance
                    // objective.
                    if (record.modelChange == 0)
                        break;
                    if ((constrained
                             ? next.worst < current.worst - 1e-8
                             : next.cost < current.cost && next.worst < current.worst * 1.2))
                    {
                        p = std::move(trial);
                        current = std::move(next);
                        record.accepted = true;
                        jacobianValid = false;
                        ++diagnostic.acceptedSteps;
                        lambda = std::max(1e-10, lambda * 0.3);
                        if (current.worst < best.worst ||
                            (current.worst == best.worst && current.cost < best.cost))
                        {
                            best = current;
                            bestPatch = p;
                        }
                        break;
                    }
                }
                catch (const Failure &)
                {
                }
                for (auto &x : step)
                    x *= 0.5;
            }
        if (record.accepted)
        {
            failures = 0;
            smallChanges =
                record.modelChange < options.smallModelChangeTolerance ? smallChanges + 1 : 0;
        }
        else
        {
            ++failures;
            smallChanges = 0;
            lambda *= 10;
            trust *= 0.5;
        }
        phaseBestCost = std::min(phaseBestCost, current.cost);
        history.push_back({best.worst, phaseBestCost});
        if (history.size() > size_t(options.stagnationWindow) + 1)
            history.pop_front();
        const bool fullWindow = history.size() == size_t(options.stagnationWindow) + 1;
        const bool stalled =
            fullWindow && history.front().first - best.worst <=
                              options.worstImprovementAbsTol +
                                  options.worstImprovementRelTol * history.front().first;
        const bool costStalled = fullWindow && history.front().second - phaseBestCost <=
                                                   options.costImprovementRelTol *
                                                       std::max(history.front().second, 1e-30);
        record.worst = current.worst;
        record.cost = current.cost;
        record.bestWorst = best.worst;
        record.trustAfter = trust;
        diagnostic.lastStepNorm = record.stepNorm;
        if (options.recordIterationTrace)
            diagnostic.trace.push_back(record);
        // Every failure path reaches these checks, including a failed constrained
        // solve.
        if (best.worst <= 0.75)
        {
            diagnostic.stopReason = RefinementStopReason::MarginReached;
            break;
        }
        if (trust <= options.minTrustRadius)
        {
            diagnostic.stopReason = RefinementStopReason::TrustTooSmall;
            break;
        }
        if (stalled && best.worst < 1)
        {
            diagnostic.stopReason = RefinementStopReason::SampleFeasibleStagnation;
            break;
        }
        if (!constrained)
        {
            if ((stalled && costStalled) || failures >= options.maxConsecutiveFailures ||
                diagnostic.warmupIterations >= options.maxWarmupIterations)
                enterConstrained();
        }
        else
        {
            if (failures >= options.maxConsecutiveFailures)
            {
                diagnostic.stopReason = record.stepSolved ? RefinementStopReason::NoAcceptedStep
                                                          : RefinementStopReason::SolverFailure;
                break;
            }
            if (stalled)
            {
                diagnostic.stopReason = smallChanges >= options.stagnationWindow
                                            ? RefinementStopReason::SmallStep
                                            : RefinementStopReason::Stagnation;
                break;
            }
        }
    }
    if (best.worst <= 0.75)
        diagnostic.stopReason = RefinementStopReason::MarginReached;
    p = std::move(bestPatch);
    diagnostic.bestWorst = best.worst;
    diagnostic.bestCost = best.cost;
    diagnostic.finalTrust = trust;
    RefinementResult result{best.worst < diagnostic.initialWorst, diagnostic.stopReason};
    // Always materialize the best sampled candidate; success still requires full
    // validation.
    MakeNative(p);
    for (unsigned i = 0; i < 4; ++i)
    {
        unsigned dim = p.mapping.swappedUV ? 1 : 0;
        std::array<double, 2> uv{};
        uv[dim] = p.domain[2 * dim + (bool(i % 2) != p.mapping.reversedU ? 1 : 0)];
        uv[1 - dim] = p.domain[2 * (1 - dim) + (bool(i / 2) != p.mapping.reversedV ? 1 : 0)];
        p.mapping.analyticCorners[i] = sggk::Point2D(uv[0], uv[1]);
    }
    statistics.refinementRuns.push_back(std::move(diagnostic));
    return result;
}
} // namespace bspline2analytic::detail
