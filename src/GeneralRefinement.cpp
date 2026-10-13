#include "GeneralInternal.h"
namespace bspline2analytic::detail
{
RefinementProblem MotherRefinementProblem(const SampleCache &samples, const ConversionOptions &o)
{
    RefinementProblem problem;
    problem.evaluate = [&samples, &o](const AnalyticPatch &p, EvaluationWorkspace &w, ResidualBuffer &r) {
        ++w.evaluations;
        if (r.constraints.size() != samples.size()) ++w.resizes;
        r.constraints.resize(samples.size());
        r.leastSquares.resize(3 * samples.size());
        r.cost = r.worst = 0;
        for (size_t i = 0; i < samples.size(); ++i)
        {
            V d = (samples[i].point - ProjectToMother(p.parameters, samples[i].point)) /
                  o.validation.distanceTolerance;
            r.leastSquares[3 * i] = d.x;
            r.leastSquares[3 * i + 1] = d.y;
            r.leastSquares[3 * i + 2] = d.z;
            r.constraints[i] = Norm(d);
            r.cost += Dot(d, d);
            r.worst = std::max(r.worst, r.constraints[i]);
        }
        if (!std::isfinite(r.cost) || !std::isfinite(r.worst))
            throw Failure(Status::VerificationFailed, "Nonfinite general residual");
    };
    problem.increment = [](const AnalyticPatch &p, const std::vector<double> &d, double scale) {
        // The scratch chart never contributes residuals or output trims. Only
        // MotherSurface indices are active, so chart endpoints stay fixed.
        auto next = Increment(p, d, scale);
        if (next.parameters.type == sggk::SurfaceType::Torus &&
            !(next.parameters.majorRadius > next.parameters.radius))
            throw Failure(Status::InvalidInput, "Trial torus left the regular ring family");
        return next;
    };
    auto evaluate = problem.evaluate;
    auto increment = problem.increment;
    problem.jacobian = [evaluate, increment](const AnalyticPatch &p, const ActiveParameters &active,
        double scale, EvaluationWorkspace &w, RefinementDiagnostics &diagnostic) {
        const auto start = std::chrono::steady_clock::now();
        w.jacobian.assign(w.current.leastSquares.size(), std::vector<double>(active.indices.size()));
        w.constraintJacobian.assign(w.current.constraints.size(), std::vector<double>(active.indices.size()));
        ++w.resizes;
        std::vector<double> delta(ParametersFor(p.parameters.type, ParameterLevel::FullPatch).indices.size(), 0);
        bool ok = true;
        for (size_t j = 0; j < active.indices.size(); ++j)
        {
            constexpr double h = 2e-6;
            bool plus = false, minus = false;
            delta[active.indices[j]] = h;
            try { evaluate(increment(p, delta, scale), w, w.plus); plus = true; }
            catch (const Failure &) {}
            delta[active.indices[j]] = -h;
            try { evaluate(increment(p, delta, scale), w, w.minus); minus = true; }
            catch (const Failure &) {}
            delta[active.indices[j]] = 0;
            diagnostic.fallbackResidualEvaluations += unsigned(plus) + unsigned(minus);
            if (!plus && !minus) { ok = false; break; }
            const auto &a = plus ? w.plus : w.current;
            const auto &b = minus ? w.minus : w.current;
            double denominator = plus && minus ? 2 * h : h;
            for (size_t row = 0; row < w.jacobian.size(); ++row)
            {
                double value = (a.leastSquares[row] - b.leastSquares[row]) / denominator;
                w.jacobian[row][j] = value; ok &= std::isfinite(value);
            }
            for (size_t row = 0; row < w.constraintJacobian.size(); ++row)
            {
                // Zero distance uses a zero subgradient for its inactive constraint;
                // the three vector residual derivatives remain available to LM.
                double value = w.current.constraints[row] < 1e-8 ? 0 :
                    (a.constraints[row] - b.constraints[row]) / denominator;
                w.constraintJacobian[row][j] = value; ok &= std::isfinite(value);
            }
        }
        diagnostic.fallbackSampleBlocks += static_cast<unsigned>(w.current.constraints.size());
        diagnostic.jacobianMs += std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
        return ok;
    };
    problem.materialize = [](AnalyticPatch &p) { MakeNative(p); p.sdkDomain.reset(); };
    return problem;
}
}
