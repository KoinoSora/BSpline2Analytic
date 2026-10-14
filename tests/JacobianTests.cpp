#include "ApproximateCases.h"
#include "SpecializedTimingChecks.h"
#include "ApproximateInternal.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
static unsigned checks = 0;
static void Require(bool ok, const std::string &message)
{
    ++checks;
    if (!ok)
        throw std::runtime_error(message);
}
static double MaxDifference(const std::vector<double> &a, const std::vector<double> &b)
{
    double r = 0;
    for (size_t i = 0; i < a.size(); ++i)
        r = std::max(r, std::abs(a[i] - b[i]));
    return r;
}
static void Derivatives(const tests::ReferenceCase &reference, bool transformed,
                        std::ostream &kernel)
{
    auto s = tests::RefitPolynomial(*reference.spline, 4);
    if (transformed)
    {
        s->ExchangeUV();
        s->UReverse();
        s->VReverse();
    }
    auto exact = Recognize(*reference.spline);
    Require(exact.status == Status::Success, "Reference recognition");
    auto p = RecoverApproxPatch(*s, exact.parameters);
    auto samples = SampleSurface(*s, {});
    auto full = ParametersFor(p.parameters.type, ParameterLevel::FullPatch);
    double scale = SampleScale(p, samples);
    std::vector<double> perturb(full.indices.size());
    for (size_t j = 0; j < perturb.size(); ++j)
        perturb[j] = .002 * std::sin(.83 * (j + 1));
    p = Increment(p, perturb, scale);
    // Include stable finite-edge and finite-corner projections, in addition to
    // the interior and boundary correspondence blocks from the source cache.
    auto f = FrameOf(p.parameters);
    for (auto uv :
         {std::array<double, 2>{p.domain[0] - .13, p.domain[2] + .37 * (p.domain[3] - p.domain[2])},
          std::array<double, 2>{p.domain[0] - .17, p.domain[2] - .11}})
    {
        Sample sample;
        sample.point = f.global(CanonicalPoint(p.parameters, uv[0], uv[1]));
        samples.points.push_back(sample);
    }
    samples.RebuildIndices();
    ConversionOptions o;
    EvaluationWorkspace w;
    RefinementDiagnostics diagnostic;
    EvaluateResiduals(p, samples, o.validation, w, w.current);
    Require(BuildJacobian(p, samples, o, full, scale, w, diagnostic),
            "Analytic linearization failed");
    Require(diagnostic.analyticSampleBlocks > 0 && diagnostic.fallbackSampleBlocks == 0,
            "Regular case unexpectedly fell back");
    double worstJ = 0, worstC = 0;
    // Multi-step central differences, independently evaluating the existing
    // value path (not the derivative formulas).
    for (double h : {2e-6, 4e-7, 8e-8})
    {
        double ej = 0, ec = 0;
        for (size_t j = 0; j < full.indices.size(); ++j)
        {
            std::vector<double> d(full.indices.size(), 0);
            EvaluationWorkspace scratch;
            d[j] = h;
            EvaluateResiduals(Increment(p, d, scale), samples, o.validation, scratch, scratch.plus);
            d[j] = -h;
            EvaluateResiduals(Increment(p, d, scale), samples, o.validation, scratch,
                              scratch.minus);
            double diffJ = 0, diffC = 0, normJ = 1, normC = 1;
            for (size_t i = 0; i < w.jacobian.size(); ++i)
            {
                double fd =
                    (scratch.plus.leastSquares[i] - scratch.minus.leastSquares[i]) / (2 * h);
                diffJ = std::max(diffJ, std::abs(fd - w.jacobian[i][j]));
                normJ = std::max(normJ, std::max(std::abs(fd), std::abs(w.jacobian[i][j])));
            }
            for (size_t i = 0; i < w.constraintJacobian.size(); ++i)
            {
                double fd = (scratch.plus.constraints[i] - scratch.minus.constraints[i]) / (2 * h);
                diffC = std::max(diffC, std::abs(fd - w.constraintJacobian[i][j]));
                normC =
                    std::max(normC, std::max(std::abs(fd), std::abs(w.constraintJacobian[i][j])));
            }
            ej = std::max(ej, diffJ / normJ);
            ec = std::max(ec, diffC / normC);
        }
        worstJ = ej;
        worstC = ec;
        std::cout << reference.name << " derivative h=" << h << " J=" << ej << " C=" << ec
                  << std::endl;
    }
    Require(worstJ < 5e-7 && worstC < 5e-7, "Analytic/finite-difference mismatch");
    std::vector<double> direction(full.indices.size());
    for (size_t j = 0; j < direction.size(); ++j)
        direction[j] = .3 * std::cos(.71 * (j + 1));
    double previous = 0;
    for (double h : {1e-4, 5e-5, 2.5e-5})
    {
        auto d = direction;
        for (auto &x : d)
            x *= h;
        EvaluationWorkspace trial;
        EvaluateResiduals(Increment(p, d, scale), samples, o.validation, trial, trial.trial);
        auto prediction = w.current.leastSquares;
        for (size_t i = 0; i < prediction.size(); ++i)
            for (size_t j = 0; j < d.size(); ++j)
                prediction[i] += w.jacobian[i][j] * d[j];
        double error = MaxDifference(prediction, trial.trial.leastSquares);
        if (previous)
            Require(error < previous * .30, "Directional Taylor remainder not quadratic");
        previous = error;
    }
    // Trial evaluations must not poison the current linearization geometry.
    auto saved = w.jacobian;
    EvaluateResiduals(Increment(p, perturb, scale), samples, o.validation, w, w.trial);
    Require(BuildJacobian(p, samples, o, full, scale, w, diagnostic),
            "Linearization after rejected trial failed");
    for (size_t i = 0; i < saved.size(); ++i)
        Require(MaxDifference(saved[i], w.jacobian[i]) < 1e-12, "Stale trial geometry reused");
    for (auto level : {ParameterLevel::SizeOrOffset, ParameterLevel::MotherSurface})
    {
        auto active = ParametersFor(p.parameters.type, level);
        Require(BuildJacobian(p, samples, o, active, scale, w, diagnostic),
                "Reduced linearization failed");
        for (size_t i = 0; i < saved.size(); ++i)
            for (size_t j = 0; j < active.indices.size(); ++j)
                Require(std::abs(saved[i][active.indices[j]] - w.jacobian[i][j]) < 1e-8,
                        "Reduced derivative seed differs from full seed");
    }
    if (!transformed)
        for (auto mode : {JacobianMode::AnalyticWithFallback, JacobianMode::FiniteDifference})
        {
            ConversionOptions timed;
            timed.refinement.jacobianMode = mode;
            RefinementDiagnostics warmup, measured;
            Require(BuildJacobian(p, samples, timed, full, scale, w, warmup),
                    "Kernel warmup failed");
            unsigned initial = w.evaluations;
            for (unsigned repeat = 0; repeat < 20; ++repeat)
                Require(BuildJacobian(p, samples, timed, full, scale, w, measured),
                        "Kernel benchmark failed");
            kernel << reference.name << ','
                   << (mode == JacobianMode::FiniteDifference ? "finite_difference" : "analytic")
                   << ',' << full.indices.size() << ',' << samples.size() << ','
                   << measured.jacobianMs / 20 << ',' << (w.evaluations - initial) / 20.0 << '\n';
        }
}
static void Nonsmooth(const tests::ReferenceCase &reference)
{
    auto s = tests::RefitPolynomial(*reference.spline, 2);
    auto exact = Recognize(*reference.spline);
    auto p = RecoverApproxPatch(*s, exact.parameters);
    auto samples = SampleSurface(*s, {});
    auto f = FrameOf(p.parameters);
    Sample edge;
    edge.point =
        f.global(CanonicalPoint(p.parameters, p.domain[0], (p.domain[2] + p.domain[3]) * .5));
    samples.points.push_back(edge);
    samples.RebuildIndices();
    ConversionOptions o;
    EvaluationWorkspace w;
    RefinementDiagnostics d;
    auto active = ParametersFor(p.parameters.type, ParameterLevel::FullPatch);
    double scale = SampleScale(p, samples);
    EvaluateResiduals(p, samples, o.validation, w, w.current);
    Require(BuildJacobian(p, samples, o, active, scale, w, d), "Nonsmooth fallback failed");
    Require(d.fallbackSampleBlocks > 0 && d.analyticSampleBlocks > 0 &&
                d.fallbackResidualEvaluations > 0 && d.nonsmoothConstraintRows > 0,
            "Fallback/zero-subgradient diagnostics missing");
    auto row = w.jacobian.size() - 3;
    auto fallback = w.jacobian;
    o.refinement.jacobianMode = JacobianMode::FiniteDifference;
    Require(BuildJacobian(p, samples, o, active, scale, w, d), "Difference reference failed");
    for (size_t i = row; i < row + 3; ++i)
        Require(MaxDifference(fallback[i], w.jacobian[i]) < 1e-8,
                "Local fallback differs from full difference");
}
static void OneSided(const tests::ReferenceCase &reference)
{
    auto s = tests::RefitPolynomial(*reference.spline, 2);
    auto p = RecoverApproxPatch(*s, Recognize(*reference.spline).parameters);
    p.domain[1] = p.domain[0] + 2 * Pi;
    auto samples = SampleSurface(*s, {});
    auto f = FrameOf(p.parameters);
    Sample point;
    point.point =
        f.global(CanonicalPoint(p.parameters, p.domain[0], (p.domain[2] + p.domain[3]) * .5));
    samples.points.push_back(point);
    samples.RebuildIndices();
    ActiveParameters active{{6}}; // Cylinder umax: the positive step exceeds one turn.
    ConversionOptions o;
    EvaluationWorkspace w;
    RefinementDiagnostics d;
    EvaluateResiduals(p, samples, o.validation, w, w.current);
    Require(BuildJacobian(p, samples, o, active, 1, w, d), "One-sided analytic fallback failed");
    Require(d.fallbackSampleBlocks > 0 && d.fallbackResidualEvaluations == 2,
            "Fallback must count only base and valid negative-side evaluations");
    auto analytic = w.jacobian;
    o.refinement.jacobianMode = JacobianMode::FiniteDifference;
    unsigned before = w.evaluations;
    Require(BuildJacobian(p, samples, o, active, 1, w, d) && w.evaluations - before == 1,
            "Invalid side was evaluated or one-sided difference failed");
    for (size_t row = analytic.size() - 3; row < analytic.size(); ++row)
        Require(MaxDifference(analytic[row], w.jacobian[row]) < 1e-8,
                "One-sided local/full fallback mismatch");
}
static void Perturbations(const std::vector<tests::ReferenceCase> &references, const char *path)
{
    std::ofstream csv(path);
    Require(bool(csv), "Cannot open perturbation CSV");
    csv << "case,mode,run,classify_ms,optimize_ms,validate_ms,total_ms,jacobian_ms,iterations,"
           "jacobians,value_evaluations,analytic_blocks,fallback_blocks,distance_bound,tangent_deg,"
           "normal_deg\n"
        << std::setprecision(17);
    for (unsigned k = 0; k < 5; ++k)
    {
        auto fixture = tests::CurvedPlateauFixture(references[k]);
        Require(!fixture.spline->IsRational(), "Perturbation fixture is rational");
        auto known = ValidateConversion(*fixture.spline, fixture.feasible);
        std::cout << fixture.name << " known feasible=" << known.Passed()
                  << " bounds=" << known.geometry.upperBound << ","
                  << known.boundaryTangent.upperBound << "," << known.boundaryNormal.upperBound
                  << std::endl;
        Require(known.Passed(), fixture.name + " has no certified known feasible solution");
        Require(!ValidateConversion(*fixture.spline, fixture.initial).Passed(),
                fixture.name + " perturbation did not require refinement");
        for (auto mode : {JacobianMode::AnalyticWithFallback, JacobianMode::FiniteDifference})
        {
            ConversionOptions o;
            o.refinement.jacobianMode = mode;
            // Exercise Jacobians on the plane too; the default direct solver is
            // separately covered by the Plateau regression.
            o.refinement.enablePlaneOffsetSolver = false;
            for (int run = -1; run < 3; ++run)
            {
                auto r = ConvertApproximate(*fixture.spline, o);
                std::cout << fixture.name << " "
                          << (mode == JacobianMode::FiniteDifference ? "FD" : "analytic")
                          << " run=" << run << " " << ConversionStatusName(r.status)
                          << " bounds=" << r.validation.geometry.upperBound << ","
                          << r.validation.boundaryTangent.upperBound << ","
                          << r.validation.boundaryNormal.upperBound << " message=" << r.message
                          << " iterations=" << r.statistics.optimizationIterations
                          << " opt=" << r.timing.optimization.elapsedMs << std::endl;
                if (!r.validation.Passed())
                {
                    std::cout << "states " << CheckStateName(r.validation.geometry.state) << ","
                              << CheckStateName(r.validation.boundaryTangent.state) << ","
                              << CheckStateName(r.validation.boundaryNormal.state)
                              << " samples=" << r.statistics.sampleCount
                              << " pieces=" << r.statistics.validationPieces << std::endl;
                    for (const auto &d : r.statistics.refinementRuns)
                        std::cout << "level=" << int(d.level) << " initial=" << d.initialWorst
                                  << " best=" << d.bestWorst
                                  << " stop=" << RefinementStopReasonName(d.stopReason)
                                  << " round=" << d.roundIndex << std::endl;
                    if (r.patch)
                    {
                        std::cout << "origin diff="
                                  << Norm(Get(r.patch->parameters.origin) -
                                          Get(fixture.feasible.parameters.origin))
                                  << " radius diff="
                                  << r.patch->parameters.radius - fixture.feasible.parameters.radius
                                  << std::endl;
                        for (unsigned j = 0; j < 4; ++j)
                            std::cout << "domain " << j << " " << r.patch->domain[j] << " feasible "
                                      << fixture.feasible.domain[j] << std::endl;
                    }
                }
                Require(r.status == ConversionStatus::Success && r.validation.Passed(),
                        fixture.name + " perturbed conversion failed");
                Require(r.patch->parameters.type == fixture.initial.parameters.type,
                        "Wrong perturbed type");
                Require(r.validation.geometry.upperBound < 1e-5 &&
                            r.validation.boundaryTangent.upperBound < .1 &&
                            r.validation.boundaryNormal.upperBound < .1,
                        "Tolerance not strict");
                double jm = 0;
                unsigned jc = 0, ve = 0, ab = 0, fb = 0;
                for (const auto &d : r.statistics.refinementRuns)
                {
                    jm += d.jacobianMs;
                    jc += d.jacobianEvaluations;
                    ve += d.residualEvaluations;
                    ab += d.analyticSampleBlocks;
                    fb += d.fallbackSampleBlocks;
                }
                Require(jc > 0 && r.statistics.optimizationIterations > 0,
                        "Jacobian not exercised");
                Require(mode == JacobianMode::FiniteDifference || ab > 0,
                        "Analytic mode not exercised");
                const auto &t = r.timing;
                tests::CheckSpecializedTiming(t);
                Require(std::abs(t.classification.elapsedMs + t.optimization.elapsedMs +
                                 t.validation.elapsedMs + t.otherMs - t.totalMs) < 1e-6,
                        "Stage timing mismatch");
                Require(jm <= t.optimization.elapsedMs + 1e-6,
                        "Jacobian timer exceeds optimization timer");
                if (run >= 0)
                    csv << fixture.name << ','
                        << (mode == JacobianMode::FiniteDifference ? "finite_difference"
                                                                   : "analytic")
                        << ',' << run << ',' << t.classification.elapsedMs << ','
                        << t.optimization.elapsedMs << ',' << t.validation.elapsedMs << ','
                        << t.totalMs << ',' << jm << ',' << r.statistics.optimizationIterations
                        << ',' << jc << ',' << ve << ',' << ab << ',' << fb << ','
                        << r.validation.geometry.upperBound << ','
                        << r.validation.boundaryTangent.upperBound << ','
                        << r.validation.boundaryNormal.upperBound << '\n';
            }
        }
        ConversionOptions disabled;
        disabled.refinement.maxIterations = 0;
        auto unoptimized = fixture.initial;
        ConversionStatistics stopped;
        auto sourceSamples = SampleSurface(*fixture.spline, {});
        RefineCandidate(*fixture.spline, unoptimized, sourceSamples, disabled, stopped);
        Require(stopped.optimizationIterations == 0 &&
                    !ValidateConversion(*fixture.spline, unoptimized).Passed(),
                "Zero iteration budget changed the invalid candidate");
        ValidationOptions budget;
        budget.maxValidationPieces = 1;
        Require(!ValidateConversion(*fixture.spline, fixture.feasible, budget).Passed(),
                "Insufficient certificate budget falsely accepted");
    }
    Require(bool(csv), "Cannot write perturbation CSV");
}
int main(int argc, char **argv)
{
    sggk::init();
    int result = 1;
    try
    {
        auto cases = tests::MakeCases();
        std::ofstream kernel(std::string(argc > 1 ? argv[1] : "perturbation.csv") +
                             ".jacobian.csv");
        Require(bool(kernel), "Cannot open kernel CSV");
        kernel << "case,mode,dimension,samples,jacobian_ms,value_evaluations\n"
               << std::setprecision(17);
        for (unsigned k = 0; k < 5; ++k)
        {
            Derivatives(cases[k], false, kernel);
            Derivatives(cases[k], true, kernel);
        }
        Nonsmooth(cases[0]);
        OneSided(cases[1]);
        Perturbations(cases, argc > 1 ? argv[1] : "perturbation.csv");
        std::cout << "Jacobian/perturbation checks passed: " << checks << std::endl;
        result = 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << std::endl;
    }
    sggk::fini();
    return result;
}
