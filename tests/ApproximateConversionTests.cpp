#include "ApproximateCases.h"
#include "ApproximateInternal.h"
#include <Foundation/init.h>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
static unsigned checks = 0;
static void Require(bool condition, const std::string &message)
{
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}
static void Timing(const ConversionResult &r)
{
    const auto &t = r.timing;
    Require(t.totalMs >= 0 && t.otherMs >= -1e-9, "Negative timing");
    double sum =
        t.classification.elapsedMs + t.optimization.elapsedMs + t.validation.elapsedMs + t.otherMs;
    Require(std::abs(sum - t.totalMs) < 1e-8 * std::max(1.0, t.totalMs),
            "Timing partition mismatch");
    if (!t.optimization.calls)
        Require(t.optimization.elapsedMs == 0, "Unused optimization time");
}

using namespace bspline2analytic::detail;
static void Print(const std::string &name, const ConversionResult &r)
{
    std::cout << name << " " << ConversionStatusName(r.status) << " " << r.message
              << " bounds=" << r.validation.geometry.upperBound << ","
              << r.validation.boundaryTangent.upperBound << ","
              << r.validation.boundaryNormal.upperBound
              << " pieces=" << r.statistics.validationPieces
              << " iterations=" << r.statistics.optimizationIterations
              << " times=" << r.timing.classification.elapsedMs << ","
              << r.timing.optimization.elapsedMs << "," << r.timing.validation.elapsedMs
              << std::endl;
}
static sggk::BSplineSurfacePtr PlaneFixture(unsigned mode)
{
    sggk::Point3DMatrix cp(4, std::vector<sggk::Point3D>(4));
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            cp[i][j] = sggk::Point3D(i * .001 / 3, j / 3.0, 0);
    if (mode == 1)
        for (unsigned j = 1; j < 3; ++j)
            cp[1][j] = sggk::Point3D(.001 / 3, j / 3.0, 8e-6);
    if (mode == 2)
        for (unsigned i = 1; i < 3; ++i)
            cp[i][0] = sggk::Point3D(i * .001 / 3, 8e-6, 0);
    if (mode == 3)
        for (unsigned i = 1; i < 3; ++i)
            for (unsigned j = 1; j < 3; ++j)
                cp[i][j] = sggk::Point3D(i * .001 / 3, j / 3.0, 1e-3);
    return std::make_shared<sggk::BSplineSurface>(3, 3, cp, sggk::RealArray{0, 1},
                                                  sggk::RealArray{0, 1}, sggk::UIntArray{4, 4},
                                                  sggk::UIntArray{4, 4});
}
static void RejectionTests()
{
    AnalyticPatch p;
    p.parameters.type = sggk::SurfaceType::Plane;
    p.domain = {0, .001, 0, 1};
    MakeNative(p);
    auto n = ValidateConversion(*PlaneFixture(1), p);
    Require(n.geometry.state == CheckState::Passed, "Normal fixture distance not passed");
    Require(n.boundaryTangent.state == CheckState::Passed, "Normal fixture tangent not passed");
    Require(n.boundaryNormal.state == CheckState::Exceeded, "Normal-only violation missed");
    auto t = ValidateConversion(*PlaneFixture(2), p);
    Require(t.geometry.state == CheckState::Passed, "Tangent fixture distance not passed");
    Require(t.boundaryTangent.state == CheckState::Exceeded, "Tangent-only violation missed");
    Require(t.boundaryNormal.state == CheckState::Passed, "Tangent fixture normal not passed");
    auto d = ValidateConversion(*PlaneFixture(3), p);
    Require(d.geometry.state == CheckState::Exceeded && d.geometry.hasWitness,
            "Interior bulge missed");
    Require(d.geometry.upperBound >= 5e-4, "Failed report lost the region upper bound");
    ValidationOptions budget;
    budget.maxValidationPieces = 1;
    auto b = ValidateConversion(*PlaneFixture(0), p, budget);
    Require(!b.Passed() && b.boundaryNormal.state == CheckState::Unresolved,
            "Budget exhaustion falsely accepted");
    ConversionOptions invalid;
    invalid.sampling.boundaryPointsPerSpan = 0;
    auto bad = ConvertApproximate(*PlaneFixture(0), invalid);
    Timing(bad);
    Require(bad.status == ConversionStatus::InvalidInput && bad.timing.classification.calls == 0,
            "Invalid options accepted");
    invalid = {};
    invalid.validation.distanceTolerance = 0;
    bad = ConvertApproximate(*PlaneFixture(0), invalid);
    Timing(bad);
    Require(bad.status == ConversionStatus::InvalidInput, "Invalid distance tolerance accepted");
}
static void IntervalTests()
{
    Require(IntervalEnvironmentSupported(), "Unsupported floating-point environment");
    auto invalid = I(INFINITY) - I(INFINITY);
    Require(invalid.lo == -INFINITY && invalid.hi == INFINITY,
            "Indeterminate interval subtraction narrowed");
    BPoly nonfinite(1, 1, I(std::numeric_limits<double>::quiet_NaN()));
    Require(nonfinite.range().hi == INFINITY, "NaN polynomial range narrowed");
    auto contains = [](I a, double x) { return a.lo <= x && x <= a.hi; };
    for (int i = -100; i <= 100; ++i)
    {
        double x = i * .071;
        Require(contains(Sin(I(x)), std::sin(x)), "Sin enclosure");
        Require(contains(Cos(I(x)), std::cos(x)), "Cos enclosure");
        Require(contains(Atan(I(x)), std::atan(x)), "Atan enclosure");
        double y = i * .01;
        Require(contains(Acos(I(y)), std::acos(y)), "Acos enclosure");
    }
    sggk::Point3DMatrix cp(6, std::vector<sggk::Point3D>(5));
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 5; ++j)
            cp[i][j] = sggk::Point3D(i * .13, j * .21, std::sin(i + .7 * j) * .05);
    sggk::BSplineSurface s(3, 3, cp, sggk::RealArray{0, .2, .75, 1}, sggk::RealArray{0, .63, 1},
                           sggk::UIntArray{4, 1, 1, 4}, sggk::UIntArray{4, 1, 4});
    auto cells = ExtractCells(s);
    Require(cells.size() == 6, "Extraction cell count");
    for (const auto &cell : cells)
        for (double u : {.13, .47, .81})
            for (double v : {.17, .53, .89})
            {
                auto q = Range(Evaluate(cell.p, u, v));
                auto actual = Get(s.CalcPoint(cell.u0 + (cell.u1 - cell.u0) * u,
                                              cell.v0 + (cell.v1 - cell.v0) * v));
                // SDK evaluation and parameter interpolation also round; this allowance
                // is for the independent comparison, never for the certificate.
                auto close = [](I a, double x) {
                    double slack =
                        32 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(x));
                    return x >= a.lo - slack && x <= a.hi + slack;
                };
                Require(close(q.x, actual.x) && close(q.y, actual.y) && close(q.z, actual.z),
                        "Extracted cell disagrees with SDK");
            }
}

int main()
{
    sggk::init();
    int code = 0;
    try
    {
        IntervalTests();
        RejectionTests();
        auto cases = tests::MakeCases();
        for (unsigned k = 0; k < 5; ++k)
        {
            auto s = tests::RefitPolynomial(*cases[k].spline, 12);
            Require(!s->IsRational(), "Fixture is rational");
            auto r = ConvertApproximate(*s);
            Timing(r);
            std::cout << cases[k].name << " " << ConversionStatusName(r.status) << " " << r.message
                      << " geometry=" << CheckStateName(r.validation.geometry.state) << ":"
                      << r.validation.geometry.upperBound
                      << " tangent=" << CheckStateName(r.validation.boundaryTangent.state) << ":"
                      << r.validation.boundaryTangent.upperBound
                      << " normal=" << CheckStateName(r.validation.boundaryNormal.state) << ":"
                      << r.validation.boundaryNormal.upperBound
                      << " pieces=" << r.statistics.validationPieces
                      << " times=" << r.timing.classification.elapsedMs << ","
                      << r.timing.optimization.elapsedMs << "," << r.timing.validation.elapsedMs
                      << std::endl;
            Require(r.status == ConversionStatus::Success, cases[k].name + " conversion failed");
            Require(r.patch->parameters.type == cases[k].analytic->SurfaceType(),
                    "Wrong analytic type");
            Require(r.validation.Passed(), "Success without all certificates");

            auto f = FrameOf(r.patch->parameters);
            for (double a : {0.0, .3, 1.0})
                for (double b : {0.0, .7, 1.0})
                {
                    double u = r.patch->domain[0] + (r.patch->domain[1] - r.patch->domain[0]) * a;
                    double v = r.patch->domain[2] + (r.patch->domain[3] - r.patch->domain[2]) * b;
                    Require(Norm(f.global(CanonicalPoint(r.patch->parameters, u, v)) -
                                 Get(r.patch->surface->CalcPoint(u, v))) < 1e-10,
                            "Native analytic parameter mismatch");
                }
            auto perturbed = r.patch->parameters;
            perturbed.origin = Point(Get(perturbed.origin) + f.z * 1e-4 + f.x * 2e-5);
            if (k != 0 && k != 2)
                perturbed.radius += 5e-5;
            auto candidate = RecoverApproxPatch(*s, perturbed);
            ValidationOptions rejectBudget;
            rejectBudget.maxValidationPieces = 1000;
            auto before = ValidateConversion(*s, candidate, rejectBudget);
            Require(!before.Passed(), "Perturbation unexpectedly passes");
            auto cache = SampleSurface(*s, {});
            ConversionStatistics stats;
            RefineCandidate(*s, candidate, cache, {}, stats);
            auto after = ValidateConversion(*s, candidate);
            std::cout << "refine " << cases[k].name
                      << " iterations=" << stats.optimizationIterations
                      << " passed=" << after.Passed() << " bounds=" << after.geometry.upperBound
                      << "," << after.boundaryTangent.upperBound << ","
                      << after.boundaryNormal.upperBound << std::endl;
            Require(stats.optimizationIterations > 0, "Optimizer not exercised");
            Require(after.Passed(), cases[k].name + " perturbed candidate not recovered");
            s->ExchangeUV();
            s->UReverse();
            s->VReverse();
            s->InsertUKnots(s->MinParamU() + (s->MaxParamU() - s->MinParamU()) * .371);
            auto transformed = ConvertApproximate(*s);
            Timing(transformed);
            Print(cases[k].name + " transformed", transformed);
            Require(transformed.status == ConversionStatus::Success,
                    cases[k].name + " reversed/swapped/inserted failed");
            Require(transformed.patch->parameters.type == cases[k].analytic->SurfaceType(),
                    "Transformed type mismatch");
        }
        for (unsigned k : {5u, 6u, 11u, 12u, 13u})
        {
            auto s = tests::RefitPolynomial(*cases[k].spline, 16);
            auto r = ConvertApproximate(*s);
            Timing(r);
            Print(cases[k].name, r);
            Require(r.status == ConversionStatus::Success,
                    cases[k].name + " extended chart failed");
            Require(r.patch->parameters.type == cases[k].analytic->SurfaceType(),
                    "Extended chart type mismatch");
        }
        // Known feasible plane at half the interior plateau height.

        auto bump = tests::PlateauFixture(*cases[0].spline);
        auto exact = Recognize(*cases[0].spline);
        auto z = Get(exact.parameters.axis);
        auto feasible = RecoverApproxPatch(*bump, exact.parameters);
        feasible.parameters.origin = Point(Get(feasible.parameters.origin) + z * 8e-6);
        MakeNative(feasible);
        Require(ValidateConversion(*bump, feasible).Passed(),
                "Known feasible plateau fixture invalid");
        auto optimized = ConvertApproximate(*bump);
        Timing(optimized);
        Print("plateau", optimized);
        Require(optimized.status == ConversionStatus::Success, "Plateau conversion failed");
        Require(optimized.timing.optimization.calls > 0 &&
                    optimized.statistics.optimizationIterations > 0,
                "Full conversion did not time optimizer");
        Require(optimized.timing.validation.calls >= 2, "Validation retry not timed");
        ConversionOptions disabled;
        disabled.refinement.maxIterations = 0;
        auto failed = ConvertApproximate(*bump, disabled);
        Timing(failed);
        Require(failed.status != ConversionStatus::Success, "Unoptimized plateau falsely accepted");
        Require(failed.timing.optimization.calls > 0 && failed.timing.validation.calls >= 2,
                "Failed attempt timings lost");
        auto rational = ConvertApproximate(*cases[1].spline);
        Timing(rational);
        Require(rational.status == ConversionStatus::Unsupported,
                "Rational input accepted by approximate API");
        std::cout << "Approximate assertions: " << checks << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << std::endl;
        code = 1;
    }
    sggk::fini();
    return code;
}
