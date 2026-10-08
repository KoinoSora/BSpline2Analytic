#include "AnalyticReconstruction.h"
#include "Numerics.h"
#include "ReferenceCases.h"
#include <Foundation/init.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

using namespace bspline2analytic;
static unsigned checks = 0;
static void Require(bool condition, const std::string &message)
{
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}
static void Check(const bspline2analytic::tests::ReferenceCase &reference,
                  const sggk::BSplineSurface &s, const std::string &suffix)
{
    auto r = Recognize(s);
    Require(r.status == Status::Success,
            reference.name + suffix + ": " + StatusName(r.status) + " " + r.message);
    Require(r.parameters.type == reference.analytic->SurfaceType(),
            reference.name + suffix + ": wrong analytic type");
    Require(r.domain.has_value() && r.surface != nullptr, "Missing output geometry");
    Require(r.validation.level == VerificationLevel::BernsteinBound &&
                r.validation.surfaceDistanceBound <= 1e-5,
            "Missing distance bound");
    const std::array<sggk::Point3D, 4> corners{
        s.CalcPoint(s.MinParamU(), s.MinParamV()), s.CalcPoint(s.MaxParamU(), s.MinParamV()),
        s.CalcPoint(s.MinParamU(), s.MaxParamV()), s.CalcPoint(s.MaxParamU(), s.MaxParamV())};
    for (unsigned i = 0; i < 4; ++i)
    {
        const auto &uv = r.mapping.analyticCorners[i];
        Require(detail::Norm(detail::Get(corners[i]) - detail::Get(r.surface->CalcPoint(uv))) <=
                    1e-5,
                reference.name + suffix + ": recovered corner mismatch");
        Require(uv.X() >= r.domain->IntervalU().Min() - 1e-10 &&
                    uv.X() <= r.domain->IntervalU().Max() + 1e-10 &&
                    uv.Y() >= r.domain->IntervalV().Min() - 1e-10 &&
                    uv.Y() <= r.domain->IntervalV().Max() + 1e-10,
                reference.name + suffix + ": corner outside recovered range");
    }
    auto coordinate = [&](unsigned corner, unsigned dimension) {
        return dimension ? r.mapping.analyticCorners[corner].Y()
                         : r.mapping.analyticCorners[corner].X();
    };
    unsigned dimensionU = r.mapping.swappedUV ? 1 : 0, dimensionV = 1 - dimensionU;
    Require(r.mapping.reversedU == (coordinate(1, dimensionU) < coordinate(0, dimensionU)),
            "U orientation mismatch");
    Require(r.mapping.reversedV == (coordinate(2, dimensionV) < coordinate(0, dimensionV)),
            "V orientation mismatch");
    // Independent point projection through SGK tests the returned SDK object's frame and
    // parameters.
    for (unsigned i = 0; i <= 6; ++i)
        for (unsigned j = 0; j <= 6; ++j)
        {
            auto p = s.CalcPoint(s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 6,
                                 s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * j / 6);
            auto uv = r.surface->CalcParam(p);
            auto q = r.surface->CalcPoint(uv);
            Require(detail::Norm(detail::Get(p) - detail::Get(q)) <= 1e-5,
                    reference.name + suffix + ": SGK reconstructed surface mismatch");
        }
}
static int Run(bool benchmark, const std::string &csvPath)
{
    auto references = bspline2analytic::tests::MakeCases();
    if (benchmark)
    {
        std::vector<sggk::BSplineSurfacePtr> data;
        for (unsigned variant = 0; variant < 4; ++variant)
            for (unsigned n = 0; n < 5; ++n)
            {
                auto s = std::make_shared<sggk::BSplineSurface>(*references[n].spline);
                const std::array<unsigned, 4> targetU{5, 8, 10, 10}, targetV{10, 10, 15, 20};
                unsigned insertionsU =
                    targetU[variant] - static_cast<unsigned>(s->ControlPoints().size());
                unsigned insertionsV =
                    targetV[variant] - static_cast<unsigned>(s->ControlPoints().front().size());
                for (unsigned k = 1; k <= insertionsU; ++k)
                {
                    // Avoid existing analytic-conversion knots (notably the midpoint): SGK
                    // InsertKnots specifies target multiplicity rather than an increment.
                    double fraction = double(k) / (insertionsU + std::sqrt(2.0));
                    s->InsertUKnots(s->MinParamU() + (s->MaxParamU() - s->MinParamU()) * fraction);
                }
                for (unsigned k = 1; k <= insertionsV; ++k)
                {
                    double fraction = double(k) / (insertionsV + std::sqrt(2.0));
                    s->InsertVKnots(s->MinParamV() + (s->MaxParamV() - s->MinParamV()) * fraction);
                }
                data.push_back(s);
            }
        std::ofstream csv;
        if (!csvPath.empty())
        {
            csv.open(csvPath);
            Require(bool(csv), "Cannot create benchmark CSV");
            csv << "index,type,controls,elapsed_ms,sampled_distance,bound\n"
                << std::setprecision(17);
        }
        std::array<std::vector<double>, 5> elapsed;
        std::array<std::size_t, 5> minControls{1000, 1000, 1000, 1000, 1000}, maxControls{};
        for (unsigned i = 0; i < 1000; ++i)
        {
            unsigned k = i % 5;
            auto s = std::make_shared<sggk::BSplineSurface>(*data[i % data.size()]);
            double scale = 0.5 + 0.013 * (i % 97);
            auto controls = s->ControlPoints();
            // Every input has distinct placement and scale; preparation is outside the timed call.
            for (unsigned u = 0; u < controls.size(); ++u)
                for (unsigned v = 0; v < controls[u].size(); ++v)
                {
                    const auto &p = controls[u][v];
                    s->SetControlPoint(sggk::Point3D(p.X() * scale + i * 0.31,
                                                     p.Y() * scale - i * 0.17,
                                                     p.Z() * scale + i * 0.07),
                                       u, v);
                }
            if (i % 2)
                s->ExchangeUV();
            if (i % 3 == 0)
                s->UReverse();
            if (i % 7 == 0)
                s->VReverse();
            const auto count = controls.size() * controls.front().size();
            const std::array<std::size_t, 4> expectedCounts{50, 80, 150, 200};
            Require(count == expectedCounts[(i % data.size()) / 5],
                    "Benchmark controls differ from target count: " + std::to_string(count));
            auto start = std::chrono::steady_clock::now();
            auto r = Recognize(*s);
            auto end = std::chrono::steady_clock::now();
            Require(r.status == Status::Success &&
                        r.parameters.type == references[k].analytic->SurfaceType(),
                    "Batch recognition failed: " + r.message);
            double ms = std::chrono::duration<double, std::milli>(end - start).count();
            elapsed[k].push_back(ms);
            minControls[k] = std::min(minControls[k], count);
            maxControls[k] = std::max(maxControls[k], count);
            if (csv)
                csv << i << ',' << references[k].name << ',' << count << ',' << ms << ','
                    << r.validation.maxSampledDistance << ',' << r.validation.surfaceDistanceBound
                    << '\n';
            if (i % 50 == 0)
            {
                auto repeat = Recognize(*s);
                Require(repeat.status == r.status && repeat.parameters.type == r.parameters.type &&
                            repeat.validation.surfaceDistanceBound ==
                                r.validation.surfaceDistanceBound,
                        "Repeated recognition changed result");
            }
        }
        for (unsigned k = 0; k < 5; ++k)
        {
            auto t = elapsed[k];
            std::sort(t.begin(), t.end());
            std::cout << references[k].name << ": controls=" << minControls[k] << '-'
                      << maxControls[k]
                      << " mean_ms=" << std::accumulate(t.begin(), t.end(), 0.0) / t.size()
                      << " p95_ms=" << t[t.size() * 95 / 100] << " max_ms=" << t.back() << '\n';
        }
        std::cout << "Batch passed: 1000 distinct patches, 20 deterministic repeat checks. This is "
                     "not a leak audit.\n";
        return 0;
    }
    unsigned failures = 0;
    for (const auto &reference : references)
    {
        try
        {
            Check(reference, *reference.spline, "");
            auto modified = std::make_shared<sggk::BSplineSurface>(*reference.spline);
            for (double t : {0.23, 0.51, 0.77})
            {
                modified->InsertUKnots(modified->MinParamU() +
                                       (modified->MaxParamU() - modified->MinParamU()) * t);
                modified->InsertVKnots(modified->MinParamV() +
                                       (modified->MaxParamV() - modified->MinParamV()) * t);
            }
            Check(reference, *modified, "/inserted");
            modified->ExchangeUV();
            modified->UReverse();
            modified->VReverse();
            Check(reference, *modified, "/swapped-reversed");
            auto elevated = bspline2analytic::tests::ElevateBezierPatch(*reference.spline);
            Check(reference, *elevated, "/elevated");
            std::cout << "PASS " << reference.name << '\n';
        }
        catch (const std::exception &e)
        {
            ++failures;
            std::cerr << "FAIL " << e.what() << '\n';
        }
    }
    auto changed = std::make_shared<sggk::BSplineSurface>(*references[3].spline);
    changed->InsertUKnots((changed->MinParamU() + changed->MaxParamU()) / 2);
    changed->InsertVKnots((changed->MinParamV() + changed->MaxParamV()) / 2);
    auto p = changed->ControlPoints()[1][1];
    p.SetZ(p.Z() + 0.1);
    changed->SetControlPoint(p, 1, 1);
    Require(AnalyzeBoundaries(*changed).status == Status::Success,
            "Interior perturbation changed boundaries unexpectedly");
    Require(Recognize(*changed).status != Status::Success,
            "Interior perturbation was not rejected");
    RecognitionOptions invalid;
    invalid.distanceTolerance = -1;
    Require(Recognize(*references[0].spline, invalid).status == Status::InvalidInput,
            "Invalid tolerance accepted");
    sggk::SphereSurface sphere(sggk::Ucs3D(), 3);
    auto pole = sphere.ToBSpline(sggk::UVRange(sggk::Interval(0.2, 1.2), sggk::Interval(0, 0.7)));
    Require(Recognize(*pole).status == Status::Degenerate,
            "Collapsed pole boundary must report Degenerate");
    auto ellipse = std::make_shared<sggk::BSplineSurface>(*references[1].spline);
    auto controls = ellipse->ControlPoints();
    for (unsigned i = 0; i < controls.size(); ++i)
        for (unsigned j = 0; j < controls[i].size(); ++j)
        {
            auto q = controls[i][j];
            q.SetX(q.X() * 2);
            ellipse->SetControlPoint(q, i, j);
        }
    Require(Recognize(*ellipse).status != Status::Success, "Elliptical boundaries accepted");
    RecognitionOptions sampled;
    sampled.requireBound = false;
    auto sampledResult = Recognize(*references[3].spline, sampled);
    Require(sampledResult.status == Status::Success &&
                sampledResult.validation.level == VerificationLevel::Sampled &&
                !std::isfinite(sampledResult.validation.surfaceDistanceBound),
            "Sample-only verification reported a global bound");
    RecognitionOptions limited;
    limited.maxValidationPieces = 1;
    Require(Recognize(*references[4].spline, limited).status == Status::VerificationFailed,
            "Validation budget exhaustion accepted");
    std::cout << checks << " assertions, " << failures << " failed reference cases\n";
    return failures ? 1 : 0;
}
int main(int argc, char **argv)
{
    sggk::init();
    int code = 1;
    try
    {
        code = Run(argc > 1 && std::string(argv[1]) == "--benchmark", argc > 2 ? argv[2] : "");
    }
    catch (const std::exception &e)
    {
        std::cerr << "TEST ERROR: " << e.what() << '\n';
    }
    sggk::fini();
    return code;
}
