#include "ApproximateCases.h"
#include <Foundation/init.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
using namespace bspline2analytic;
static void Summary(const std::string &name, std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    auto quantile = [&](double q) {
        return values[static_cast<size_t>(std::ceil(q * values.size())) - 1];
    };
    std::cout << name
              << " mean=" << std::accumulate(values.begin(), values.end(), 0.0) / values.size()
              << " p50=" << quantile(.5) << " p95=" << quantile(.95) << " max=" << values.back()
              << " ms\n";
}
static int Run(const char *path)
{
    auto references = tests::MakeCases();
    struct Case
    {
        std::string name;
        sggk::BSplineSurfacePtr surface;
        ConversionOptions options;
        bool success;
    };
    std::vector<Case> cases;
    for (unsigned k = 0; k < 5; ++k)
        cases.push_back(
            {references[k].name, tests::RefitPolynomial(*references[k].spline, 12), {}, true});
    auto bump = tests::PlateauFixture(*references[0].spline);
    cases.push_back({"plateau", bump, {}, true});
    ConversionOptions disabled;
    disabled.refinement.maxIterations = 0;
    cases.push_back({"plateau-disabled", bump, disabled, false});
    std::ofstream csv(path);
    if (!csv)
        throw std::runtime_error("Cannot open benchmark CSV");
    csv << "case,run,controls,status,classify_ms,optimize_ms,validate_ms,other_ms,total_ms,"
           "classify_calls,optimize_calls,validate_calls,candidates,iterations,rounds,samples,"
           "pieces,distance_bound,tangent_deg_bound,normal_deg_bound\n";
    csv << std::setprecision(17);
    std::cout << std::fixed << std::setprecision(4);
    for (const auto &c : cases)
    {
        // SDK startup, fixture generation, warm-up, console/CSV I/O are outside timers.
        auto warm = ConvertApproximate(*c.surface, c.options);
        if ((warm.status == ConversionStatus::Success) != c.success)
            throw std::runtime_error(c.name + " warm-up status mismatch");
        std::array<std::vector<double>, 4> stages;
        for (unsigned run = 0; run < 5; ++run)
        {
            auto r = ConvertApproximate(*c.surface, c.options);
            if ((r.status == ConversionStatus::Success) != c.success)
                throw std::runtime_error(c.name + " benchmark status mismatch");
            const auto &t = r.timing;
            const auto &s = r.statistics;
            stages[0].push_back(t.classification.elapsedMs);
            stages[1].push_back(t.optimization.elapsedMs);
            stages[2].push_back(t.validation.elapsedMs);
            stages[3].push_back(t.totalMs);
            csv << c.name << ',' << run << ','
                << c.surface->ControlPoints().size() * c.surface->ControlPoints().front().size()
                << ',' << ConversionStatusName(r.status) << ',' << t.classification.elapsedMs << ','
                << t.optimization.elapsedMs << ',' << t.validation.elapsedMs << ',' << t.otherMs
                << ',' << t.totalMs << ',' << t.classification.calls << ',' << t.optimization.calls
                << ',' << t.validation.calls << ',' << s.candidateCount << ','
                << s.optimizationIterations << ',' << s.refinementRounds << ',' << s.sampleCount
                << ',' << s.validationPieces << ',' << r.validation.geometry.upperBound << ','
                << r.validation.boundaryTangent.upperBound << ','
                << r.validation.boundaryNormal.upperBound << '\n';
        }
        for (unsigned i = 0; i < 4; ++i)
            Summary(c.name + " " +
                        std::array<const char *, 4>{"classification", "optimization", "validation",
                                                    "total"}[i],
                    stages[i]);
    }
    if (!csv)
        throw std::runtime_error("Benchmark CSV write failed");
    std::cout << "35 measured conversions; 30 success, 5 intentionally unsuccessful.\n";
    return 0;
}
int main(int argc, char **argv)
{
    sggk::init();
    int result = 1;
    try
    {
        result = Run(argc > 1 ? argv[1] : "approximate-benchmark.csv");
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
    sggk::fini();
    return result;
}
