#include "GeneralDataset.h"
#include "GeneralInternal.h"
#include <Foundation/init.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace bspline2analytic::detail
{
// The runner extracts and renames the baseline implementation from Git.
std::vector<GeneralSeed> ClassifyGeneralBeforeDedup(SurfaceSampleContext &, const GeneralConversionOptions &);
}
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
static std::vector<double> Values(const std::vector<GeneralSeed> &seeds)
{
    std::vector<double> values;
    for (const auto &seed : seeds)
    {
        const auto &p = seed.parameters;
        const auto &s = seed.score;
        values.push_back(double(p.type));
        for (V v : {Get(p.origin), Get(p.axis), Get(p.xDirection)})
            for (double x : {v.x, v.y, v.z}) values.push_back(x);
        for (double x : {p.radius, p.majorRadius, p.halfAngle, double(s.type),
                        double(s.degreesOfFreedom), s.interiorMax, s.interiorRms,
                        s.boundaryMax, s.normalMax}) values.push_back(x);
    }
    return values;
}
int main(int argc, char **argv)
{
    std::cerr << "Retired analytic-UV fallback dataset entry. Use GeneralRationalDatasetTests and rational-domain-v2.\n";
    return 1;

    sggk::init();
    int code = 1;
    try
    {
        if (argc != 3) throw std::runtime_error("Expected dataset directory and output CSV");
        std::vector<GeneralFixture> fixtures;
        for (unsigned id = 0; id < 500; ++id)
        {
            auto fixture = ReadGeneralFixture(std::filesystem::path(argv[1]) / GeneralFixtureFilename(id));
            if (fixture.id != id || fixture.group != GeneralFixtureGroup::Feasible)
                throw std::runtime_error("Expected main positive fixture");
            fixtures.push_back(std::move(fixture));
        }
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("Cannot open benchmark CSV");
        out << std::setprecision(17)
            << "repeat,id,source_type,old_first,before_ms,after_ms,candidates,raw_candidates,duplicate_candidates,score_evaluations,identical\n";
        GeneralConversionOptions options;
        // Keep the historical budget/mode explicit. The strict comparison also
        // reports subsequent algorithm changes beyond the old dedup baseline.
        options.seedAttemptsPerType = 64;
        options.coarseScorePointCount = 0;
        for (unsigned repeat = 0; repeat < 3; ++repeat)
        {
            double beforeTotal = 0, afterTotal = 0;
            for (const auto &fixture : fixtures)
            {
                std::vector<GeneralSeed> before, after;
                GeneralClassificationDiagnostics diagnostic;
                double beforeMs = 0, afterMs = 0;
                auto run = [&](bool old) {
                    SurfaceSampleContext context(*fixture.input, options.fitting.sampling.maxSamples);
                    auto start = std::chrono::steady_clock::now();
                    auto seeds = old ? ClassifyGeneralBeforeDedup(context, options)
                                     : ClassifyGeneral(context, options, &diagnostic);
                    double ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - start).count();
                    if (old) { before = std::move(seeds); beforeMs = ms; }
                    else { after = std::move(seeds); afterMs = ms; }
                };
                bool oldFirst = (fixture.id + repeat) % 2 == 0;
                run(oldFirst); run(!oldFirst);
                if (Values(before) != Values(after))
                    throw std::runtime_error("Ranked candidates differ at " + std::to_string(fixture.id));
                beforeTotal += beforeMs; afterTotal += afterMs;
                out << repeat << ',' << fixture.id << ',' << fixture.type << ',' << oldFirst << ','
                    << beforeMs << ',' << afterMs << ',' << after.size() << ','
                    << diagnostic.rawCandidateCount << ',' << diagnostic.duplicateCandidateCount << ','
                    << diagnostic.scoreEvaluationCount << ",1\n";
            }
            out.flush();
            std::cout << "repeat=" << repeat << " before_ms=" << beforeTotal
                      << " after_ms=" << afterTotal << "; all 500 ranked candidate lists identical" << std::endl;
        }
        if (!out) throw std::runtime_error("Benchmark CSV write failed");
        code = 0;
    }
    catch (const std::exception &e) { std::cerr << e.what() << std::endl; }
    sggk::fini();
    return code;
}
