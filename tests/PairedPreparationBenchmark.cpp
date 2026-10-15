#include "Dataset500.h"
#include "PairedValidationInternal.h"
#include <Foundation/init.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace bspline2analytic
{
// Linked only by BenchmarkPairedPreparation.py's isolated comparison build.
PairedValidationReport ValidatePairedConversionBaseline(const sggk::BSplineSurface &,
                                                        const AnalyticPatch &,
                                                        const PairedValidationOptions &);
std::shared_ptr<void> MakePairedBaselineCache();
PairedValidationReport ValidatePairedConversionBaselineCached(const sggk::BSplineSurface &,
                                                              const AnalyticPatch &,
                                                              const PairedValidationOptions &,
                                                              void *);
} // namespace bspline2analytic
using namespace bspline2analytic;
using namespace bspline2analytic::tests;
using namespace bspline2analytic::detail;
namespace
{
void Require(bool b, const char *message)
{
    if (!b)
        throw std::runtime_error(message);
}
void Check(const PairedValidationReport &r)
{
    Require(r.Passed() && r.certificate.geometry.upperBound < 1e-5 &&
                r.certificate.boundaryTangent.upperBound < .1 &&
                r.certificate.boundaryNormal.upperBound < .1,
            "Paired certificate failed");
    const auto &t = r.certificate.timing;
    double sum = t.preparation.elapsedMs + t.surfaceBounds.elapsedMs + t.subdivision.elapsedMs +
                 t.boundaryChecks.elapsedMs + t.reporting.elapsedMs;
    Require(std::abs(sum - t.totalMs) < 1e-7 * std::max(1.0, t.totalMs), "Timing mismatch");
}
void Same(const PairedValidationReport &a, const PairedValidationReport &b)
{
    Require(
        a.certificate.geometry.upperBound == b.certificate.geometry.upperBound &&
            a.certificate.boundaryTangent.upperBound == b.certificate.boundaryTangent.upperBound &&
            a.certificate.boundaryNormal.upperBound == b.certificate.boundaryNormal.upperBound &&
            a.commonCells == b.commonCells && a.surfaceNodes == b.surfaceNodes &&
            a.boundaryNodes == b.boundaryNodes && a.uMap.source == b.uMap.source &&
            a.uMap.target == b.uMap.target && a.vMap.source == b.vMap.source &&
            a.vMap.target == b.vMap.target,
        "Cached/uncached or repeated results differ");
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init();
    int code = 0;
    try
    {
        Require(argc == 4, "Expected data directory, output CSV, repetitions");
        unsigned repetitions = static_cast<unsigned>(std::stoul(argv[3]));
        Require(repetitions > 0 && repetitions <= 20, "Invalid repetition count");
        std::ofstream csv(argv[2]);
        Require(bool(csv), "Cannot open comparison CSV");
        csv << std::setprecision(17);
        csv << "id,run,controls,degree_u,degree_v,cells,cache_fill_ms,baseline_cache_fill_ms";
        for (auto name : {"baseline", "new", "cached", "baseline_cached"})
            for (auto field : {"wall_ms",
                               "prepare_ms",
                               "extract_ms",
                               "local_ms",
                               "cache_ms",
                               "surface_ms",
                               "boundary_ms",
                               "extractions",
                               "cache_hits",
                               "u_matrices",
                               "v_matrices",
                               "geometry",
                               "tangent",
                               "normal",
                               "compose_ms",
                               "product_ms",
                               "difference_ms",
                               "fast_bound_ms",
                               "tight_bound_ms",
                               "fast_passes",
                               "tight_checks",
                               "source_center_ms",
                               "source_boundary_ms",
                               "target_transform_ms",
                               "target_boundary_ms",
                               "edge_compose_ms",
                               "edge_algebra_ms",
                               "edge_bounds_ms",
                               "edge_subdivision_ms",
                               "source_boundary_builds",
                               "boundary_segments"})
                csv << ',' << name << '_' << field;
        csv << '\n';
        unsigned cases = 0;
        // Each row compares all four modes on the same actual initial candidate.
        for (unsigned id = 301; id < 400; id += 2)
        {
            auto c = ReadDatasetCase(std::filesystem::path(argv[1]) / CaseFilename(id));
            if (c.collapsedMask)
                continue;
            Require(!c.input->IsRational(), "Wrong input representation");
            Require(BoundPair(*c.source, *c.input, 1e-6).passed, "Fixture certificate failed");
            auto candidates = ClassifyApproximate(*c.input, {});
            Require(!candidates.empty(), "No actual candidate");
            const auto &p = candidates.front();
            PairedValidationContext cache;
            auto baselineCache = MakePairedBaselineCache();
            std::array<PairedValidationReport, 4> reports, reference;
            std::array<double, 4> wall{};
            auto run = [&](unsigned mode) {
                auto start = std::chrono::steady_clock::now();
                reports[mode] = mode == 0   ? ValidatePairedConversionBaseline(*c.input, p, {})
                                : mode == 1 ? ValidatePairedConversion(*c.input, p)
                                : mode == 2 ? ValidatePairedConversion(*c.input, p, cache)
                                            : ValidatePairedConversionBaselineCached(
                                                  *c.input, p, {}, baselineCache.get());
                wall[mode] = std::chrono::duration<double, std::milli>(
                                 std::chrono::steady_clock::now() - start)
                                 .count();
                Check(reports[mode]);
            };
            for (unsigned mode = 0; mode < 4; ++mode)
                run(mode);
            Require(reports[2].sourceExtractions == 1 && !reports[2].sourceCacheHits,
                    "Initial cache fill was not measured");
            const double fillMs = wall[2], baselineFillMs = wall[3];
            reference = reports;
            Same(reference[1], reference[2]);
            constexpr unsigned order[8][4] = {{0, 1, 2, 3}, {1, 2, 3, 0}, {2, 3, 0, 1},
                                              {3, 0, 1, 2}, {3, 2, 1, 0}, {2, 1, 0, 3},
                                              {1, 0, 3, 2}, {0, 3, 2, 1}};
            for (unsigned repeat = 0; repeat < repetitions; ++repeat)
            {
                for (unsigned mode : order[(cases + repeat) % 8])
                    run(mode);
                for (unsigned mode = 0; mode < 4; ++mode)
                    Same(reports[mode], reference[mode]);
                Same(reports[1], reports[2]);
                Same(reports[0], reports[3]);
                Require(reports[3].sourceCacheHits == 1 && reports[3].sourceExtractions == 0,
                        "Baseline cache was not exercised");
                Require(reports[1].sourceExtractions == 1 && reports[2].sourceExtractions == 0 &&
                            reports[2].sourceCacheHits == 1 && reports[2].sourceExtractionMs == 0,
                        "Invalid cache behavior");
                csv << id << ',' << repeat + 1 << ','
                    << c.input->ControlPoints().size() * c.input->ControlPoints()[0].size() << ','
                    << c.input->DegreeU() << ',' << c.input->DegreeV() << ','
                    << reports[1].commonCells << ',' << fillMs << ',' << baselineFillMs;
                for (unsigned mode = 0; mode < 4; ++mode)
                {
                    const auto &r = reports[mode];
                    const auto &t = r.certificate.timing;
                    csv << ',' << wall[mode] << ',' << t.preparation.elapsedMs << ','
                        << r.sourceExtractionMs << ',' << r.localTransformMs << ','
                        << r.sourceCacheMs << ',' << t.surfaceBounds.elapsedMs << ','
                        << t.boundaryChecks.elapsedMs << ',' << r.sourceExtractions << ','
                        << r.sourceCacheHits << ',' << r.uExtractionMatrices << ','
                        << r.vExtractionMatrices << ',' << r.certificate.geometry.upperBound << ','
                        << r.certificate.boundaryTangent.upperBound << ','
                        << r.certificate.boundaryNormal.upperBound << ',' << r.surfaceComposeMs
                        << ',' << r.surfaceProductMs << ',' << r.surfaceDifferenceMs << ','
                        << r.surfaceFastBoundMs << ',' << r.surfaceTightBoundMs << ','
                        << r.fastBoundPasses << ',' << r.tightBoundChecks << ',' << r.sourceCenterMs
                        << ',' << r.sourceBoundaryMs << ',' << r.targetTransformMs << ','
                        << r.targetBoundaryMs << ',' << r.edgeComposeMs << ',' << r.edgeAlgebraMs
                        << ',' << r.edgeBoundsMs << ',' << r.edgeSubdivisionMs << ','
                        << r.sourceBoundaryBuilds << ',' << r.boundarySegments;
                }
                csv << '\n';
            }
            ++cases;
        }
        Require(cases == 40, "Incomplete sphere coverage");
        std::cout << "PASS: " << cases << " cases, " << repetitions << " repeats, 4 modes"
                  << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << std::endl;
        code = 1;
    }
    sggk::fini();
    return code;
}
