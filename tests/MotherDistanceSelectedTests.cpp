#include "Dataset500.h"
#include "MotherDistanceInternal.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace
{
void Require(bool ok, const std::string &message)
{
    if (!ok)
        throw std::runtime_error(message);
}
std::vector<RecoveredCandidate> Classify(paired::AnalysisContext &context)
{
    context.classification.emplace();
    std::vector<RecoveredCandidate> recovered;
    ClassifyApproximate(context.samples.source, {}, nullptr, true, &context.samples, &recovered,
                        &*context.classification);
    Require(!recovered.empty(), "No recovered candidate");
    return recovered;
}
void Audit(const DatasetCase &c, const AnalyticPatch &p, const MotherDistanceValidationReport &r)
{
    Require(r.Passed(), r.message);
    const auto &t = r.timing;
    Require(std::abs(t.preparation.elapsedMs + t.surfaceBounds.elapsedMs + t.subdivision.elapsedMs +
                     t.boundaryChecks.elapsedMs + t.reporting.elapsedMs - t.totalMs) <
                1e-7 * std::max(1., t.totalMs),
            "Timer mismatch");
    SurfaceSampleContext samples(*c.input);
    for (unsigned i = 0; i <= 16; ++i)
        for (unsigned j = 0; j <= 16; ++j)
        {
            double u =
                c.input->MinParamU() + (c.input->MaxParamU() - c.input->MinParamU()) * i / 16;
            double v =
                c.input->MinParamV() + (c.input->MaxParamV() - c.input->MinParamV()) * j / 16;
            double d = Distance(
                p.parameters,
                Point(samples.At(u, v, sggk::DerivType::NONE, sggk::DerivType::NONE, false).point));
            Require(d < 1e-5 && d <= r.motherDistance.upperBound + 1e-11,
                    "Independent point audit exceeds bound");
        }
}
#ifdef BSPLINE_MOTHER_STREAMING
void ExportExtraction(const mother::SourceCache &s, const std::filesystem::path &file)
{
    std::ofstream o(file);
    Require(bool(o), "Cannot open extraction audit");
    o << std::setprecision(17);
    auto vector = [&](V p) { o << '[' << p.x << ',' << p.y << ',' << p.z << ']'; };
    auto array = [&](const auto &values) {
        o << '[';
        bool first = true;
        for (auto v : values)
        {
            if (!first)
                o << ',';
            first = false;
            o << v;
        }
        o << ']';
    };
    o << "{\"p\":" << s.p << ",\"q\":" << s.q << ",\"origin\":";
    vector(s.origin);
    o << ",\"knots_u\":";
    array(s.knotsU);
    o << ",\"knots_v\":";
    array(s.knotsV);
    o << ",\"mults_u\":";
    array(s.multsU);
    o << ",\"mults_v\":";
    array(s.multsV);
    o << ",\"controls\":[";
    for (size_t i = 0; i < s.controls.size(); ++i)
    {
        if (i)
            o << ',';
        o << '[';
        for (size_t j = 0; j < s.controls[i].size(); ++j)
        {
            if (j)
                o << ',';
            vector(s.controls[i][j]);
        }
        o << ']';
    }
    o << "],\"cells\":[";
    for (size_t id = 0; id < s.cells.size(); ++id)
    {
        if (id)
            o << ',';
        Require(bool(s.cells[id]), "Audit encountered missing source cell");
        const auto &c = *s.cells[id];
        o << "{\"first_u\":" << s.us[id / s.vs.size()].first
          << ",\"first_v\":" << s.vs[id % s.vs.size()].first << ",\"domain\":";
        array(c.domain);
        o << ",\"error\":";
        vector(c.error);
        o << ",\"coefficients\":[";
        for (size_t k = 0; k < c.x.c.size(); ++k)
        {
            if (k)
                o << ',';
            vector({c.x.c[k], c.y.c[k], c.z.c[k]});
        }
        o << "]}";
    }
    o << "]}";
}
#endif
void Run(const std::filesystem::path &dir, unsigned id, const std::filesystem::path &file,
         unsigned runs)
{
    auto c = ReadDatasetCase(dir / CaseFilename(id));
    CheckDataset500Contract(c, id);
    auto certificate = BoundPair(*c.source, *c.input, c.refitted ? 1e-6 : 1e-9);
    Require(certificate.passed, "Input source certificate failed");
    // Process/SDK and code warmup; each measured cold run still has a fresh cache.
    {
        paired::AnalysisContext ctx(*c.input);
        auto candidates = Classify(ctx);
        Audit(c, candidates.front().patch, mother::ValidatePrepared(ctx, candidates.front()));
    }
    std::ofstream csv(file);
    Require(bool(csv), "Cannot open selected CSV");
    csv << std::setprecision(17);
    csv << "id,type,run,cache,candidates,status,classification_ms,validation_ms,algorithm_ms,"
           "preparation_ms,position_ms,subdivision_ms,boundary_ms,reporting_ms,source_extraction_"
           "ms,source_boundary_ms,target_preparation_ms,sample_ms,residual_build_ms,residual_bound_"
           "ms,boundary_position_ms,boundary_angle_ms,surface_nodes,boundary_nodes,mapping_"
           "attempts,source_cache_hits,mother_bound,edge_bound,tangent_u_bound_deg,tangent_v_bound_deg,"
           "joint_tangent_bound_deg,joint_sine_squared_bound,joint_tangent_failure,source_"
           "certificate_bound,source_cells,extracted_cells,cached_cell_hits,tightened_cells,"
           "boundary_fallbacks,extraction_plan_ms,source_contraction_ms,source_envelope_ms,"
           "extraction_error\n";
    for (unsigned run = 0; run < runs; ++run)
    {
        paired::AnalysisContext ctx(*c.input);
        auto start = std::chrono::steady_clock::now();
        auto candidates = Classify(ctx);
        double classification =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        for (unsigned warm = 0; warm < 2; ++warm)
        {
            auto r = mother::ValidatePrepared(ctx, candidates.front());
            const auto &t = r.timing;
            Audit(c, candidates.front().patch, r);
#ifdef BSPLINE_MOTHER_STREAMING
            if (run == 0 && warm == 0)
                ExportExtraction(*ctx.motherSource, file.string() + ".extraction.json");
#endif
            Require(r.sourceCacheHits == warm && r.sourceExtractions == 1 - warm,
                    "Unexpected cache mode");
            csv << id << ',' << TypeName(c.truth.type) << ',' << run << ','
                << (warm ? "warm" : "cold") << ',' << candidates.size() << ','
                << PairedValidationStatusName(r.status) << ',';
            if (!warm)
                csv << classification;
            csv << ',' << t.totalMs << ',';
            if (!warm)
                csv << classification + t.totalMs;
            csv << ',' << t.preparation.elapsedMs << ',' << t.surfaceBounds.elapsedMs << ','
                << t.subdivision.elapsedMs << ',' << t.boundaryChecks.elapsedMs << ','
                << t.reporting.elapsedMs << ',' << r.sourceExtractionMs << ',' << r.sourceBoundaryMs
                << ',' << r.targetPreparationMs << ',' << r.sampleMs << ',' << r.residualBuildMs
                << ',' << r.residualBoundMs << ',' << r.boundaryPositionMs << ','
                << r.boundaryAngleMs << ',' << r.surfaceNodes << ',' << r.boundaryNodes << ','
                << r.mappingAttempts << ',' << r.sourceCacheHits << ','
                << r.motherDistance.upperBound << ',' << r.boundaryPosition.upperBound << ','
                << r.boundaryTangentU.upperBound << ',' << r.boundaryTangentV.upperBound << ','
                << r.boundaryJointTangent.upperBound << ',' << r.jointTangentSineSquaredUpperBound
                << ',' << r.jointTangentFailure << ',' << certificate.upper;
#ifdef BSPLINE_MOTHER_STREAMING
            csv << ',' << r.sourceCells << ',' << r.extractedCells << ',' << r.cachedCellHits << ','
                << r.tightenedCells << ',' << r.boundaryFallbacks << ',' << r.extractionPlanMs
                << ',' << r.sourceContractionMs << ',' << r.sourceEnvelopeMs << ','
                << r.maxExtractionError;
#else
            csv << ",,,,,,,,,";
#endif
            csv << '\n';
            csv.flush();
        }
    }
    std::cout << "selected " << id << " runs=" << runs << " PASSED\n";
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init();
    int result = 0;
    try
    {
        Require(argc == 5, "Usage: dataset id csv runs");
        Run(argv[1], static_cast<unsigned>(std::stoul(argv[2])), argv[3],
            static_cast<unsigned>(std::stoul(argv[4])));
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << '\n';
        result = 1;
    }
    sggk::fini();
    return result;
}
