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
struct Classified
{
    RecoveredCandidate candidate;
    MotherDistanceValidationReport report;
    ClassificationDiagnostics diagnostics;
    double classificationMs = 0;
    unsigned attempts = 0;
};
void Accumulate(MotherDistanceValidationReport &sum, const MotherDistanceValidationReport &r)
{
    sum.timing.totalMs += r.timing.totalMs;
    sum.timing.preparation.elapsedMs += r.timing.preparation.elapsedMs;
    sum.timing.surfaceBounds.elapsedMs += r.timing.surfaceBounds.elapsedMs;
    sum.timing.subdivision.elapsedMs += r.timing.subdivision.elapsedMs;
    sum.timing.boundaryChecks.elapsedMs += r.timing.boundaryChecks.elapsedMs;
    sum.timing.reporting.elapsedMs += r.timing.reporting.elapsedMs;
    sum.sourceExtractionMs += r.sourceExtractionMs;
    sum.sourceBoundaryMs += r.sourceBoundaryMs;
    sum.targetPreparationMs += r.targetPreparationMs;
    sum.mappingMs += r.mappingMs;
    sum.sampleMs += r.sampleMs;
    sum.residualBuildMs += r.residualBuildMs;
    sum.residualBoundMs += r.residualBoundMs;
    sum.boundaryPositionMs += r.boundaryPositionMs;
    sum.boundaryAngleMs += r.boundaryAngleMs;
    sum.extractionPlanMs += r.extractionPlanMs;
    sum.sourceContractionMs += r.sourceContractionMs;
    sum.sourceEnvelopeMs += r.sourceEnvelopeMs;
    sum.sourceExtractions += r.sourceExtractions;
    sum.sourceCacheHits += r.sourceCacheHits;
    sum.surfaceNodes += r.surfaceNodes;
    sum.boundaryNodes += r.boundaryNodes;
    sum.mappingAttempts += r.mappingAttempts;
    sum.extractedCells += r.extractedCells;
    sum.cachedCellHits += r.cachedCellHits;
    sum.tightenedCells += r.tightenedCells;
    sum.boundaryFallbacks += r.boundaryFallbacks;
    sum.tangentResidualFits += r.tangentResidualFits;
    sum.tangentResidualSplits += r.tangentResidualSplits;
}
Classified ClassifyAndValidate(paired::AnalysisContext &context, std::ostream *attempts = nullptr,
                              unsigned id = 0, unsigned run = 0)
{
    Classified out;
    MotherDistanceValidationReport sum;
    auto start = std::chrono::steady_clock::now();
    context.classification.emplace();
    auto search = std::make_unique<IsoClassificationSearch>(context.samples, ConversionOptions{}, out.diagnostics);
    out.classificationMs += std::chrono::duration<double,std::milli>(
        std::chrono::steady_clock::now()-start).count();
    for (;;)
    {
        start = std::chrono::steady_clock::now();
        auto candidate = search->Next(*context.classification);
        out.classificationMs += std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count();
        Require(bool(candidate), "Incremental classifier exhausted all candidates");
        ++out.attempts;
        auto r = mother::ValidatePrepared(context,*candidate);
        if (attempts)
        {
            *attempts << id << ',' << run << ',' << out.attempts << ','
                << PairedValidationStatusName(r.status) << ',' << out.classificationMs << ','
                << r.timing.totalMs << ',' << r.timing.preparation.elapsedMs << ','
                << r.timing.surfaceBounds.elapsedMs << ',' << r.timing.subdivision.elapsedMs << ','
                << r.timing.boundaryChecks.elapsedMs << ',' << r.timing.reporting.elapsedMs << ','
                << r.sourceCacheHits << ',' << r.motherDistance.upperBound << ','
                << r.boundaryPosition.upperBound << ',' << r.boundaryJointTangent.upperBound << ','
                << out.diagnostics.classificationStages << '\n';
            attempts->flush();
        }
        Accumulate(sum,r);
        if (!r.Passed()) continue;
        out.candidate = *candidate; out.report = r;
        // Geometry and certificates belong to the successful attempt; costs
        // and work counters include every failed attempt in this same context.
        out.report.timing = sum.timing;
        out.report.sourceExtractionMs = sum.sourceExtractionMs;
        out.report.sourceBoundaryMs = sum.sourceBoundaryMs;
        out.report.targetPreparationMs = sum.targetPreparationMs;
        out.report.mappingMs = sum.mappingMs;
        out.report.sampleMs = sum.sampleMs;
        out.report.residualBuildMs = sum.residualBuildMs;
        out.report.residualBoundMs = sum.residualBoundMs;
        out.report.boundaryPositionMs = sum.boundaryPositionMs;
        out.report.boundaryAngleMs = sum.boundaryAngleMs;
        out.report.extractionPlanMs = sum.extractionPlanMs;
        out.report.sourceContractionMs = sum.sourceContractionMs;
        out.report.sourceEnvelopeMs = sum.sourceEnvelopeMs;
        out.report.sourceExtractions = sum.sourceExtractions;
        out.report.sourceCacheHits = sum.sourceCacheHits;
        out.report.surfaceNodes = sum.surfaceNodes;
        out.report.boundaryNodes = sum.boundaryNodes;
        out.report.mappingAttempts = sum.mappingAttempts;
        out.report.extractedCells = sum.extractedCells;
        out.report.cachedCellHits = sum.cachedCellHits;
        out.report.tightenedCells = sum.tightenedCells;
        out.report.boundaryFallbacks = sum.boundaryFallbacks;
        out.report.tangentResidualFits = sum.tangentResidualFits;
        out.report.tangentResidualSplits = sum.tangentResidualSplits;
        start = std::chrono::steady_clock::now();
        search.reset();
        out.classificationMs += std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count();
        return out;
    }
}
struct AuditResult
{
    double position = 0, edge = 0, tangentU = 0, tangentV = 0, joint = 0;
};
double Map(const PairedParameterMap &m, double t)
{
    if (t == m.source.back()) return m.target.back();
    auto it = std::upper_bound(m.source.begin(), m.source.end(), t);
    Require(it != m.source.begin() && it != m.source.end(), "Audit map coverage");
    size_t i = static_cast<size_t>(it - m.source.begin() - 1);
    return m.target[i] + (t - m.source[i]) / (m.source[i + 1] - m.source[i]) *
        (m.target[i + 1] - m.target[i]);
}
AuditResult Audit(const DatasetCase &c, const AnalyticPatch &p, const MotherDistanceValidationReport &r)
{
    Require(r.Passed(), r.message);
    AuditResult out;
    const auto &timing = r.timing;
    Require(std::abs(timing.preparation.elapsedMs + timing.surfaceBounds.elapsedMs + timing.subdivision.elapsedMs +
                     timing.boundaryChecks.elapsedMs + timing.reporting.elapsedMs - timing.totalMs) <
                1e-7 * std::max(1., timing.totalMs),
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
            out.position = std::max(out.position, d);
        }
    // Independent analytic directions at the returned mapped points, rather
    // than reusing the validator's rational derivative numerators or products.
    const auto chart = paired::BuildChart(p, false);
    const auto frame = FrameOf(p.parameters);
    const auto &s = *c.input;
    for (unsigned edge = 0; edge < 4; ++edge)
        for (unsigned i = 0; i <= 16; ++i)
        {
            const double u = edge < 2 ? (edge ? s.MaxParamU() : s.MinParamU()) :
                s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 16;
            const double v = edge >= 2 ? (edge == 3 ? s.MaxParamV() : s.MinParamV()) :
                s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * i / 16;
            double q[2] = {Map(r.uMap, u), Map(r.vMap, v)}, angle[2]{};
            if (p.mapping.reversedU) q[0] = 1 - q[0];
            if (p.mapping.reversedV) q[1] = 1 - q[1];
            if (p.mapping.swappedUV) std::swap(q[0], q[1]);
            for (unsigned d = 0; d < 2; ++d)
            {
                unsigned n = static_cast<unsigned>(chart.angles[d].size() - 1);
                double x = std::clamp(q[d], 0., 1.) * n;
                unsigned k = std::min(n - 1, static_cast<unsigned>(x));
                double t = x - k, a = chart.angles[d][k], b = chart.angles[d][k + 1];
                angle[d] = chart.angular[d] ? a + 2 * std::atan2(t * std::sin((b - a) / 2),
                    1 - t + t * std::cos((b - a) / 2)) : a + t * (b - a);
            }
            const double cu = std::cos(angle[0]), su = std::sin(angle[0]);
            const double cv = std::cos(angle[1]), sv = std::sin(angle[1]);
            const auto &a = p.parameters;
            V du{}, dv{};
            switch (a.type)
            {
            case sggk::SurfaceType::Plane: du = {1, 0, 0}; dv = {0, 1, 0}; break;
            case sggk::SurfaceType::Cylinder:
                du = {-a.radius * su, a.radius * cu, 0}; dv = {0, 0, 1}; break;
            case sggk::SurfaceType::Cone: {
                double k = std::tan(a.halfAngle);
                du = {-angle[1] * k * su, angle[1] * k * cu, 0}; dv = {k * cu, k * su, 1}; break;
            }
            case sggk::SurfaceType::Sphere:
                du = {-a.radius * sv * su, a.radius * sv * cu, 0};
                dv = {a.radius * cv * cu, a.radius * cv * su, a.radius * sv}; break;
            case sggk::SurfaceType::Torus:
                du = {-(a.majorRadius + a.radius * cv) * su, (a.majorRadius + a.radius * cv) * cu, 0};
                dv = {-a.radius * sv * cu, -a.radius * sv * su, a.radius * cv}; break;
            default: throw std::runtime_error("Unsupported audit mother");
            }
            if (p.mapping.swappedUV) std::swap(du, dv);
            du = frame.vector(du) * (p.mapping.reversedU ? -1. : 1.);
            dv = frame.vector(dv) * (p.mapping.reversedV ? -1. : 1.);
            const auto &sample = samples.At(u, v);
            auto angleOf = [](V x, V y) {
                Require(Norm(x) > 0 && Norm(y) > 0, "Degenerate audit tangent");
                return std::atan2(Norm(Cross(x, y)), Dot(x, y));
            };
            double au = angleOf(sample.du, du), av = angleOf(sample.dv, dv);
            double joint = std::asin(std::min(1., std::sqrt(std::pow(std::sin(au), 2) +
                std::pow(std::sin(av), 2)))) * 180 / Pi;
            double distance = Norm(sample.point - frame.global(CanonicalPoint(a, angle[0], angle[1])));
            Require(std::isfinite(distance) && distance < 1e-5 &&
                distance <= r.boundaryPosition.upperBound + 1e-10, "Independent edge-position audit");
            Require(std::isfinite(au) && std::isfinite(av) && au < Pi / 2 && av < Pi / 2 &&
                au * 180 / Pi <= r.boundaryTangentU.upperBound + 1e-7 &&
                av * 180 / Pi <= r.boundaryTangentV.upperBound + 1e-7 &&
                joint < .1 && joint <= r.boundaryJointTangent.upperBound + 1e-7,
                "Independent joint U/V tangent audit");
            out.edge = std::max(out.edge, distance);
            out.tangentU = std::max(out.tangentU, au * 180 / Pi);
            out.tangentV = std::max(out.tangentV, av * 180 / Pi);
            out.joint = std::max(out.joint, joint);
        }
    return out;
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
// Exercise continuation after success using the same representative geometry.
// This is branch coverage, never included in production latency measurements.
void Coverage(const std::filesystem::path &dir, unsigned id, const std::filesystem::path &file)
{
    auto c = ReadDatasetCase(dir / CaseFilename(id));
    CheckDataset500Contract(c,id);
    Require(BoundPair(*c.source,*c.input,c.refitted?1e-6:1e-9).passed,"Input certificate failed");
    paired::AnalysisContext context(*c.input);
    context.classification.emplace();
    ClassificationDiagnostics diagnostics;
    IsoClassificationSearch search(context.samples,{},diagnostics);
    unsigned count=0,passed=0;
    std::ofstream csv(file);
    Require(bool(csv),"Cannot open continuation coverage CSV");
    csv << "id,attempt,stage,status,fit_passes,geometric_iterations,sampling_refinements\n";
    while (auto candidate=search.Next(*context.classification))
    {
        Require(++count<=512,"Unbounded candidate continuation");
        Require(candidate->sourceStamp==context.classification->stamp,"Stale sampling version");
        for (unsigned e=0;e<4;++e)
        {
            Require(candidate->trace.paths[e].size()==context.classification->boundaryParameters[e].size(),
                    "Trace/sample revision mismatch");
            for (size_t i=0;i<context.classification->boundaryParameters[e].size();++i)
            {
                double t=context.classification->boundaryParameters[e][i];
                double u=e<2?(e==0?c.input->MinParamU():c.input->MaxParamU()):t;
                double v=e<2?t:(e==2?c.input->MinParamV():c.input->MaxParamV());
                V direct=Get(c.input->CalcPoint(u,v)),cached=context.classification->boundaries[e][i];
                Require(Norm(direct-cached)<=512*std::numeric_limits<double>::epsilon()*(1+Norm(direct)),
                        "Isocurve/surface evaluation disagreement");
            }
        }
        auto report=mother::ValidatePrepared(context,*candidate);
        if (report.Passed()) { Audit(c,candidate->patch,report); ++passed; }
        csv << id << ',' << count << ',' << diagnostics.classificationStages << ','
            << PairedValidationStatusName(report.status) << ',' << diagnostics.fitPasses << ','
            << diagnostics.geometricIterations << ',' << diagnostics.samplingRefinements << '\n';
    }
    Require(passed>0 && diagnostics.classificationStages==6 && diagnostics.geometricIterations>0 &&
            diagnostics.samplingRefinements>0,"Missing refinement/continuation coverage");
    std::cout << "continuation " << id << " candidates=" << count << " certified=" << passed
              << " iterations=" << diagnostics.geometricIterations << " refinements="
              << diagnostics.samplingRefinements << " PASSED\n";
}
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
        auto classified = ClassifyAndValidate(ctx);
        Audit(c, classified.candidate.patch, classified.report);
    }
    std::ofstream attempts(file.string()+".attempts.csv");
    Require(bool(attempts), "Cannot open candidate attempts CSV");
    attempts << std::setprecision(17) << "id,run,attempt,status,classification_cumulative_ms,validation_ms,preparation_ms,position_ms,subdivision_ms,boundary_ms,reporting_ms,source_cache_hits,mother_bound,edge_bound,joint_bound_deg,classification_stages\n";
    std::ofstream csv(file);
    Require(bool(csv), "Cannot open selected CSV");
    csv << std::setprecision(17);
    csv << "id,type,run,cache,candidates,status,classification_ms,validation_ms,algorithm_ms,"
           "preparation_ms,position_ms,subdivision_ms,boundary_ms,reporting_ms,source_extraction_"
           "ms,source_boundary_ms,target_preparation_ms,mapping_ms,sample_ms,residual_build_ms,residual_bound_"
           "ms,boundary_position_ms,boundary_angle_ms,surface_nodes,boundary_nodes,mapping_"
           "attempts,source_cache_hits,mother_bound,edge_bound,tangent_u_bound_deg,tangent_v_bound_deg,"
           "joint_tangent_bound_deg,joint_sine_squared_bound,joint_tangent_failure,source_"
           "certificate_bound,source_cells,extracted_cells,cached_cell_hits,tightened_cells,"
           "boundary_fallbacks,extraction_plan_ms,source_contraction_ms,source_envelope_ms,"
           "extraction_error,audit_position_max,audit_edge_max,audit_tangent_u_max_deg,"
           "audit_tangent_v_max_deg,audit_joint_max_deg,degree_u,degree_v,controls,"
           "tangent_residual_fits,tangent_residual_splits,validation_attempts,class_sampling_ms,class_algebraic_ms,class_geometric_ms,class_scoring_ms,class_recovery_ms,class_other_ms,curve_extractions,curve_point_evaluations,score_point_evaluations,point_cache_hits,fit_passes,geometric_iterations,fast_circle_fits,sampling_refinements,classification_stages,raw_candidates,distance_evaluations,shared_mother_scores,duplicate_candidates,patch_recovery_attempts\n";
    for (unsigned run = 0; run < runs; ++run)
    {
        paired::AnalysisContext ctx(*c.input);
        auto classified = ClassifyAndValidate(ctx,&attempts,id,run);
        double classification = classified.classificationMs;
        for (unsigned warm = 0; warm < 2; ++warm)
        {
            auto r = warm ? mother::ValidatePrepared(ctx, classified.candidate) : classified.report;
            const auto &t = r.timing;
            const auto audit = Audit(c, classified.candidate.patch, r);
#ifdef BSPLINE_MOTHER_STREAMING
            if (run == 0 && warm == 0)
                ExportExtraction(*ctx.motherSource, file.string() + ".extraction.json");
#endif
            Require(r.sourceCacheHits == (warm ? 1 : classified.attempts-1) && r.sourceExtractions == 1 - warm,
                    "Unexpected cache mode");
            csv << id << ',' << TypeName(c.truth.type) << ',' << run << ','
                << (warm ? "warm" : "cold") << ',' << classified.attempts << ','
                << PairedValidationStatusName(r.status) << ',';
            if (!warm)
                csv << classification;
            csv << ',' << t.totalMs << ',';
            if (!warm)
                csv << classification + t.totalMs;
            csv << ',' << t.preparation.elapsedMs << ',' << t.surfaceBounds.elapsedMs << ','
                << t.subdivision.elapsedMs << ',' << t.boundaryChecks.elapsedMs << ','
                << t.reporting.elapsedMs << ',' << r.sourceExtractionMs << ',' << r.sourceBoundaryMs
                << ',' << r.targetPreparationMs << ',' << r.mappingMs << ',' << r.sampleMs << ',' << r.residualBuildMs
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
            csv << ',' << audit.position << ',' << audit.edge << ',' << audit.tangentU << ','
                << audit.tangentV << ',' << audit.joint << ',' << c.input->DegreeU() << ','
                << c.input->DegreeV() << ','
                << c.input->ControlPoints().size() * c.input->ControlPoints().front().size() << ','
                << r.tangentResidualFits << ',' << r.tangentResidualSplits << ',' << (warm ? 1 : classified.attempts);
            if (!warm)
            {
                const auto &d=classified.diagnostics;
                double stages=d.extractionSamplingMs+d.algebraicFitMs+d.geometricFitMs+
                    d.candidateScoringMs+d.patchRecoveryMs;
                Require(stages <= classification+1e-9, "Overlapping classification timers");
                csv << ',' << d.extractionSamplingMs << ',' << d.algebraicFitMs << ','
                    << d.geometricFitMs << ',' << d.candidateScoringMs << ',' << d.patchRecoveryMs
                    << ',' << classification-stages << ',' << d.curveExtractions << ','
                    << d.curvePointEvaluations << ',' << d.sourcePointEvaluations << ','
                    << d.pointCacheHits << ',' << d.fitPasses << ',' << d.geometricIterations << ','
                    << d.fastCircleFits << ',' << d.samplingRefinements << ',' << d.classificationStages
                    << ',' << d.rawCandidateCount << ',' << d.distanceEvaluations << ','
                    << d.sharedMotherScores << ',' << d.duplicateCandidates << ',' << d.patchRecoveryAttempts;
            }
            else for (unsigned i=0;i<20;++i) csv << ',';
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
        Require(argc == 5, "Usage: dataset id csv runs|coverage");
        if (std::string(argv[4]) == "coverage")
            Coverage(argv[1],static_cast<unsigned>(std::stoul(argv[2])),argv[3]);
        else Run(argv[1], static_cast<unsigned>(std::stoul(argv[2])), argv[3],
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
