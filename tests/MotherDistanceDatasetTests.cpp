/**
 * @file MotherDistanceDatasetTests.cpp
 * @brief 当前等参固定数据的母面验证回归驱动。
 * 使用真实候选和冷上下文，有理输入单列精确识别；独立审计与文件输出不进入算法计时。
 */
#include "Dataset500.h"
#include "ExtremeDataset.h"
#include "MotherDistanceInternal.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace
{
using Row = std::map<std::string, std::string>;
const std::vector<std::string> Columns = {
    "id", "type", "output_type", "same_type", "rational", "collapsed", "degree_u", "degree_v", "controls", "run", "backend",
    "status", "accepted", "expected_outcome", "audit", "candidates", "validation_ran", "optimization_ran",
    "classification_compute_ms", "classification_validation_ms", "classification_ms",
    "optimization_compute_ms", "optimization_validation_ms", "algorithm_ms", "api_wall_ms",
    "classification_validation_calls", "preparation_ms", "position_ms", "subdivision_ms", "boundary_ms",
    "reporting_ms", "source_extraction_ms", "source_boundary_ms", "target_preparation_ms", "mapping_ms",
    "sample_ms", "residual_build_ms", "residual_bound_ms", "boundary_position_ms", "boundary_angle_ms",
    "surface_nodes", "boundary_position_nodes", "boundary_nodes", "surface_depth", "mapping_attempts",
    "source_cells", "extracted_cells", "cached_cell_hits", "source_cache_hits", "tightened_cells",
    "boundary_fallbacks", "extraction_plan_ms", "source_contraction_ms", "source_envelope_ms",
    "mother_bound", "edge_bound", "tangent_u_bound_deg", "tangent_v_bound_deg", "joint_tangent_bound_deg",
    "joint_sine_squared_bound", "joint_tangent_failure", "source_certificate_bound",
    "qualified", "qualification_message", "tangent_residual_fits", "tangent_residual_splits",
    "class_sampling_ms", "class_algebraic_ms", "class_geometric_ms", "class_scoring_ms", "class_recovery_ms", "class_other_ms",
    "curve_extractions", "curve_point_evaluations", "score_point_evaluations", "point_cache_hits", "fit_passes",
    "geometric_iterations", "fast_circle_fits", "sampling_refinements", "classification_stages", "raw_candidates",
    "distance_evaluations", "shared_mother_scores", "duplicate_candidates", "patch_recovery_attempts", "audit_position_max",
    "audit_edge_max", "audit_tangent_u_max", "audit_tangent_v_max", "audit_joint_tangent_max", "message"};
void Require(bool ok, const std::string &s) { if (!ok) throw std::runtime_error(s); }
double Ms(std::chrono::steady_clock::time_point start)
{ return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count(); }
template<class T> void Set(Row &r, const std::string &key, T value)
{ std::ostringstream s; s << std::setprecision(17) << value; r[key] = s.str(); }
void Write(std::ostream &s, const Row &r)
{
    for (size_t i = 0; i < Columns.size(); ++i)
    {
        if (i) s << ',';
        auto it = r.find(Columns[i]);
        if (it == r.end()) continue;
        s << '"';
        for (char c : it->second) { if (c == '"') s << '"'; s << (c == '\n' || c == '\r' ? ' ' : c); }
        s << '"';
    }
    s << '\n'; s.flush();
}
double Map(const PairedParameterMap &m, double t)
{
    if (t == m.source.back()) return m.target.back();
    auto it = std::upper_bound(m.source.begin(), m.source.end(), t);
    Require(it != m.source.begin() && it != m.source.end(), "Audit map coverage");
    size_t i = size_t(it - m.source.begin() - 1);
    return m.target[i] + (t - m.source[i]) / (m.source[i + 1] - m.source[i]) * (m.target[i + 1] - m.target[i]);
}
// 在平移后的局部坐标中审计，避免 1e6/1e9 平移下 SDK 导数及世界坐标相减的相消；不改变算法输入、实测候选或证书。
sggk::BSplineSurface CenteredAuditSurface(const sggk::BSplineSurface &s, V anchor)
{
    sggk::BSplineSurface out(s);
    const auto &cp = s.ControlPoints();
    for (unsigned i=0;i<cp.size();++i)
        for (unsigned j=0;j<cp[i].size();++j)
            out.SetControlPoint(Point(Get(cp[i][j])-anchor),i,j);
    return out;
}
void AuditPosition(const sggk::BSplineSurface &s, const AnalyticParameters &p, double bound, Row &row)
{
    V anchor = Get(s.ControlPoints().front().front());
    auto centered = CenteredAuditSurface(s,anchor);
    auto local = p; local.origin = Point(Get(p.origin)-anchor);
    SurfaceSampleContext samples(centered);
    double maximum = 0; bool valid = true;
    for (unsigned i = 0; i <= 12; ++i)
        for (unsigned j = 0; j <= 12; ++j)
        {
            double u = s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 12;
            double v = s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * j / 12;
            double d = Distance(local, Point(samples.At(u, v, sggk::DerivType::NONE, sggk::DerivType::NONE, false).point));
            valid &= std::isfinite(d) && d < 1e-5 && d <= bound + 1e-10;
            maximum = std::max(maximum, d);
        }
    Set(row, "audit_position_max", maximum);
    Require(valid,"Independent centered mother-distance audit");
}
void AuditEdges(const sggk::BSplineSurface &s, const AnalyticPatch &p,
                const MotherDistanceValidationReport &r, Row &row)
{
    auto chart = paired::BuildChart(p, false);
    V anchor = Get(s.ControlPoints().front().front());
    auto centered = CenteredAuditSurface(s,anchor);
    auto local = p.parameters; local.origin = Point(Get(local.origin)-anchor);
    auto frame = FrameOf(local);
    SurfaceSampleContext samples(centered);
    double maxPosition = 0, maxU = 0, maxV = 0, maxJoint = 0;
    bool positionOK = true, tangentOK = true;
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        for (unsigned i = 0; i <= 16; ++i)
        {
            double u = edge < 2 ? (edge ? s.MaxParamU() : s.MinParamU())
                               : s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 16;
            double v = edge >= 2 ? (edge == 3 ? s.MaxParamV() : s.MinParamV())
                                : s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * i / 16;
            double q[2] = {Map(r.uMap, u), Map(r.vMap, v)}, angle[2]{};
            if (p.mapping.reversedU) q[0] = 1 - q[0];
            if (p.mapping.reversedV) q[1] = 1 - q[1];
            if (p.mapping.swappedUV) std::swap(q[0], q[1]);
            for (unsigned d = 0; d < 2; ++d)
            {
                unsigned n = unsigned(chart.angles[d].size() - 1);
                double x = std::clamp(q[d], 0., 1.) * n;
                unsigned k = std::min(n - 1, unsigned(x));
                double t = x - k, a = chart.angles[d][k], b = chart.angles[d][k + 1];
                angle[d] = chart.angular[d] ? a + 2 * std::atan2(t * std::sin((b - a) / 2),
                    1 - t + t * std::cos((b - a) / 2)) : a + t * (b - a);
            }
            const auto &sample = samples.At(u, v);
            V target = frame.global(CanonicalPoint(p.parameters, angle[0], angle[1]));
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
            double distance = Norm(sample.point - target);
            auto angleOf = [](V x, V y) {
                Require(Norm(x) > 0 && Norm(y) > 0, "Degenerate audit tangent");
                return std::atan2(Norm(Cross(x, y)), Dot(x, y));
            };
            double au = angleOf(sample.du, du), av = angleOf(sample.dv, dv);
            double joint = std::asin(std::min(1., std::sqrt(std::pow(std::sin(au), 2) +
                                                          std::pow(std::sin(av), 2)))) * 180 / Pi;
            positionOK &= std::isfinite(distance) && distance < 1e-5 && distance <= r.boundaryPosition.upperBound + 1e-10;
            tangentOK &= std::isfinite(au) && std::isfinite(av) && au < Pi / 2 && av < Pi / 2 &&
                au * 180 / Pi <= r.boundaryTangentU.upperBound + 1e-7 &&
                av * 180 / Pi <= r.boundaryTangentV.upperBound + 1e-7 &&
                joint < .1 && joint <= r.boundaryJointTangent.upperBound + 1e-7;
            maxPosition = std::max(maxPosition, distance);
            maxU = std::max(maxU, au * 180 / Pi); maxV = std::max(maxV, av * 180 / Pi);
            maxJoint = std::max(maxJoint, joint);
        }
    }
    Set(row, "audit_edge_max", maxPosition); Set(row, "audit_tangent_u_max", maxU);
    Set(row, "audit_tangent_v_max", maxV); Set(row, "audit_joint_tangent_max", maxJoint);
    Require(positionOK, "Independent centered paired edge-position audit");
    Require(tangentOK, "Independent centered joint U/V tangent audit");
}
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
    sum.boundaryPositionNodes += r.boundaryPositionNodes;
    sum.surfaceDepth = std::max(sum.surfaceDepth,r.surfaceDepth);
    sum.boundaryNodes += r.boundaryNodes;
    sum.mappingAttempts += r.mappingAttempts;
    sum.extractedCells += r.extractedCells;
    sum.cachedCellHits += r.cachedCellHits;
    sum.tightenedCells += r.tightenedCells;
    sum.boundaryFallbacks += r.boundaryFallbacks;
    sum.tangentResidualFits += r.tangentResidualFits;
    sum.tangentResidualSplits += r.tangentResidualSplits;
}
void CopyCosts(MotherDistanceValidationReport &out,const MotherDistanceValidationReport &sum)
{
    out.timing = sum.timing;
    out.sourceExtractionMs = sum.sourceExtractionMs;
    out.sourceBoundaryMs = sum.sourceBoundaryMs;
    out.targetPreparationMs = sum.targetPreparationMs;
    out.mappingMs = sum.mappingMs;
    out.sampleMs = sum.sampleMs;
    out.residualBuildMs = sum.residualBuildMs;
    out.residualBoundMs = sum.residualBoundMs;
    out.boundaryPositionMs = sum.boundaryPositionMs;
    out.boundaryAngleMs = sum.boundaryAngleMs;
    out.extractionPlanMs = sum.extractionPlanMs;
    out.sourceContractionMs = sum.sourceContractionMs;
    out.sourceEnvelopeMs = sum.sourceEnvelopeMs;
    out.sourceExtractions = sum.sourceExtractions;
    out.sourceCacheHits = sum.sourceCacheHits;
    out.surfaceNodes = sum.surfaceNodes;
    out.boundaryPositionNodes = sum.boundaryPositionNodes;
    out.surfaceDepth = sum.surfaceDepth;
    out.boundaryNodes = sum.boundaryNodes;
    out.mappingAttempts = sum.mappingAttempts;
    out.extractedCells = sum.extractedCells;
    out.cachedCellHits = sum.cachedCellHits;
    out.tightenedCells = sum.tightenedCells;
    out.boundaryFallbacks = sum.boundaryFallbacks;
    out.tangentResidualFits = sum.tangentResidualFits;
    out.tangentResidualSplits = sum.tangentResidualSplits;
}
Row Measure(const DatasetCase &c, double sourceBound, unsigned run, bool strictType, std::ostream *attempts)
{
    Row row;
    Set(row, "id", c.id); row["type"] = TypeName(c.truth.type);
    Set(row, "rational", c.input->IsRational()); Set(row, "collapsed", c.collapsedMask);
    Set(row, "degree_u", c.input->DegreeU()); Set(row, "degree_v", c.input->DegreeV());
    Set(row, "controls", c.input->ControlPoints().size() * c.input->ControlPoints().front().size());
    Set(row, "run", run); Set(row, "source_certificate_bound", sourceBound);
    row["backend"] = c.input->IsRational() ? "exact" : "mother";
    row["optimization_ran"] = "0"; row["validation_ran"] = "0"; row["accepted"] = "0";
    row["audit"] = "NotApplicable"; row["expected_outcome"] = "0";
    double compute = 0, validation = 0;
    ClassificationDiagnostics diagnostics;
    MotherDistanceValidationReport sum;
    unsigned calls = 0;
    bool accepted = false;
    std::optional<AnalyticPatch> patch;
    MotherDistanceValidationReport mother;
    RecognitionResult exact;
    auto start = std::chrono::steady_clock::now();
    if (c.input->IsRational())
    {
        exact = Recognize(*c.input);
        Set(row, "api_wall_ms", Ms(start));
        compute = exact.timing.classificationCompute.elapsedMs;
        validation = exact.timing.classificationValidation.elapsedMs;
        Set(row, "classification_validation_calls", exact.timing.classificationValidation.calls);
        Set(row, "validation_ran", exact.timing.classificationValidation.calls != 0);
        row["status"] = StatusName(exact.status); row["message"] = exact.message;
        accepted = exact.status == Status::Success;
        Set(row, "mother_bound", exact.validation.surfaceDistanceBound);
    }
    else
    {
        paired::AnalysisContext context(*c.input);
        context.classification.emplace();
        std::unique_ptr<IsoClassificationSearch> search;
        auto part = std::chrono::steady_clock::now();
        try
        {
            search = std::make_unique<IsoClassificationSearch>(context.samples, ConversionOptions{}, diagnostics);
            compute += Ms(part);
            for (;;)
            {
                part = std::chrono::steady_clock::now();
                auto candidate = search->Next(*context.classification);
                compute += Ms(part);
                if (!candidate)
                {
                    if (!calls) row["status"] = "NoCandidate";
                    row["message"] += "; incremental search exhausted";
                    break;
                }
                ++calls;
                patch = candidate->patch;
                mother = mother::ValidatePrepared(context, *candidate);
                Accumulate(sum, mother);
                row["validation_ran"] = "1";
                row["status"] = PairedValidationStatusName(mother.status); row["message"] = mother.message;
                accepted = mother.Passed();
                if (attempts)
                {
                    Row attempt;
                    Set(attempt,"id",c.id); Set(attempt,"run",run); Set(attempt,"candidates",calls);
                    attempt["status"] = row["status"]; attempt["message"] = row["message"];
                    Set(attempt,"accepted",accepted); attempt["output_type"] = TypeName(patch->parameters.type);
                    Set(attempt,"classification_stages",diagnostics.classificationStages);
                    Set(attempt,"classification_compute_ms",compute);
                    Set(attempt,"classification_validation_ms",mother.timing.totalMs);
                    Set(attempt,"preparation_ms",mother.timing.preparation.elapsedMs);
                    Set(attempt,"position_ms",mother.timing.surfaceBounds.elapsedMs);
                    Set(attempt,"subdivision_ms",mother.timing.subdivision.elapsedMs);
                    Set(attempt,"boundary_ms",mother.timing.boundaryChecks.elapsedMs);
                    Set(attempt,"reporting_ms",mother.timing.reporting.elapsedMs);
                    Set(attempt,"source_cache_hits",mother.sourceCacheHits);
                    Set(attempt,"mother_bound",mother.motherDistance.upperBound);
                    Set(attempt,"edge_bound",mother.boundaryPosition.upperBound);
                    Set(attempt,"joint_tangent_bound_deg",mother.boundaryJointTangent.upperBound);
                    attempt["joint_tangent_failure"] = mother.jointTangentFailure;
                    Write(*attempts,attempt);
                }
                if (accepted) break;
                part = std::chrono::steady_clock::now();
            }
        }
        catch (const Failure &e)
        {
            compute += Ms(part); row["status"] = "ClassificationFailed";
            row["message"] = std::string(StatusName(e.status)) + ": " + e.what();
        }
        part = std::chrono::steady_clock::now();
        search.reset();
        compute += Ms(part);
        Set(row, "api_wall_ms", Ms(start)); // 包含测试日志耗时，不是算法计时。
        Set(row, "candidates", calls);
        Set(row, "classification_validation_calls", calls);
        validation = sum.timing.totalMs;
        CopyCosts(mother,sum);
        if (row["validation_ran"] == "1")
        {
            const auto &t = mother.timing;
            Set(row, "preparation_ms", t.preparation.elapsedMs); Set(row, "position_ms", t.surfaceBounds.elapsedMs);
            Set(row, "subdivision_ms", t.subdivision.elapsedMs); Set(row, "boundary_ms", t.boundaryChecks.elapsedMs);
            Set(row, "reporting_ms", t.reporting.elapsedMs);
#define FIELD(name, value) Set(row, name, mother.value)
            FIELD("source_extraction_ms", sourceExtractionMs); FIELD("source_boundary_ms", sourceBoundaryMs);
            FIELD("target_preparation_ms", targetPreparationMs); FIELD("mapping_ms", mappingMs);
            FIELD("sample_ms", sampleMs); FIELD("residual_build_ms", residualBuildMs); FIELD("residual_bound_ms", residualBoundMs);
            FIELD("boundary_position_ms", boundaryPositionMs); FIELD("boundary_angle_ms", boundaryAngleMs);
            FIELD("surface_nodes", surfaceNodes); FIELD("boundary_position_nodes", boundaryPositionNodes);
            FIELD("boundary_nodes", boundaryNodes); FIELD("surface_depth", surfaceDepth); FIELD("mapping_attempts", mappingAttempts);
            FIELD("source_cells", sourceCells); FIELD("extracted_cells", extractedCells); FIELD("cached_cell_hits", cachedCellHits);
            FIELD("source_cache_hits", sourceCacheHits); FIELD("tightened_cells", tightenedCells); FIELD("boundary_fallbacks", boundaryFallbacks);
            FIELD("extraction_plan_ms", extractionPlanMs); FIELD("source_contraction_ms", sourceContractionMs);
            FIELD("source_envelope_ms", sourceEnvelopeMs); FIELD("mother_bound", motherDistance.upperBound);
            FIELD("edge_bound", boundaryPosition.upperBound);
            FIELD("tangent_u_bound_deg", boundaryTangentU.upperBound);
            FIELD("tangent_v_bound_deg", boundaryTangentV.upperBound);
            FIELD("joint_tangent_bound_deg", boundaryJointTangent.upperBound);
            FIELD("joint_sine_squared_bound", jointTangentSineSquaredUpperBound);
            FIELD("joint_tangent_failure", jointTangentFailure);
            FIELD("tangent_residual_fits", tangentResidualFits); FIELD("tangent_residual_splits", tangentResidualSplits);
#undef FIELD
        }
    }
    if (!c.input->IsRational())
    {
#define DIAG(name, field) Set(row, name, diagnostics.field)
        DIAG("class_sampling_ms", extractionSamplingMs); DIAG("class_algebraic_ms", algebraicFitMs);
        DIAG("class_geometric_ms", geometricFitMs); DIAG("class_scoring_ms", candidateScoringMs);
        DIAG("class_recovery_ms", patchRecoveryMs);
        Set(row,"class_other_ms",compute-diagnostics.extractionSamplingMs-diagnostics.algebraicFitMs-
            diagnostics.geometricFitMs-diagnostics.candidateScoringMs-diagnostics.patchRecoveryMs);
        DIAG("curve_extractions", curveExtractions); DIAG("curve_point_evaluations", curvePointEvaluations);
        DIAG("score_point_evaluations", sourcePointEvaluations); DIAG("point_cache_hits", pointCacheHits);
        DIAG("fit_passes", fitPasses); DIAG("geometric_iterations", geometricIterations);
        DIAG("fast_circle_fits", fastCircleFits); DIAG("sampling_refinements", samplingRefinements);
        DIAG("classification_stages", classificationStages); DIAG("raw_candidates", rawCandidateCount);
        DIAG("distance_evaluations", distanceEvaluations); DIAG("shared_mother_scores", sharedMotherScores);
        DIAG("duplicate_candidates", duplicateCandidates); DIAG("patch_recovery_attempts", patchRecoveryAttempts);
#undef DIAG
    }
    if (patch || accepted)
    {
        auto outputType = c.input->IsRational() ? exact.parameters.type : patch->parameters.type;
        row["output_type"] = TypeName(outputType); Set(row,"same_type",outputType == c.truth.type);
    }
    Set(row, "classification_compute_ms", compute); Set(row, "classification_validation_ms", validation);
    Set(row, "classification_ms", compute + validation); Set(row, "algorithm_ms", compute + validation);
    Set(row, "accepted", accepted);
    // 所有测试真值读取与独立审计均位于算法计时外。
    if (accepted)
    {
        try
        {
            Require(!c.collapsedMask, "Collapsed fixture accepted");
            if (c.input->IsRational())
            {
                Require(exact.surface && exact.domain && (!strictType || exact.parameters.type == c.truth.type) &&
                    exact.validation.level == VerificationLevel::BernsteinBound &&
                    exact.validation.surfaceDistanceBound < 1e-5, "Exact output contract");
                AuditPosition(*c.input, exact.parameters, exact.validation.surfaceDistanceBound, row);
            }
            else
            {
                Require(patch && (!strictType || patch->parameters.type == c.truth.type) && mother.motherDistance.upperBound < 1e-5 &&
                    mother.boundaryPosition.upperBound < 1e-5 && mother.boundaryJointTangent.upperBound < .1,
                    "Mother output contract");
                Require(mother.sourceExtractions == 1 &&
                    mother.extractedCells == mother.sourceCells, "Cold source/full-domain coverage");
                AuditPosition(*c.input, patch->parameters, mother.motherDistance.upperBound, row);
                AuditEdges(*c.input, *patch, mother, row);
            }
            row["audit"] = "Passed"; row["expected_outcome"] = "1";
        }
        catch (const std::exception &e) { row["audit"] = "Failed"; row["message"] += std::string("; audit: ") + e.what(); }
    }
    else if (c.collapsedMask)
    {
        const auto &status = row["status"];
        bool expected = c.input->IsRational() ? status == "Degenerate" : status == "NoCandidate" ||
            (status == "ClassificationFailed" && row["message"].find("Degenerate:") == 0);
        Set(row, "expected_outcome", expected);
    }
    return row;
}
int Run(const std::filesystem::path &dir, unsigned id, const std::filesystem::path &file, unsigned runs)
{
    const bool standard = dir.filename() == "dataset500";
    const bool business = dir.filename() == "business-extreme200";
    Require(standard || business || dir.filename() == "extreme200", "Not a specialized dataset");
    Require(id < (standard ? 500u : 200u) && runs > 0 && runs <= 100, "Invalid case or run count");
    auto c = ReadDatasetCase(dir / (standard ? CaseFilename(id) : ExtremeFilename(id)));
    std::cout << "BEGIN_SOURCE " << id << std::endl;
    double bound = INFINITY;
    bool qualified = false;
    std::string qualification;
    if (standard)
    {
        CheckDataset500Contract(c,id);
        auto certificate = BoundPair(*c.source,*c.input,c.refitted ? 1e-6 : 1e-9);
        bound = certificate.upper; qualified = certificate.passed;
        qualification = qualified ? "Dataset500 contract and continuous source bound" : "source certificate unresolved";
    }
    else
    {
        auto specs = business ? BusinessExtremeSpecs() : ExtremeSpecs();
        Require(c.id == id && c.pair == id/2 && c.type == specs[id/2].type && c.refitted == bool(id%2),
                "Extreme fixture identity mismatch");
        auto q = QualifyExtreme(c,business);
        bound = q.pair.upper; qualified = q.passed; qualification = q.message;
    }
    std::cout << "END_SOURCE qualified=" << qualified << std::endl;
    std::ofstream out(file), attempts(file.string()+".attempts.csv");
    Require(bool(out) && bool(attempts), "Cannot open CSV");
    for (auto stream : {&out,&attempts})
    {
        for (size_t i = 0; i < Columns.size(); ++i) *stream << (i ? "," : "") << Columns[i];
        *stream << '\n'; stream->flush();
    }
    // 不做逐案例预热；报告三次全新上下文的运行，包括第一次。
    unsigned unexpected = 0;
    for (unsigned run = 0; run < runs; ++run)
    {
        std::cout << "BEGIN_ALGORITHM run=" << run << std::endl;
        auto row = Measure(c,bound,run,standard,&attempts);
        Set(row,"qualified",qualified); row["qualification_message"] = qualification;
        if (!qualified) row["expected_outcome"] = "0";
        unexpected += row["expected_outcome"] != "1";
        Write(out,row);
        std::cout << "END_ALGORITHM run=" << run << " status=" << row["status"]
                  << " ms=" << row["algorithm_ms"] << std::endl;
    }
    std::cout << "case=" << id << " runs=" << runs << " unexpected=" << unexpected << '\n';
    return unexpected ? 2 : 0;
}
}
int main(int argc, char **argv)
{
    sggk::init(); int result = 1;
    try
    {
        Require(argc == 6 && std::string(argv[1]) == "--case", "Usage: --case dataset id csv runs");
        result = Run(argv[2], unsigned(std::stoul(argv[3])), argv[4], unsigned(std::stoul(argv[5])));
    }
    catch (const std::exception &e) { std::cerr << "FAIL " << e.what() << '\n'; }
    sggk::fini(); return result;
}
