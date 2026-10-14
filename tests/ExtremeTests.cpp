#include "ApproximateInternal.h"
#include "ExtremeDataset.h"
#include "SurfaceConversion.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace
{
std::string Quote(std::string text)
{
    std::string out = "\"";
    for (char c : text)
    {
        if (c == '\"')
            out += '\"';
        out += (c == '\n' || c == '\r') ? ' ' : c;
    }
    return out + '\"';
}
double DirectionAngle(V a, V b)
{
    a = Unit(a);
    b = Unit(b);
    return std::atan2(Norm(Cross(a, b)), std::abs(Dot(a, b))) * 180 / Pi;
}
struct Audit
{
    double distance = 0, tangent = 0, normal = 0;
    bool valid = true;
};
std::vector<double> Parameters(const sggk::RealArray &knots, unsigned subdivisions)
{
    std::vector<double> out;
    for (size_t i = 1; i < knots.size(); ++i)
        for (unsigned j = 0; j < subdivisions; ++j)
            out.push_back(knots[i - 1] + (knots[i] - knots[i - 1]) * j / subdivisions);
    out.push_back(knots.back());
    return out;
}
std::array<V, 3> EvaluateSource(const sggk::BSplineSurface &s, unsigned e, double t)
{
    double u = e < 2 ? (e == 0 ? 0 : 1) : t, v = e < 2 ? t : (e == 2 ? 0 : 1);
    sggk::Point3D p;
    auto d = s.CalcPartialD1(u, v, u == 1 ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT,
                             v == 1 ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT, p);
    return {Get(p), Get(d[e < 2 ? 1 : 0]), Cross(Get(d[0]), Get(d[1]))};
}
Audit AuditSource(const DatasetCase &c)
{
    Audit a;
    auto us = Parameters(c.input->KnotsU(), 4), vs = Parameters(c.input->KnotsV(), 4);
    for (double u : us)
        for (double v : vs)
            a.distance =
                std::max(a.distance, ExtremeMotherDistance(c.truth, Get(c.input->CalcPoint(u, v))));
    for (unsigned e = 0; e < 4; ++e)
        for (double t : Parameters(e < 2 ? c.input->KnotsV() : c.input->KnotsU(), 32))
        {
            auto input = EvaluateSource(*c.input, e, t), source = EvaluateSource(*c.source, e, t);
            a.tangent = std::max(a.tangent, DirectionAngle(input[1], source[1]));
            a.normal = std::max(a.normal, DirectionAngle(input[2], source[2]));
        }
    return a;
}
V PCurveTangent(const GeneralAnalyticPatch &patch, unsigned e, double t)
{
    const auto &edge = patch.PCurves()[e];
    auto uv = edge.CalcPoint(t);
    double u = uv.X(), v = uv.Y(), c = std::cos(u), s = std::sin(u), cv = std::cos(v),
           sv = std::sin(v);
    const auto &p = patch.Parameters();
    V du, dv;
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        du = {1, 0, 0};
        dv = {0, 1, 0};
        break;
    case sggk::SurfaceType::Cylinder:
        du = {-p.radius * s, p.radius * c, 0};
        dv = {0, 0, 1};
        break;
    case sggk::SurfaceType::Cone:
        du = {-v * std::tan(p.halfAngle) * s, v * std::tan(p.halfAngle) * c, 0};
        dv = {std::tan(p.halfAngle) * c, std::tan(p.halfAngle) * s, 1};
        break;
    case sggk::SurfaceType::Sphere:
        du = {-p.radius * sv * s, p.radius * sv * c, 0};
        dv = {p.radius * cv * c, p.radius * cv * s, p.radius * sv};
        break;
    case sggk::SurfaceType::Torus:
        du = {-(p.majorRadius + p.radius * cv) * s, (p.majorRadius + p.radius * cv) * c, 0};
        dv = {-p.radius * sv * c, -p.radius * sv * s, p.radius * cv};
        break;
    default:
        throw std::runtime_error("Unknown output type");
    }
    V z = Unit(Get(p.axis)), x = Unit(Get(p.xDirection) - z * Dot(z, Get(p.xDirection))),
      y = Cross(z, x);
    V local = du * (edge.end.X() - edge.start.X()) + dv * (edge.end.Y() - edge.start.Y());
    return x * local.x + y * local.y + z * local.z;
}
Audit AuditOutput(const sggk::BSplineSurface &s, const GeneralAnalyticPatch &p)
{
    Audit a;
    for (double u : Parameters(s.KnotsU(), 4))
        for (double v : Parameters(s.KnotsV(), 4))
            a.distance =
                std::max(a.distance, ExtremeMotherDistance(p.Parameters(), Get(s.CalcPoint(u, v))));
    for (unsigned e = 0; e < 4; ++e)
        for (double t : Parameters(e < 2 ? s.KnotsV() : s.KnotsU(), 32))
        {
            auto source = EvaluateSource(s, e, t);
            auto sq = [&](double f) {
                V d = Get(p.CalcBoundaryPoint(BoundaryId(e), f)) - source[0];
                return Dot(d, d);
            };
            // Bracket nearest fraction independently of the source parameter. The
            // rational angle map is not affine in the source spline parameter.
            unsigned best = 0;
            double value = sq(0);
            for (unsigned j = 1; j <= 32; ++j)
            {
                double z = sq(double(j) / 32);
                if (z < value)
                {
                    value = z;
                    best = j;
                }
            }
            double lo = double(best ? best - 1 : 0) / 32, hi = double(std::min(32u, best + 1)) / 32;
            constexpr double g = .6180339887498948482;
            double x = hi - g * (hi - lo), y = lo + g * (hi - lo), fx = sq(x), fy = sq(y);
            for (unsigned j = 0; j < 60; ++j)
            {
                if (fx < fy)
                {
                    hi = y;
                    y = x;
                    fy = fx;
                    x = hi - g * (hi - lo);
                    fx = sq(x);
                }
                else
                {
                    lo = x;
                    x = y;
                    fx = fy;
                    y = lo + g * (hi - lo);
                    fy = sq(y);
                }
            }
            double f = (lo + hi) / 2;
            for (double endpoint : {0., 1.})
                if (sq(endpoint) < sq(f))
                    f = endpoint;
            V projected = Get(p.CalcBoundaryPoint(BoundaryId(e), f));
            a.distance = std::max(a.distance, Norm(projected - source[0]));
            a.tangent = std::max(a.tangent, DirectionAngle(source[1], PCurveTangent(p, e, f)));
            a.normal =
                std::max(a.normal,
                         DirectionAngle(source[2], ExtremeMotherNormal(p.Parameters(), projected)));
        }
    return a;
}
std::string FailureStage(const SurfaceConversionResult &r)
{
    if (r.Succeeded())
        return "success";
    auto status = std::visit(
        [](const auto &x) -> std::string {
            using T = std::decay_t<decltype(x)>;
            if constexpr (std::is_same_v<T, RecognitionResult>)
                return StatusName(x.status);
            else
                return ConversionStatusName(x.status);
        },
        r.details);
    return status;
}
int Run(const std::filesystem::path &dir, const std::filesystem::path &csvPath, unsigned first,
        unsigned last)
{
    if (first >= last || last > 200)
        throw std::runtime_error("Expected 0 <= first < last <= 200");
    std::ofstream csv(csvPath);
    if (!csv)
        throw std::runtime_error("Cannot write result CSV");
    auto fitPath = csvPath.parent_path() / (csvPath.stem().string() + "-fits.csv");
    std::ofstream fits(fitPath);
    if (!fits) throw std::runtime_error("Cannot write curve fit diagnostics");
    fits << std::setprecision(17)
         << "id,path,direction,parameter,boundary,edge,samples,candidate_tolerance,numeric_floor,"
            "transverse_ratio,line_error,circle_error,radius,sweep_angle,sagitta,has_line,has_circle,circle_state\n";
    csv << std::setprecision(17)
        << "id,group,type,source_kind,qualified,controls,rational,pair_bound,source_distance_"
           "sampled,source_tangent_sampled_deg,source_normal_sampled_deg,status,path,output_type,"
           "classify_ms,optimize_ms,classify_optimize_ms,validate_ms,total_ms,geometry_bound,"
           "boundary_bound,output_distance_sampled,output_tangent_sampled_deg,output_normal_"
           "sampled_deg,angle_evidence,iterations,candidates,pieces,attempts,message\n";
    unsigned qualified = 0, success = 0, failed = 0;
    auto specs = ExtremeSpecs();
    for (unsigned id = first; id < last; ++id)
    {
        const auto &spec = specs[id / 2];
        ExtremeQualification q;
        Audit source, output;
        output.valid = false;
        std::string status = "generation_failed", message, attempts, path, outputType,
                    angleEvidence = "not_checked";
        ConversionTiming timing;
        double gb = INFINITY, bb = INFINITY;
        unsigned iterations = 0, candidates = 0, pieces = 0;
        bool rational = false;
        try
        {
            auto c = ReadDatasetCase(dir / ExtremeFilename(id));
            if (c.id != id || c.pair != id / 2 || c.type != spec.type || c.refitted != bool(id % 2))
                throw std::runtime_error("Fixture identity mismatch");
            q = QualifyExtreme(c);
            rational = c.input->IsRational();
            message = q.message;
            if (q.passed)
            {
                ++qualified;
                source = AuditSource(c);
                // Only the surface is supplied. Truth, radii, groups and source
                // representation labels never reach Convert's candidate logic.
                std::cout << "BEGIN_CONVERT " << id << std::endl;
                auto r = Convert(*c.input);
                std::cout << "END_CONVERT " << id << std::endl;
                timing = r.Timing();
                path = ConversionPathName(r.path);
                outputType = TypeName(r.Type());
                message = r.Message();
                status = FailureStage(r);
                for (const auto &attempt : r.attempts)
                {
                    attempts +=
                        std::string(ConversionPathName(attempt.path)) + ":" +
                        (attempt.skipped ? "skipped" : (attempt.succeeded ? "success" : "failed")) +
                        ":" + attempt.message + " | ";
                    for (const auto &d : attempt.isoCurveFits)
                        fits << id << ',' << ConversionPathName(attempt.path) << ',' << d.fixedDirection
                             << ',' << d.fixedParameter << ',' << d.boundary << ',' << unsigned(d.edge)
                             << ',' << d.sampleCount << ',' << d.candidateTolerance << ',' << d.numericFloor
                             << ',' << d.transverseRatio << ',' << d.lineError << ',' << d.circleError
                             << ',' << d.radius << ',' << d.sweepAngle << ',' << d.sagitta << ',' << d.hasLine
                             << ',' << d.hasCircle << ',' << CircleFitStateName(d.circleState) << '\n';
                }
                fits.flush();
                gb = r.pcurveValidation.distanceUpperBound;
                bb = r.pcurveValidation.boundaryHausdorffUpperBound;
                pieces = r.pcurveValidation.piecesVisited;
                if (auto p = std::get_if<ConversionResult>(&r.details))
                {
                    iterations = p->statistics.optimizationIterations;
                    candidates = p->statistics.candidateCount;
                    if (p->validation.Passed())
                        angleEvidence = "specialized_continuous_report";
                }
                if (auto p = std::get_if<GeneralConversionResult>(&r.details))
                {
                    iterations = p->statistics.optimizationIterations;
                    candidates = static_cast<unsigned>(p->candidates.size());
                }
                if (r.Succeeded())
                {
                    output = AuditOutput(*c.input, *r.patch);
                    output.valid = true;
                    if (angleEvidence == "not_checked")
                        angleEvidence = "sampled_only";
                    bool auditOK =
                        output.distance < 1e-5 && output.tangent < .1 && output.normal < .1;
                    status = !auditOK                   ? "success_with_sampled_violation"
                             : r.Type() == c.truth.type ? "same_type_success"
                                                        : "alternate_type_success";
                    if (auditOK)
                        ++success;
                }
            }
            else
                status = "qualification_unresolved";
        }
        catch (const std::exception &e)
        {
            message = e.what();
            status = q.passed ? "test_exception" : "generation_failed";
        }
        if (status != "same_type_success" && status != "alternate_type_success")
            ++failed;
        csv << id << ',' << spec.group << ',' << TypeName(spec.truth.type) << ','
            << (id % 2 ? "polynomial" : "analytic_nurbs") << ',' << q.passed << ',' << q.controls
            << ',' << rational << ',' << q.pair.upper << ',' << source.distance << ','
            << source.tangent << ',' << source.normal << ',' << status << ',' << path << ','
            << outputType << ',' << timing.classification.elapsedMs << ','
            << timing.optimization.elapsedMs << ','
            << timing.classification.elapsedMs + timing.optimization.elapsedMs << ','
            << timing.validation.elapsedMs << ',' << timing.totalMs << ',' << gb << ',' << bb << ','
            << (output.valid ? output.distance : INFINITY) << ','
            << (output.valid ? output.tangent : INFINITY) << ','
            << (output.valid ? output.normal : INFINITY) << ',' << angleEvidence << ','
            << iterations << ',' << candidates << ',' << pieces << ',' << Quote(attempts) << ','
            << Quote(message) << '\n';
        csv.flush();
        std::cout << id << " " << spec.group << " " << status
                  << " cls+opt=" << timing.classification.elapsedMs + timing.optimization.elapsedMs
                  << " ms" << std::endl;
    }
    std::cout << "qualified=" << qualified << " accepted_with_sampled_audit=" << success
              << " failed_or_unresolved=" << failed << std::endl;
    if (!csv)
        throw std::runtime_error("Result CSV write failed");
    return failed ? 1 : 0;
}
int Modules(const std::filesystem::path &dir, const std::filesystem::path &csvPath)
{
    std::ofstream csv(csvPath);
    if (!csv)
        throw std::runtime_error("Cannot write module CSV");
    csv << std::setprecision(17)
        << "test,type,fixture,initial_worst,final_worst,iterations,analytic_blocks,fallback_blocks,"
           "jacobian_best_relative_error,geometry,tangent,normal,passed,message\n";
    unsigned failed = 0;
    // New-data-only biased-seed probes: ordinary translated-base fixtures,
    // angle-gate cone, and the high-aspect plane. No old regression input.
    for (unsigned id : {197u, 159u, 105u, 165u, 171u})
    {
        bool passed = false;
        std::string message;
        double before = INFINITY, after = INFINITY, jacError = INFINITY;
        unsigned iterations = 0, ab = 0, fb = 0;
        sggk::SurfaceType type = sggk::SurfaceType::Unknown;
        try
        {
            auto c = ReadDatasetCase(dir / ExtremeFilename(id));
            type = c.truth.type;
            auto p = RecoverApproxPatch(*c.input, c.truth);
            if (type == sggk::SurfaceType::Plane)
                p.parameters.origin =
                    Point(Get(p.parameters.origin) + Get(p.parameters.axis) * 2e-5);
            else if (type == sggk::SurfaceType::Cone)
                p.parameters.halfAngle += 2e-5 / std::max(1., p.domain[2]);
            else
                p.parameters.radius += 2e-5;
            MaterializePatch(p);
            auto samples = SampleSurface(*c.input, {});
            ConversionOptions o;
            EvaluationWorkspace w;
            ConversionStatistics stats;
            EvaluateResiduals(p, samples, o.validation, w, w.current);
            before = w.current.worst;
            auto active = ParametersFor(type, ParameterLevel::SizeOrOffset);
            double scale = SampleScale(p, samples);
            RefinementDiagnostics d;
            bool built = BuildJacobian(p, samples, o, active, scale, w, d);
            ab = d.analyticSampleBlocks;
            fb = d.fallbackSampleBlocks;
            if (built)
            {
                for (double h : {1e-6, 1e-7, 1e-8})
                {
                    double delta = 0, norm = 1;
                    for (size_t k = 0; k < active.indices.size(); ++k)
                    {
                        std::vector<double> step(
                            ParametersFor(type, ParameterLevel::FullPatch).indices.size());
                        step[active.indices[k]] = h;
                        EvaluationWorkspace scratch;
                        EvaluateResiduals(Increment(p, step, scale), samples, o.validation, scratch,
                                          scratch.plus);
                        step[active.indices[k]] = -h;
                        EvaluateResiduals(Increment(p, step, scale), samples, o.validation, scratch,
                                          scratch.minus);
                        for (size_t j = 0; j < w.jacobian.size(); ++j)
                        {
                            double fd =
                                (scratch.plus.leastSquares[j] - scratch.minus.leastSquares[j]) /
                                (2 * h);
                            delta = std::max(delta, std::abs(fd - w.jacobian[j][k]));
                            norm = std::max(norm, std::abs(fd));
                        }
                    }
                    jacError = std::min(jacError, delta / norm);
                }
            }
            RefineCandidateLevel(p, samples, o, stats, ParameterLevel::SizeOrOffset, 30, false, w);
            EvaluateResiduals(p, samples, o.validation, w, w.current);
            after = w.current.worst;
            iterations = stats.optimizationIterations;
            passed = built && jacError < 1e-3 && after < 1 && after <= before;
            message = stats.refinementRuns.empty()
                          ? "no_refinement_record"
                          : RefinementStopReasonName(stats.refinementRuns.back().stopReason);
        }
        catch (const std::exception &e)
        {
            message = e.what();
        }
        if (!passed)
            ++failed;
        csv << "biased_seed," << TypeName(type) << ',' << id << ',' << before << ',' << after << ','
            << iterations << ',' << ab << ',' << fb << ',' << jacError
            << ",not_checked,not_checked,not_checked," << passed << ',' << Quote(message) << '\n';
    }
    for (unsigned mode = 0; mode < 2; ++mode)
    {
        AnalyticPatch p;
        p.parameters.type = sggk::SurfaceType::Plane;
        p.domain = {-.00255, .00255, -.00255, .00255};
        MakeNative(p);
        auto source = NormalizeDatasetSurface(*p.surface->ToBSpline(*p.sdkDomain));
        PadDatasetSurface(*source, 8, 8, 1);
        p = RecoverApproxPatch(*source, p.parameters);
        double a = .15 * Pi / 180;
        if (mode == 0)
            p.parameters.xDirection = Direction({std::cos(a), std::sin(a), 0});
        else
            p.parameters.axis = Direction({0, -std::sin(a), std::cos(a)});
        MaterializePatch(p);
        auto check = ValidateConversion(*source, p);
        bool passed = check.geometry.state == CheckState::Passed &&
                      (mode == 0 ? check.boundaryTangent.state : check.boundaryNormal.state) ==
                          CheckState::Exceeded &&
                      !check.Passed();
        if (!passed)
            ++failed;
        csv << (mode == 0 ? "tangent_rejection" : "normal_rejection")
            << ",Plane,tiny_square,0,0,0,0,0,0," << CheckStateName(check.geometry.state) << ','
            << CheckStateName(check.boundaryTangent.state) << ','
            << CheckStateName(check.boundaryNormal.state) << ',' << passed << ','
            << Quote(check.message) << '\n';
    }
    if (!csv)
        throw std::runtime_error("Module CSV write failed");
    std::cout << "module failures=" << failed << std::endl;
    return failed ? 1 : 0;
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init();
    int status = 2;
    try
    {
        if (argc == 3 && std::string(argv[1]) == "--generate")
        {
            GenerateExtremeDataset(argv[2]);
            status = 0;
        }
        else if ((argc == 4 || argc == 6) && std::string(argv[1]) == "--run")
            status = Run(argv[2], argv[3], argc == 6 ? std::stoul(argv[4]) : 0,
                         argc == 6 ? std::stoul(argv[5]) : 200);
        else if (argc == 4 && std::string(argv[1]) == "--modules")
            status = Modules(argv[2], argv[3]);
        else
            throw std::runtime_error("Usage: --generate dir | --run dir result.csv [first last] | "
                                     "--modules dir result.csv");
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
    }
    sggk::fini();
    return status;
}
