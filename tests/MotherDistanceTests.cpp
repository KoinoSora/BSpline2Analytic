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
unsigned checks = 0;
void Require(bool ok, const std::string &message)
{
    ++checks;
    if (!ok) throw std::runtime_error(message);
}
double Ms(std::chrono::steady_clock::time_point t)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
}
std::string Quote(std::string s)
{
    std::string out = "\"";
    for (char c : s) { if (c == '"') out += '"'; out += c; }
    return out + '"';
}
void Timing(const ValidationTiming &t)
{
    double sum = t.preparation.elapsedMs + t.surfaceBounds.elapsedMs + t.subdivision.elapsedMs +
        t.boundaryChecks.elapsedMs + t.reporting.elapsedMs;
    Require(t.totalMs >= 0 && std::abs(sum - t.totalMs) < 1e-7 * std::max(1., t.totalMs),
        "Exclusive five-stage timer mismatch");
}
std::vector<RecoveredCandidate> Classify(paired::AnalysisContext &context)
{
    std::vector<RecoveredCandidate> recovered;
    context.classification.emplace();
    auto list = ClassifyApproximate(context.samples.source, {}, nullptr, true, &context.samples,
        &recovered, &*context.classification);
    Require(list.size() == recovered.size(), "Candidate/trace count differs");
    return recovered;
}
void Spot(const sggk::BSplineSurface &s, const AnalyticPatch &p,
          const MotherDistanceValidationReport &r)
{
    SurfaceSampleContext samples(s);
    for (unsigned i = 0; i <= 12; ++i)
        for (unsigned j = 0; j <= 12; ++j)
        {
            double u = s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 12;
            double v = s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * j / 12;
            double d = Distance(p.parameters, Point(samples.At(u, v, sggk::DerivType::NONE,
                sggk::DerivType::NONE, false).point));
            Require(d < 1e-5 && d <= r.motherDistance.upperBound + 1e-11,
                "Independent point audit exceeds mother bound");
        }
}
sggk::BSplineSurfacePtr Plane(unsigned mode)
{
    sggk::Point3DMatrix cp(4, std::vector<sggk::Point3D>(4));
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
        {
            double x = mode == 4 ? double(i) / 6 + double(i * (i ? i - 1 : 0)) / 12 : double(i) / 3;
            cp[i][j] = sggk::Point3D(x * (mode == 4 ? 1 : .001), double(j) / 3, 0);
        }
    if (mode == 1)
        for (unsigned j = 1; j < 3; ++j) cp[1][j] = sggk::Point3D(.001 / 3, double(j) / 3, 8e-6);
    if (mode == 2)
        for (unsigned i = 1; i < 3; ++i) cp[i][0] = sggk::Point3D(double(i) * .001 / 3, 8e-6, 0);
    if (mode == 3)
        for (unsigned i = 1; i < 3; ++i)
            for (unsigned j = 1; j < 3; ++j) cp[i][j] = sggk::Point3D(double(i) * .001 / 3, double(j) / 3, 1e-3);
    return std::make_shared<sggk::BSplineSurface>(3, 3, cp, sggk::RealArray{0, 1},
        sggk::RealArray{0, 1}, sggk::UIntArray{4, 4}, sggk::UIntArray{4, 4});
}
void Unit(const std::filesystem::path &dir)
{
    MotherDistanceValidationOptions affine;
    affine.boundary.enablePiecewiseMapping = false;
    AnalyticPatch plane;
    plane.parameters.type = sggk::SurfaceType::Plane;
    plane.domain = {0, .001, 0, 1};
    auto exact = ValidateMotherDistanceConversion(*Plane(0), plane);
    Require(exact.Passed(), "Exact plane: " + exact.message); Timing(exact.timing);
    auto wrong = plane; wrong.domain[0] += .1; wrong.domain[1] += .1;
    auto rectangle = ValidateMotherDistanceConversion(*Plane(0), wrong, affine);
    Require(!rectangle.Passed() && rectangle.motherDistance.state == CheckState::Passed &&
        rectangle.boundaryPosition.state == CheckState::Exceeded, "Wrong rectangle on same mother accepted");
    for (unsigned mode = 1; mode <= 3; ++mode)
    {
        auto r = ValidateMotherDistanceConversion(*Plane(mode), plane, affine);
        Require(!r.Passed(), "Negative plane fixture accepted"); Timing(r.timing);
        if (mode == 1) Require(r.boundaryNormal.state == CheckState::Exceeded, "Normal defect missed");
        if (mode == 2) Require(r.boundaryTangent.state == CheckState::Exceeded, "Tangent defect missed");
        if (mode == 3) Require(r.motherDistance.state == CheckState::Exceeded, "Interior defect missed");
    }
    wrong = plane; wrong.normalSign = -1;
    Require(!ValidateMotherDistanceConversion(*Plane(0), wrong, affine).Passed(), "Reversed normal accepted");
    auto warped = Plane(4); plane.domain = {0, 1, 0, 1};
    auto mapped = ValidateMotherDistanceConversion(*warped, plane);
    Require(mapped.Passed() && mapped.mappingAttempts == 2 && mapped.uMap.source.size() > 2,
        "Nonlinear monotone correspondence failed: " + mapped.message);
    for (unsigned flags = 1; flags < 8; ++flags)
    {
        auto oriented = Plane(4);
        if (flags & 1) oriented->UReverse();
        if (flags & 2) oriented->VReverse();
        if (flags & 4) oriented->ExchangeUV();
        oriented = NormalizeDatasetSurface(*oriented);
        auto patch = RecoverApproxPatch(*oriented, plane.parameters);
        auto result = ValidateMotherDistanceConversion(*oriented, patch);
        Require(result.Passed(), "Swapped/reversed correspondence failed: " + result.message);
        Timing(result.timing);
    }
    {
        auto translated = Plane(0);
        V delta{1e6, -2e6, 3e6};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                translated->SetControlPoint(Point(Get(translated->ControlPoints()[i][j]) + delta), i, j);
        auto patch = plane; patch.parameters.origin = Point(delta); patch.domain[1] = .001;
        auto result = ValidateMotherDistanceConversion(*translated, patch);
        Require(result.Passed() && result.sourceEnclosureRadius > 1e-12,
            "Large translation lost source enclosure or certificate");
    }
    {
        // Coefficients exceed tolerance, but the entire cubic is below it:
        // z(t)=4.5e-5*t*(1-t)*(1-2*t). A root-only bound must fail, whereas
        // residual subdivision can prove the bound without rebuilding products.
        auto s = Plane(0);
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                s->SetControlPoint(sggk::Point3D(double(i) * .001 / 3, double(j) / 3,
                    i == 1 ? 1.5e-5 : i == 2 ? -1.5e-5 : 0), i, j);
        PairedValidationReport prep;
        auto source = paired::PrepareSource(*s, prep);
        for (unsigned depth : {0u, 12u})
        {
            MotherDistanceValidationReport r;
            RecognitionOptions options; options.maxSubdivisionDepth = depth;
            ValidationTimer timer(r.timing);
            bool passed = mother::ValidatePosition(source, plane.parameters, options, r, timer);
            timer.Finish(); Timing(r.timing);
            Require(passed == (depth != 0), "Residual subdivision/depth-limit contract failed");
            if (passed)
                Require(r.surfaceNodes == 3 && r.surfaceDepth == 1 &&
                    r.timing.subdivision.elapsedMs > 0 && r.motherDistance.upperBound < 1e-5,
                    "Residual split accounting failed");
            else Require(!std::isfinite(r.motherDistance.upperBound), "Partial bound presented as complete");
        }
    }
    auto invalid = affine; invalid.position.requireBound = false;
    Require(ValidateMotherDistanceConversion(*warped, plane, invalid).status ==
        PairedValidationStatus::InvalidInput, "Sample-only mode accepted");
    invalid = affine; invalid.position.distanceTolerance = NAN;
    Require(ValidateMotherDistanceConversion(*warped, plane, invalid).status ==
        PairedValidationStatus::InvalidInput, "Invalid tolerance accepted");
    for (unsigned id : {1u, 101u, 201u, 301u, 391u, 401u, 489u, 497u})
    {
        auto c = ReadDatasetCase(dir / CaseFilename(id)); CheckDataset500Contract(c, id);
        paired::AnalysisContext context(*c.input);
        auto recovered = Classify(context);
        Require(!recovered.empty(), "Representative classification failed");
        auto r = mother::ValidatePrepared(context, recovered.front());
        Require(r.Passed(), "Representative " + std::to_string(id) + ": " + r.message);
        Require(r.sourceEnclosureRadius > 0 && r.sourceExtractions == 1, "Source enclosure missing");
        Timing(r.timing); Spot(*c.input, recovered.front().patch, r);
        auto repeat = mother::ValidatePrepared(context, recovered.front());
        Require(repeat.Passed() && repeat.sourceCacheHits == 1 && repeat.sourceExtractions == 0,
            "Immutable source cache reuse failed");
        if (id >= 301)
        {
            auto displaced = recovered.front(); displaced.patch.parameters.radius += .01;
            auto bad = mother::ValidatePrepared(context, displaced);
            Require(!bad.Passed(), "Changing candidate reused stale residual");
            MotherDistanceValidationOptions budget;
            budget.position.maxValidationPieces = 1;
            auto limited = mother::ValidatePrepared(context, recovered.front(), budget);
            Require(limited.status == PairedValidationStatus::BudgetExceeded &&
                !std::isfinite(limited.motherDistance.upperBound), "Incomplete budget result accepted");
            Require(ValidateMotherDistanceConversion(*c.source, recovered.front().patch).status ==
                PairedValidationStatus::Unsupported, "Rational source accepted");
            auto singular = recovered.front().patch;
            if (id < 400) singular.domain[2] = 0;
            else singular.parameters.majorRadius = singular.parameters.radius;
            Require(ValidateMotherDistanceConversion(*c.input, singular).status ==
                PairedValidationStatus::Unsupported, "Singular analytic chart accepted");
        }
        std::cout << "representative " << id << " preparation_ms=" << r.timing.preparation.elapsedMs
            << " position_ms=" << r.timing.surfaceBounds.elapsedMs
            << " subdivision_ms=" << r.timing.subdivision.elapsedMs
            << " boundary_ms=" << r.timing.boundaryChecks.elapsedMs
            << " reporting_ms=" << r.timing.reporting.elapsedMs
            << " total_ms=" << r.timing.totalMs << std::endl;
    }
    std::cout << "Unit checks=" << checks << " PASSED" << std::endl;
}
void Batch(const std::filesystem::path &dir, const std::filesystem::path &file)
{
    std::ofstream csv(file); Require(bool(csv), "Cannot open CSV");
    csv << std::setprecision(17);
    csv << "id,type,collapsed,run,backend,candidates,status,classification_ms,validation_ms,algorithm_ms,"
        "preparation_ms,position_ms,subdivision_ms,boundary_ms,reporting_ms,source_extraction_ms,"
        "source_boundary_ms,target_preparation_ms,mapping_ms,sample_ms,residual_build_ms,residual_bound_ms,"
        "boundary_position_ms,boundary_angle_ms,surface_nodes,boundary_position_nodes,boundary_nodes,"
        "surface_depth,mapping_attempts,source_cache_hits,mother_or_paired_bound,boundary_position_bound,"
        "tangent_bound,normal_bound,source_certificate_bound,message\n";
    // Both current validators get independent cold source caches and the same
    // real first candidate. Alternate their order to reduce systematic bias.
    for (unsigned warm : {301u, 401u})
    {
        auto c = ReadDatasetCase(dir / CaseFilename(warm));
        for (bool third : {false, true})
        {
            paired::AnalysisContext a(*c.input); auto r = Classify(a);
            Require(!r.empty(), "Warmup candidate missing");
            if (third) (void)mother::ValidatePrepared(a, r.front());
            else (void)paired::ValidatePrepared(a, r.front());
        }
    }
    unsigned successes = 0, skipped = 0;
    for (unsigned id = 301; id < 500; id += 2)
    {
        auto c = ReadDatasetCase(dir / CaseFilename(id)); CheckDataset500Contract(c, id);
        auto bound = BoundPair(*c.source, *c.input, c.refitted ? 1e-6 : 1e-9);
        Require(bound.passed, "Input independent source certificate failed");
        for (unsigned run = 0; run < 3; ++run)
        {
            std::optional<AnalyticPatch> first;
            for (unsigned order = 0; order < 2; ++order)
            {
                bool third = (order + run + id / 2) % 2 == 0;
                paired::AnalysisContext context(*c.input);
                std::vector<RecoveredCandidate> recovered;
                std::string status, message;
                auto start = std::chrono::steady_clock::now();
                try { recovered = Classify(context); }
                catch (const Failure &e) { status = "ClassificationFailed"; message = e.what(); }
                double classification = Ms(start);
                MotherDistanceValidationReport m; PairedValidationReport p; ValidationTiming t;
                if (recovered.empty())
                {
                    if (status.empty()) status = "NoCandidate";
                    Require(c.collapsedMask != 0, "Regular fixture lost classification"); ++skipped;
                }
                else
                {
                    if (first)
                        Require(SameMotherParameters(first->parameters, recovered.front().patch.parameters) &&
                            first->domain == recovered.front().patch.domain &&
                            first->mapping.swappedUV == recovered.front().patch.mapping.swappedUV &&
                            first->mapping.reversedU == recovered.front().patch.mapping.reversedU &&
                            first->mapping.reversedV == recovered.front().patch.mapping.reversedV,
                            "Compared validators received different candidates");
                    else first = recovered.front().patch;
                    if (third)
                    {
                        m = mother::ValidatePrepared(context, recovered.front()); t = m.timing;
                        status = PairedValidationStatusName(m.status); message = m.message;
                        Require(m.Passed(), "Mother regular coverage failed " + std::to_string(id) + ": " + message);
                        Require(m.sourceExtractions == 1 && m.sourceCacheHits == 0, "Mother benchmark warm source cache");
                        Spot(*c.input, recovered.front().patch, m);
                    }
                    else
                    {
                        p = paired::ValidatePrepared(context, recovered.front()); t = p.certificate.timing;
                        status = PairedValidationStatusName(p.status); message = p.certificate.message;
                        Require(p.Passed(), "Paired regular coverage failed " + std::to_string(id));
                        Require(p.sourceExtractions == 1 && p.sourceCacheHits == 0, "Paired benchmark warm source cache");
                    }
                    Timing(t); ++successes;
                }
                csv << id << ',' << TypeName(c.truth.type) << ',' << c.collapsedMask << ',' << run << ','
                    << (third ? "mother" : "paired") << ',' << recovered.size() << ',' << status << ','
                    << classification << ',' << t.totalMs << ',' << classification + t.totalMs << ','
                    << t.preparation.elapsedMs << ',' << t.surfaceBounds.elapsedMs << ','
                    << t.subdivision.elapsedMs << ',' << t.boundaryChecks.elapsedMs << ',' << t.reporting.elapsedMs << ','
                    << (third ? m.sourceExtractionMs : p.sourceExtractionMs) << ','
                    << (third ? m.sourceBoundaryMs : p.sourceBoundaryMs) << ','
                    << (third ? m.targetPreparationMs : p.nurbsMs + p.targetTransformMs + p.targetBoundaryMs) << ','
                    << (third ? m.mappingMs : p.mappingMs) << ',' << m.sampleMs << ',' << m.residualBuildMs << ','
                    << m.residualBoundMs << ',' << m.boundaryPositionMs << ',' << m.boundaryAngleMs << ','
                    << (third ? m.surfaceNodes : p.surfaceNodes) << ',' << m.boundaryPositionNodes << ','
                    << (third ? m.boundaryNodes : p.boundaryNodes) << ',' << m.surfaceDepth << ','
                    << (third ? m.mappingAttempts : p.mappingAttempts) << ','
                    << (third ? m.sourceCacheHits : p.sourceCacheHits) << ','
                    << (third ? m.motherDistance.upperBound : p.certificate.geometry.upperBound) << ','
                    << m.boundaryPosition.upperBound << ','
                    << (third ? m.boundaryTangent.upperBound : p.certificate.boundaryTangent.upperBound) << ','
                    << (third ? m.boundaryNormal.upperBound : p.certificate.boundaryNormal.upperBound) << ','
                    << bound.upper << ',' << Quote(message) << '\n';
                csv.flush();
            }
        }
        std::cout << "case " << id << " completed" << std::endl;
    }
    Require(successes == 540 && skipped == 60, "Sphere/torus coverage totals changed");
    std::cout << "Batch successes=" << successes << " no_candidate=" << skipped << " PASSED" << std::endl;
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init(); int result = 0;
    try
    {
        Require(argc >= 3, "Usage: --unit directory | --sphere-torus directory csv");
        if (std::string(argv[1]) == "--unit" && argc == 3) Unit(argv[2]);
        else if (std::string(argv[1]) == "--sphere-torus" && argc == 4) Batch(argv[2], argv[3]);
        else throw std::runtime_error("Invalid test arguments");
    }
    catch (const std::exception &e) { std::cerr << "FAIL " << e.what() << std::endl; result = 1; }
    sggk::fini(); return result;
}
