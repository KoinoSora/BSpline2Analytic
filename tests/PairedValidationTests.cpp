#include "ApproximateCases.h"
#include "Dataset500.h"
#include "PairedValidationInternal.h"
#include <Foundation/init.h>
#include <cfenv>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace
{
unsigned checks = 0;
void Require(bool b, const std::string &s)
{
    ++checks;
    if (!b)
        throw std::runtime_error(s);
}
std::string Quote(std::string s)
{
    std::string q = "\"";
    for (char c : s)
    {
        if (c == '\"')
            q += '\"';
        q += c;
    }
    return q + '\"';
}
void Timing(const PairedValidationReport &r)
{
    const auto &t = r.certificate.timing;
    double sum = t.preparation.elapsedMs + t.surfaceBounds.elapsedMs + t.subdivision.elapsedMs +
                 t.boundaryChecks.elapsedMs + t.reporting.elapsedMs;
    Require(std::isfinite(t.totalMs) && t.totalMs >= 0 &&
                std::abs(sum - t.totalMs) < 1e-7 * std::max(1.0, t.totalMs),
            "Exclusive paired timing mismatch");
    Require(r.nurbsMs + r.mappingMs <= t.preparation.elapsedMs + 1e-5,
            "Nested preparation timings exceed preparation");
    if (r.Passed())
        Require(r.certificate.Passed() && std::isfinite(r.certificate.geometry.upperBound),
                "Incomplete certificate accepted");
}
sggk::BSplineSurfacePtr Plane(unsigned mode)
{
    sggk::Point3DMatrix cp(4, std::vector<sggk::Point3D>(4));
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
        {
            // mode 4: x=(u+u^2)/2, represented exactly in degree three.
            double x =
                mode == 4 ? (double(i) / 6 + double(i * (i ? i - 1 : 0)) / 12) : double(i) / 3;
            if (mode == 5)
                x = double(i) / 3 +
                    .1 * (double(i) / 3 - double(i * (i ? i - 1 : 0)) / 6) * double(j) / 3;
            cp[i][j] = sggk::Point3D(x * (mode >= 4 ? 1 : .001), double(j) / 3, 0);
        }
    if (mode == 1)
        for (unsigned j = 1; j < 3; ++j)
            cp[1][j] = sggk::Point3D(.001 / 3, double(j) / 3, 8e-6);
    if (mode == 2)
        for (unsigned i = 1; i < 3; ++i)
            cp[i][0] = sggk::Point3D(double(i) * .001 / 3, 8e-6, 0);
    if (mode == 3)
        for (unsigned i = 1; i < 3; ++i)
            for (unsigned j = 1; j < 3; ++j)
                cp[i][j] = sggk::Point3D(double(i) * .001 / 3, double(j) / 3, 1e-3);
    return std::make_shared<sggk::BSplineSurface>(3, 3, cp, sggk::RealArray{0, 1},
                                                  sggk::RealArray{0, 1}, sggk::UIntArray{4, 4},
                                                  sggk::UIntArray{4, 4});
}
void Unit()
{
    AnalyticPatch patch;
    patch.parameters.type = sggk::SurfaceType::Plane;
    patch.domain = {0, .001, 0, 1};
    MaterializePatch(patch);
    auto exact = ValidatePairedConversion(*Plane(0), patch);
    Timing(exact);
    Require(exact.Passed(), "Exact plane: " + exact.certificate.message);
    for (unsigned i = 1; i < 4; ++i)
    {
        PairedValidationOptions o;
        o.enablePiecewiseMapping = false;
        auto r = ValidatePairedConversion(*Plane(i), patch, o);
        Timing(r);
        Require(!r.Passed(), "Missed distance/tangent/normal violation " + std::to_string(i));
        Require(i != 1 || r.certificate.boundaryNormal.state == CheckState::Exceeded,
                "Normal-only violation not detected");
        Require(i != 2 || r.certificate.boundaryTangent.state == CheckState::Exceeded,
                "Tangent-only violation not detected");
        Require(i != 3 || r.certificate.geometry.state == CheckState::Exceeded,
                "Interior bulge not detected");
    }
    auto warped = Plane(4);
    patch.domain = {0, 1, 0, 1};
    MaterializePatch(patch);
    PairedValidationOptions affine;
    affine.enablePiecewiseMapping = false;
    auto noMap = ValidatePairedConversion(*warped, patch, affine);
    Require(!noMap.Passed(), "Affine correspondence should fail on warped plane");
    auto mapped = ValidatePairedConversion(*warped, patch);
    Timing(mapped);
    std::cout << "warped " << PairedValidationStatusName(mapped.status)
              << " segments=" << mapped.uMap.source.size() - 1 << " " << mapped.certificate.message
              << std::endl;
    Require(mapped.Passed() && mapped.mappingAttempts == 2 && mapped.uMap.source.size() > 2,
            "Piecewise mapping failed on warped plane");
    for (unsigned flags = 1; flags < 8; ++flags)
    {
        auto oriented = Plane(4);
        if (flags & 1)
            oriented->UReverse();
        if (flags & 2)
            oriented->VReverse();
        if (flags & 4)
            oriented->ExchangeUV();
        oriented = NormalizeDatasetSurface(*oriented);
        auto p = RecoverApproxPatch(*oriented, patch.parameters);
        auto r = ValidatePairedConversion(*oriented, p);
        Timing(r);
        Require(r.Passed(), "Reversed/swapped piecewise map failed");
    }
    auto differentEdges = Plane(5);
    Require(ValidateConversion(*differentEdges, patch).Passed(),
            "Different-speed plane is not a valid geometric fixture");
    auto conservative = ValidatePairedConversion(*differentEdges, patch);
    Timing(conservative);
    Require(!conservative.Passed() && conservative.status == PairedValidationStatus::NotCertified,
            "Incompatible separable map was incorrectly certified");
    // A second call must use changed geometry/domain, not a stale successful map.
    auto shifted = patch;
    shifted.parameters.origin = sggk::Point3D(0, 0, 1e-3);
    MaterializePatch(shifted);
    Require(!ValidatePairedConversion(*warped, shifted, affine).Passed(),
            "Stale candidate geometry accepted");
    auto smaller = patch;
    smaller.domain[1] = .9;
    MaterializePatch(smaller);
    Require(!ValidatePairedConversion(*warped, smaller, affine).Passed(),
            "Stale candidate domain accepted");
    auto flat = Plane(0);
    auto optimized = RecoverApproxPatch(*flat, patch.parameters);
    optimized.parameters.origin = sggk::Point3D(0, 0, 1e-3);
    MaterializePatch(optimized);
    auto samples = SampleSurface(*flat, {});
    ConversionStatistics optimization;
    Require(!ValidatePairedConversion(*flat, optimized, affine).Passed(),
            "Displaced candidate unexpectedly passed before refinement");
    RefineCandidate(*flat, optimized, samples, {}, optimization);
    auto afterOptimization = ValidatePairedConversion(*flat, optimized);
    Timing(afterOptimization);
    Require(optimization.optimizationIterations > 0 && afterOptimization.Passed(),
            "Independent verification of an actually optimized candidate failed");
    PairedValidationOptions budget;
    budget.validation.maxValidationPieces = 1;
    auto limited = ValidatePairedConversion(
        *Plane(0), RecoverApproxPatch(*Plane(0), patch.parameters), budget);
    Timing(limited);
    Require(limited.status == PairedValidationStatus::BudgetExceeded && !limited.Passed(),
            "Node budget not enforced");
    budget = {};
    budget.maxMappingSegments = 2;
    Require(ValidatePairedConversion(*warped, patch, budget).status ==
                PairedValidationStatus::BudgetExceeded,
            "Mapping budget not enforced");
    auto wrong = patch;
    wrong.normalSign = -1;
    Require(!ValidatePairedConversion(*warped, wrong, affine).Passed(),
            "Opposite normals accepted");
    PairedValidationOptions bad;
    bad.validation.distanceTolerance = 0;
    Require(ValidatePairedConversion(*warped, patch, bad).status ==
                PairedValidationStatus::InvalidInput,
            "Invalid tolerance accepted");
    int rounding = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    auto roundingResult = ValidatePairedConversion(*warped, patch);
    std::fesetround(rounding);
    Require(!roundingResult.Passed(), "Unsupported rounding mode accepted");
    // Verify each analytic rational chart against direct analytic evaluation,
    // including multi-span arcs, swaps/reversals and both angular directions.
    for (auto type : {sggk::SurfaceType::Plane, sggk::SurfaceType::Cylinder,
                      sggk::SurfaceType::Cone, sggk::SurfaceType::Sphere, sggk::SurfaceType::Torus})
        for (unsigned flags = 0; flags < 8; ++flags)
        {
            AnalyticPatch p;
            p.parameters.type = type;
            p.parameters.radius = 2;
            p.parameters.majorRadius = 4;
            p.parameters.halfAngle = .4;
            p.domain = {6.1, 8.3, .4, 2.1};
            p.mapping.swappedUV = bool(flags & 1);
            p.mapping.reversedU = bool(flags & 2);
            p.mapping.reversedV = bool(flags & 4);
            auto chart = paired::BuildChart(p);
            for (const auto &c : chart.cells)
            {
                double u = .31, v = .67, ou = c.surface.u0 + (c.surface.u1 - c.surface.u0) * u,
                       ov = c.surface.v0 + (c.surface.v1 - c.surface.v0) * v;
                double coordinates[2] = {p.mapping.reversedU ? 1 - ou : ou,
                                         p.mapping.reversedV ? 1 - ov : ov};
                if (p.mapping.swappedUV)
                    std::swap(coordinates[0], coordinates[1]);
                double angles[2]{};
                for (unsigned d = 0; d < 2; ++d)
                {
                    double scaled = coordinates[d] * (chart.angles[d].size() - 1);
                    unsigned k = static_cast<unsigned>(scaled);
                    double t = scaled - k, a = chart.angles[d][k], b = chart.angles[d][k + 1];
                    angles[d] = chart.angular[d]
                                    ? a + 2 * std::atan2(t * std::sin((b - a) * .5),
                                                         1 - t + t * std::cos((b - a) * .5))
                                    : a + t * (b - a);
                    Require(std::abs(paired::RationalCoordinate(chart, d, angles[d]) -
                                     coordinates[d]) < 1e-13,
                            "Angular inverse mismatch");
                }
                IV num = Range(Evaluate(c.surface.p, u, v));
                I w = Evaluate(c.surface.w, u, v).range();
                V expected = CanonicalPoint(p.parameters, angles[0], angles[1]);
                for (auto pair :
                     {std::make_pair(num.x / w, expected.x), std::make_pair(num.y / w, expected.y),
                      std::make_pair(num.z / w, expected.z)})
                    Require(pair.first.lo - 1e-12 <= pair.second &&
                                pair.second <= pair.first.hi + 1e-12,
                            "Rational chart evaluation mismatch");
            }
        }
    // Exercise the entire common-partition certificate across rational chart
    // seams, including a periodic seam and a swapped/reversed source axis.
    AnalyticPatch cylinder;
    cylinder.parameters.type = sggk::SurfaceType::Cylinder;
    cylinder.parameters.radius = 2;
    cylinder.domain = {6.1, 8.3, -.4, 1.2};
    MaterializePatch(cylinder);
    auto rational = cylinder.surface->ToBSpline(
        sggk::UVRange(sggk::Interval(6.1, 8.3), sggk::Interval(-.4, 1.2)));
    Require(bool(rational), "Multi-span reference creation failed");
    for (unsigned flags : {0u, 7u})
    {
        auto polynomial = RefitPolynomial(*rational, 12);
        if (flags & 1)
            polynomial->UReverse();
        if (flags & 2)
            polynomial->VReverse();
        if (flags & 4)
            polynomial->ExchangeUV();
        polynomial = NormalizeDatasetSurface(*polynomial);
        auto recovered = RecoverApproxPatch(*polynomial, cylinder.parameters);
        auto result = ValidatePairedConversion(*polynomial, recovered);
        Timing(result);
        Require(result.Passed() && result.targetPatches > 1,
                "Multi-span cylinder certificate failed: " + result.certificate.message);
    }
    std::cout << "Paired unit assertions: " << checks << std::endl;
}
void SpotCheck(const sggk::BSplineSurface &s, const AnalyticPatch &p,
               const PairedValidationReport &r)
{
    auto chart = paired::BuildChart(p);
    auto frame = FrameOf(p.parameters);
    auto map = [](const PairedParameterMap &m, double t) {
        if (t == m.source.back())
            return m.target.back();
        size_t i = static_cast<size_t>(std::upper_bound(m.source.begin(), m.source.end(), t) -
                                       m.source.begin() - 1);
        return m.target[i] + (m.target[i + 1] - m.target[i]) * (t - m.source[i]) /
                                 (m.source[i + 1] - m.source[i]);
    };
    // Independent double evaluation checks the certificate, never establishes it.
    SurfaceSampleContext cache(s);
    for (unsigned i = 0; i <= 8; ++i)
        for (unsigned j = 0; j <= 8; ++j)
        {
            double u =
                i == 8 ? s.MaxParamU() : s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 8;
            double v =
                j == 8 ? s.MaxParamV() : s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * j / 8;
            double q[2] = {map(r.uMap, u), map(r.vMap, v)};
            if (p.mapping.reversedU)
                q[0] = 1 - q[0];
            if (p.mapping.reversedV)
                q[1] = 1 - q[1];
            if (p.mapping.swappedUV)
                std::swap(q[0], q[1]);
            double angle[2];
            for (unsigned d = 0; d < 2; ++d)
            {
                unsigned n = static_cast<unsigned>(chart.angles[d].size() - 1);
                unsigned k = std::min(n - 1, static_cast<unsigned>(std::max(0.0, q[d]) * n));
                double t = q[d] * n - k, a = chart.angles[d][k], b = chart.angles[d][k + 1];
                angle[d] = chart.angular[d] ? a + 2 * std::atan2(t * std::sin((b - a) / 2),
                                                                 1 - t + t * std::cos((b - a) / 2))
                                            : a + t * (b - a);
            }
            V target = frame.global(CanonicalPoint(p.parameters, angle[0], angle[1]));
            const auto &sample = cache.At(u, v);
            Require(Norm(sample.point - target) <= r.certificate.geometry.upperBound + 1e-10,
                    "Sample exceeds paired position certificate");
            if (i == 0 || i == 8 || j == 0 || j == 8)
            {
                V n = frame.vector(CanonicalNormal(p.parameters, angle[0], angle[1])) *
                      double(p.normalSign);
                Require(Angle(sample.normal, n) * 180 / Pi <=
                            r.certificate.boundaryNormal.upperBound + 1e-7,
                        "Sample exceeds paired normal certificate");
            }
        }
}
int Dataset(const std::filesystem::path &dir, const std::filesystem::path &file, bool sphere)
{
    std::ofstream csv(file);
    Require(bool(csv), "Cannot open paired CSV");
    csv << std::setprecision(17);
    csv << "id,type,rational,collapsed,controls,classify_ms,candidates,old_status,new_status,old_"
           "ms,new_ms,speedup,geometry_bound,tangent_deg,normal_deg,attempts,u_segments,v_segments,"
           "common_cells,surface_nodes,boundary_nodes,prepare_ms,surface_ms,subdivide_ms,boundary_"
           "ms,report_ms,nurbs_ms,mapping_ms,message\n";
    unsigned passed = 0, notCertified = 0, skipped = 0;
    bool warmed[5]{};
    for (unsigned id = sphere ? 301u : 0u; id < (sphere ? 400u : 500u); id += sphere ? 2u : 1u)
    {
        auto c = ReadDatasetCase(dir / CaseFilename(id));
        Require(c.id == id, "Fixture identity mismatch");
        auto inputBound = BoundPair(*c.source, *c.input, c.refitted ? 1e-6 : 1e-9);
        Require(inputBound.passed, "Fixture source certificate failed");
        unsigned controls = static_cast<unsigned>(c.input->ControlPoints().size() *
                                                  c.input->ControlPoints()[0].size());
        std::string oldStatus = "NotRun", newStatus = "NotRun", message;
        ConversionValidationReport old;
        PairedValidationReport newer;
        double classifyMs = 0, oldMs = 0, newMs = 0;
        size_t candidates = 0;
        if (c.input->IsRational())
        {
            // New module intentionally supports polynomial input only.
            AnalyticPatch p;
            p.parameters = c.truth;
            p.domain = c.analyticDomain;
            auto r = ValidatePairedConversion(*c.input, p);
            Require(r.status == PairedValidationStatus::Unsupported,
                    "Rational input was unexpectedly accepted");
            newStatus = "Unsupported";
            message = r.certificate.message;
            ++skipped;
        }
        else
        {
            try
            {
                auto start = std::chrono::steady_clock::now();
                auto list = ClassifyApproximate(*c.input, ConversionOptions{});
                classifyMs = std::chrono::duration<double, std::milli>(
                                 std::chrono::steady_clock::now() - start)
                                 .count();
                candidates = list.size();
                if (list.empty())
                {
                    oldStatus = newStatus = "NoCandidate";
                    ++skipped;
                }
                else
                {
                    // First actual recovered candidate; never substitute fixture truth.
                    const auto &p = list.front();
                    if (!warmed[c.type] && !c.collapsedMask)
                    {
                        (void)ValidateConversion(*c.input, p);
                        (void)ValidatePairedConversion(*c.input, p);
                        warmed[c.type] = true;
                    }
                    auto oldRun = [&] {
                        auto t = std::chrono::steady_clock::now();
                        old = ValidateConversion(*c.input, p);
                        oldMs = std::chrono::duration<double, std::milli>(
                                    std::chrono::steady_clock::now() - t)
                                    .count();
                    };
                    auto newRun = [&] {
                        auto t = std::chrono::steady_clock::now();
                        newer = ValidatePairedConversion(*c.input, p);
                        newMs = std::chrono::duration<double, std::milli>(
                                    std::chrono::steady_clock::now() - t)
                                    .count();
                    };
                    if (c.pair % 2)
                    {
                        newRun();
                        oldRun();
                    }
                    else
                    {
                        oldRun();
                        newRun();
                    }
                    Timing(newer);
                    oldStatus = old.Passed() ? "Passed" : "NotCertified";
                    newStatus = PairedValidationStatusName(newer.status);
                    message = newer.certificate.message;
                    if (newer.Passed())
                    {
                        Require(!c.collapsedMask, "Degenerate input unexpectedly certified");
                        Require(newer.certificate.geometry.upperBound < 1e-5 &&
                                    newer.certificate.boundaryTangent.upperBound < .1 &&
                                    newer.certificate.boundaryNormal.upperBound < .1,
                                "Tolerance violated");
                        SpotCheck(*c.input, p, newer);
                        ++passed;
                    }
                    else
                        ++notCertified;
                }
            }
            catch (const Failure &e)
            {
                oldStatus = newStatus = "ClassificationFailed";
                message = e.what();
                ++skipped;
            }
        }
        const auto &r = newer.certificate;
        const auto &t = r.timing;
        csv << id << ',' << TypeName(c.truth.type) << ',' << c.input->IsRational() << ','
            << c.collapsedMask << ',' << controls << ',' << classifyMs << ',' << candidates << ','
            << oldStatus << ',' << newStatus << ',' << oldMs << ',' << newMs << ','
            << (newMs ? oldMs / newMs : 0) << ',' << r.geometry.upperBound << ','
            << r.boundaryTangent.upperBound << ',' << r.boundaryNormal.upperBound << ','
            << newer.mappingAttempts << ','
            << (newer.uMap.source.empty() ? 0 : newer.uMap.source.size() - 1) << ','
            << (newer.vMap.source.empty() ? 0 : newer.vMap.source.size() - 1) << ','
            << newer.commonCells << ',' << newer.surfaceNodes << ',' << newer.boundaryNodes << ','
            << t.preparation.elapsedMs << ',' << t.surfaceBounds.elapsedMs << ','
            << t.subdivision.elapsedMs << ',' << t.boundaryChecks.elapsedMs << ','
            << t.reporting.elapsedMs << ',' << newer.nurbsMs << ',' << newer.mappingMs << ','
            << Quote(message) << '\n';
        csv.flush();
        std::cout << id << " old=" << oldStatus << " new=" << newStatus << " ms=" << newMs
                  << std::endl;
    }
    std::cout << "Paired comparison passed=" << passed << " not_certified=" << notCertified
              << " skipped=" << skipped << std::endl;
    // Every regular polynomial fixture has an independently certified source.
    // Guard coverage as well as the correctness of individual successes.
    Require(passed == (sphere ? 40u : 280u), "Regular polynomial dataset coverage regressed");
    return 0;
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init();
    int code = 0;
    try
    {
        Require(argc == 2 || argc == 4, "Expected --unit or --sphere/--dataset directory csv");
        if (std::string(argv[1]) == "--unit")
        {
            Require(argc == 2, "Unexpected unit arguments");
            Unit();
        }
        else
        {
            Require(argc == 4 &&
                        (std::string(argv[1]) == "--sphere" || std::string(argv[1]) == "--dataset"),
                    "Invalid comparison arguments");
            code = Dataset(argv[2], argv[3], std::string(argv[1]) == "--sphere");
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << std::endl;
        code = 1;
    }
    sggk::fini();
    return code;
}
