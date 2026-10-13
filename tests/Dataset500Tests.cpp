#include "ApproximateInternal.h"
#include "Dataset500.h"
#include "IntervalBounds.h"
#include "Numerics.h"
#include "SurfaceConversion.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
using namespace bspline2analytic;
using namespace bspline2analytic::tests;
namespace
{
void Require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}
std::string Quote(std::string s)
{
    std::string out = "\"";
    for (char ch : s)
    {
        if (ch == '"')
            out += '"';
        out += ch == '\n' || ch == '\r' ? ' ' : ch;
    }
    return out + '"';
}
void Probe(const std::filesystem::path &dir)
{
    detail::BPoly x(1, 1), y(1, 1);
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j)
        {
            x.at(i, j) = detail::I(1 + 10 * j);
            y.at(i, j) = detail::I(1 + 10 * j) * detail::I(.1 + .2 * i);
        }
    auto angle = detail::AngleRange(x, y, 0);
    Require(angle.lo <= std::atan(.1) && angle.hi >= std::atan(.3) &&
                angle.hi - angle.lo < std::atan(.3) - std::atan(.1) + 1e-12,
            "Angle enclosure lost correlated radial scaling");
    Require(!detail::AngleRange(x * detail::I(-1), y, 0).finite(),
            "Angle bound accepted a nonpositive denominator");
    auto plane = ReadDatasetCase(dir / CaseFilename(0));
    auto same = BoundPair(*plane.input, *plane.input, 1e-9);
    Require(same.passed, "Identity pair certificate failed");
    auto translated = std::make_shared<sggk::BSplineSurface>(*plane.input);
    auto cp = translated->ControlPoints();
    for (unsigned i = 0; i < cp.size(); ++i)
        for (unsigned j = 0; j < cp[i].size(); ++j)
            translated->SetControlPoint(
                detail::Point(detail::Get(cp[i][j]) + detail::V{2e-6, 0, 0}), i, j);
    Require(!BoundPair(*plane.input, *translated, 1e-6).passed,
            "Pair verifier accepted excessive displacement");
    SurfaceConversionOptions bad;
    bad.polynomial.validation.distanceTolerance = -1;
    auto rejected = Convert(*plane.input, bad);
    Require(!rejected.Succeeded() && std::get<ConversionResult>(rejected.details).status ==
                                         ConversionStatus::InvalidInput,
            "Polynomial branch ignored options");
    auto cylinder = ReadDatasetCase(dir / CaseFilename(100));
    bad.rational.distanceTolerance = -1;
    rejected = Convert(*cylinder.input, bad);
    Require(!rejected.Succeeded() &&
                std::get<RecognitionResult>(rejected.details).status == Status::InvalidInput,
            "Rational branch ignored options");
}
void Inspect(const std::filesystem::path &file)
{
    using namespace detail;
    auto c = ReadDatasetCase(file);
    auto candidates = ClassifyApproximate(*c.input, ConversionOptions{});
    candidates.insert(candidates.begin(), RecoverApproxPatch(*c.input, c.truth));
    unsigned index = 0;
    for (const auto &p : candidates)
    {
        std::cout << "Candidate " << index++ << " (0 = known mother) "
                  << TypeName(p.parameters.type) << " radius=" << p.parameters.radius
                  << " half=" << p.parameters.halfAngle << std::endl;
        auto o = Get(p.parameters.origin);
        std::cout << "origin " << o.x << ' ' << o.y << ' ' << o.z << " domain ";
        for (double d : p.domain)
            std::cout << d << ' ';
        std::cout << std::endl;
        double sampled = 0;
        for (unsigned i = 0; i <= 16; ++i)
            for (unsigned j = 0; j <= 16; ++j)
            {
                V x = Get(c.input->CalcPoint(double(i) / 16, double(j) / 16));
                sampled = std::max(sampled, Norm(ProjectInteriorPoint(p, x) - x));
            }
        ValidationOptions options;
        options.maxValidationPieces = 500;
        auto r = ValidateConversion(*c.input, p, options);
        std::cout << "sampled=" << sampled << " geometry " << CheckStateName(r.geometry.state)
                  << " bound=" << r.geometry.upperBound
                  << " witness=" << r.geometry.maxSampledUpperBound
                  << " tangent=" << r.boundaryTangent.upperBound
                  << " normal=" << r.boundaryNormal.upperBound << " pieces=" << r.piecesVisited
                  << " region=";
        for (double v : r.geometry.parameterRegion)
            std::cout << v << ' ';
        std::cout << std::endl;
    }
}
int Run(const std::filesystem::path &dir, const std::filesystem::path &csvFile)
{
    unsigned files = 0;
    for (const auto &f : std::filesystem::directory_iterator(dir))
        if (f.path().extension() == ".b2a")
            ++files;
    Require(files == 500, "Dataset must contain exactly 500 cases");
    Probe(dir);
    std::ofstream csv(csvFile);
    Require(bool(csv), "Cannot open benchmark CSV");
    csv << std::setprecision(17);
    csv << "id,pair,type,source_kind,path,controls,degree_u,degree_v,spans_u,spans_v,collapsed_"
           "mask,fixture_bound,status,classify_ms,optimize_ms,validate_ms,total_ms,geometry_bound,"
           "tangent_deg,normal_deg,message\n";
    unsigned bins[5][5][2]{}, collapsed[5]{}, success = 0, degenerate = 0, failed = 0, rational = 0,
                                              polynomial = 0;
    bool warmed[5][2]{};
    std::map<std::string, std::vector<double>> times;
    for (unsigned id = 0; id < 500; ++id)
    {
        auto c = ReadDatasetCase(dir / CaseFilename(id));
        Require(c.id == id && c.pair == id / 2 && c.type == id / 100 &&
                    c.bucket == (id % 100) / 20 && c.variant == (id % 20) / 2 &&
                    c.refitted == bool(id % 2) && c.seed == 0x500b2a17u + c.pair,
                "Dataset metadata mismatch");
        unsigned controls = static_cast<unsigned>(c.input->ControlPoints().size() *
                                                  c.input->ControlPoints()[0].size());
        const unsigned lo[] = {50, 80, 110, 140, 170}, hi[] = {79, 109, 139, 169, 200};
        Require(controls >= lo[c.bucket] && controls <= hi[c.bucket],
                "Control count outside allocated band");
        Require(!c.refitted || !c.input->IsRational(),
                "Refitted fixture still has rational weights");
        unsigned expectedCollapsed = (c.type == 2 || c.type == 3) && c.variant >= 8 ? 1 : 0;
        Require(bool(c.collapsedMask) == bool(expectedCollapsed) &&
                    CollapsedMask(*c.input) == c.collapsedMask,
                "Collapsed boundary mismatch");
        // The generated boundary must be a point, not merely a short curve.
        const auto &net = c.input->ControlPoints();
        for (unsigned edge = 0; edge < 4; ++edge)
            if (c.collapsedMask & (1u << edge))
            {
                const size_t nu = net.size(), nv = net[0].size();
                auto origin = detail::Get(net[edge == 1 ? nu - 1 : 0][edge == 3 ? nv - 1 : 0]);
                for (size_t j = 0; j < (edge < 2 ? nv : nu); ++j)
                {
                    auto point = detail::Get(net[edge < 2 ? (edge == 0 ? 0 : nu - 1) : j]
                                                [edge < 2 ? j : (edge == 2 ? 0 : nv - 1)]);
                    Require(point.x == origin.x && point.y == origin.y && point.z == origin.z,
                            "Degenerate control row is not exactly collapsed");
                }
            }
        auto source = RebuildDatasetSource(c.type, c.bucket, c.variant);
        Require(BoundPair(*source, *c.source, 1e-9).passed,
                "Stored analytic source does not match generator");
        auto bound = BoundPair(*c.source, *c.input, c.refitted ? 1e-6 : 1e-9);
        Require(bound.passed && std::isfinite(c.storedBound) &&
                    c.storedBound < (c.refitted ? 1e-6 : 1e-9),
                "Fixture distance not certified");
        ++bins[c.type][c.bucket][c.refitted];
        if (c.collapsedMask)
            ++collapsed[c.type];
        unsigned path = c.input->IsRational() ? 0 : 1;
        if (!warmed[c.type][path] && !c.collapsedMask)
        {
            (void)Convert(*c.input);
            warmed[c.type][path] = true;
        }
        auto result = Convert(*c.input);
        Require(result.path == (path == 0 ? ConversionPath::RationalRecognition
                                          : ConversionPath::PolynomialApproximation),
                "Dispatch path mismatch");
        if (path == 0)
            ++rational;
        else
            ++polynomial;
        bool ok = c.collapsedMask ? result.IsDegenerate()
                                  : (result.Succeeded() && result.Type() == c.truth.type);
        if (ok)
        {
            if (c.collapsedMask)
                ++degenerate;
            else
                ++success;
        }
        else
        {
            ++failed;
            std::cout << "FAIL " << id << ' ' << result.StatusText() << ' ' << result.Message()
                      << std::endl;
        }
        const auto &t = result.Timing();
        Require(t.classification.elapsedMs >= 0 && t.optimization.elapsedMs >= 0 &&
                    t.validation.elapsedMs >= 0 && t.totalMs >= 0,
                "Negative stage timing");
        double geometry = INFINITY, tangent = INFINITY, normal = INFINITY;
        if (const auto *r = std::get_if<RecognitionResult>(&result.details))
        {
            geometry = r->validation.surfaceDistanceBound;
            Require(t.optimization.calls == 0, "Exact branch must not claim optimization");
            if (r->status == Status::Success)
                Require(r->validation.level == VerificationLevel::BernsteinBound,
                        "Exact branch lost verification contract");
        }
        else
        {
            const auto &approx = std::get<ConversionResult>(result.details);
            geometry = approx.validation.geometry.upperBound;
            tangent = approx.validation.boundaryTangent.upperBound;
            normal = approx.validation.boundaryNormal.upperBound;
            if (approx.status == ConversionStatus::Success)
                Require(approx.validation.Passed() && geometry < 1e-5 && tangent < .1 &&
                            normal < .1,
                        "Approximate success violates constraints");
        }
        csv << id << ',' << c.pair << ',' << TypeName(c.truth.type) << ','
            << (c.refitted ? "refitted" : "analytic") << ',' << ConversionPathName(result.path)
            << ',' << controls << ',' << c.input->DegreeU() << ',' << c.input->DegreeV() << ','
            << c.input->KnotsU().size() - 1 << ',' << c.input->KnotsV().size() - 1 << ','
            << c.collapsedMask << ',' << bound.upper << ',' << result.StatusText() << ','
            << t.classification.elapsedMs << ',' << t.optimization.elapsedMs << ','
            << t.validation.elapsedMs << ',' << t.totalMs << ',' << geometry << ',';
        if (path)
            csv << tangent;
        csv << ',';
        if (path)
            csv << normal;
        csv << ',' << Quote(result.Message()) << '\n';
        csv.flush();
        if (!c.collapsedMask)
            times[std::string(TypeName(c.truth.type)) + "/" + ConversionPathName(result.path)]
                .push_back(t.classification.elapsedMs);
        if (id % 20 == 19)
            std::cout << "Tested " << id + 1 << "/500; success=" << success
                      << " degenerate=" << degenerate << " failures=" << failed << std::endl;
    }
    for (unsigned type = 0; type < 5; ++type)
    {
        Require(collapsed[type] == ((type == 2 || type == 3) ? 20u : 0u),
                "Degenerate quota mismatch");
        for (unsigned b = 0; b < 5; ++b)
            for (unsigned k = 0; k < 2; ++k)
                Require(bins[type][b][k] == 10, "Dataset quota mismatch");
    }
    Require(rational == 200 && polynomial == 300, "Unexpected dispatch distribution");
    for (auto &[name, values] : times)
    {
        std::sort(values.begin(), values.end());
        std::cout << name << " classify mean="
                  << std::accumulate(values.begin(), values.end(), 0.0) / values.size()
                  << " p95=" << values[(values.size() * 95 + 99) / 100 - 1] << " ms\n";
    }
    std::cout << "500 fixtures certified; regular=" << success << " degenerate=" << degenerate
              << " unexpected=" << failed << "; rational=" << rational
              << " polynomial=" << polynomial << std::endl;
    Require(bool(csv), "CSV write failed");
    return failed ? 1 : 0;
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init();
    int result = 1;
    try
    {
        Require(argc >= 3, "Usage: Dataset500Tests --generate directory | directory results.csv");
        if (std::string(argv[1]) == "--inspect")
        {
            Inspect(argv[2]);
            result = 0;
        }
        else if (std::string(argv[1]) == "--generate")
        {
            GenerateDataset500(argv[2]);
            result = 0;
        }
        else
            result = Run(argv[1], argv[2]);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
    }
    sggk::fini();
    return result;
}
