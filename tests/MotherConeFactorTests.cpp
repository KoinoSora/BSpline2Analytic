/**
 * @file MotherConeFactorTests.cpp
 * @brief 端点因子的位置模块回归及相同源缓存上的开关对照。
 * 合成面使用多项式勾股参数化；不调用退休验证器，不把位置通过视为完整转换成功。
 */
#include "MotherConeFactor.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::detail::mother;
namespace
{
unsigned checks = 0;
void Require(bool ok, const char *message)
{ ++checks; if (!ok) throw std::runtime_error(message); }

std::array<Poly, 3> Coordinates(unsigned p, unsigned q, double small, double noise = 0,
                                bool oscillate = false)
{
    std::array<Poly, 3> out{Poly(p, q), Poly(p, q), Poly(p, q)};
    for (unsigned i = 0; i <= p; ++i)
    {
        double t = double(i) / p;
        double w = .125 * (2 * t - 1);
        double ww = .015625 * (1 - 4 * t + 4. * i * (i ? i - 1 : 0) / (p * (p - 1)));
        for (unsigned j = 0; j <= q; ++j)
        {
            double g = small + (1 - small) * j / q;
            out[0].at(i, j) = 2 * w * g;
            out[1].at(i, j) = (1 - ww + (oscillate ? (i == 1 ? 3e-5 : i == 2 ? -3e-5 : 0) : 0)) * g;
            out[2].at(i, j) = (1 + ww) * g;
            if (j == 0) out[1].at(i, j) += noise;
        }
    }
    return out;
}
sggk::BSplineSurfacePtr Surface(const std::array<Poly, 3> &p, unsigned flags = 0,
    double scale = 1, V shift = {})
{
    sggk::Point3DMatrix cp(p[0].u + 1, std::vector<sggk::Point3D>(p[0].v + 1));
    for (unsigned i = 0; i <= p[0].u; ++i)
        for (unsigned j = 0; j <= p[0].v; ++j)
        {
            V x{p[0].at(i, j), p[1].at(i, j), p[2].at(i, j)};
            // 循环坐标旋转保持精确，另由候选轴和 xDirection 描述。
            if (flags & 8) x = {x.z, x.x, x.y};
            cp[i][j] = Point(x * scale + shift);
        }
    auto s = std::make_shared<sggk::BSplineSurface>(p[0].u, p[0].v, cp,
        sggk::RealArray{0, 1}, sggk::RealArray{0, 1},
        sggk::UIntArray{p[0].u + 1, p[0].u + 1}, sggk::UIntArray{p[0].v + 1, p[0].v + 1});
    if (flags & 1) s->UReverse();
    if (flags & 2) s->VReverse();
    if (flags & 4) s->ExchangeUV();
    return s;
}
AnalyticParameters Cone(unsigned flags = 0, V shift = {})
{
    AnalyticParameters p;
    p.type = sggk::SurfaceType::Cone; p.halfAngle = Pi / 4; p.origin = Point(shift);
    if (flags & 8) { p.axis = sggk::Dir3D(1, 0, 0); p.xDirection = sggk::Dir3D(0, 1, 0); }
    return p;
}
MotherDistanceValidationReport Position(SourceCache &source, const AnalyticParameters &p,
    bool enabled = true, unsigned depth = 12, unsigned nodes = 20000)
{
    MotherDistanceValidationReport r;
    RecognitionOptions o; o.maxSubdivisionDepth = depth; o.maxValidationPieces = nodes;
    ValidationTimer timer(r.timing);
    try { r.status = ValidatePosition(source, p, o, r, timer, nullptr, enabled)
        ? PairedValidationStatus::Passed : PairedValidationStatus::NotCertified; }
    catch (const paired::Budget &) { r.status = PairedValidationStatus::BudgetExceeded; }
    timer.Finish();
    return r;
}
SourceCache Prepare(const sggk::BSplineSurface &s)
{
    SourceCache source; MotherDistanceValidationReport r; source.Initialize(s, r); return source;
}
std::array<Poly, 3> CircularCoordinates(double small)
{
    // 独立幂基到 Bernstein 转换：五次正弦/四次余弦近似的圆弧截面。
    auto choose = [](unsigned n, unsigned k) {
        double x = 1; for (unsigned j = 1; j <= k; ++j) x *= double(n + 1 - j) / j; return x;
    };
    auto power = [&](unsigned i, unsigned k) {
        double value = 0;
        for (unsigned j = 0; j <= std::min(i, k); ++j)
            value += choose(k, j) * std::pow(-2., int(j)) * choose(i, j) / choose(5, j);
        return std::pow(-.125, int(k)) * value;
    };
    std::array<Poly, 3> p{Poly(5, 3), Poly(5, 3), Poly(5, 3)};
    for (unsigned i = 0; i <= 5; ++i)
        for (unsigned j = 0; j <= 3; ++j)
        {
            double g = small + (1 - small) * j / 3;
            p[0].at(i, j) = (1 - power(i, 2) / 2 + power(i, 4) / 24) * g;
            p[1].at(i, j) = (power(i, 1) - power(i, 3) / 6 + power(i, 5) / 120) * g;
            p[2].at(i, j) = g;
        }
    return p;
}
double PointDistance(V x, const AnalyticParameters &p)
{
    x = x - Get(p.origin);
    if (p.axis.X() == 1) x = {x.y, x.z, x.x};
    double rho = std::hypot(x.x, x.y), c = std::cos(p.halfAngle), s = std::sin(p.halfAngle);
    return rho * s + x.z * c >= 0 ? std::abs(rho * c - x.z * s) : std::hypot(rho, x.z);
}
void Audit(const sggk::BSplineSurface &s, const AnalyticParameters &p, double bound)
{
    // 独立 de Casteljau 在中心化后的实际存储控制点上求值，避免 SDK 世界坐标
    // 求值舍入掩盖大平移情况下很小的几何界；测试输入均为单个夹持 Bezier 单元。
    auto cp = s.ControlPoints();
    double magnitude = 0;
    for (const auto &row : cp) for (const auto &point : row)
        magnitude = std::max(magnitude, Norm(Get(point) - Get(p.origin)));
    auto evaluate = [](std::vector<V> a, double t) {
        for (size_t n = a.size(); n > 1; --n)
            for (size_t k = 0; k + 1 < n; ++k) a[k] = a[k] * (1 - t) + a[k + 1] * t;
        return a[0];
    };
    auto local = p; local.origin = sggk::Point3D(0, 0, 0);
    double round = 128 * std::numeric_limits<double>::epsilon() * magnitude;
    for (unsigned i = 0; i <= 16; ++i)
        for (unsigned j = 0; j <= 16; ++j)
        {
            std::vector<V> rows;
            for (const auto &row : cp)
            {
                std::vector<V> points;
                for (const auto &point : row) points.push_back(Get(point) - Get(p.origin));
                rows.push_back(evaluate(std::move(points), double(j) / 16));
            }
            V x = evaluate(std::move(rows), double(i) / 16);
            Require(PointDistance(x, local) <= bound + round, "Independent point exceeds cone position bound");
        }
}
void Unit()
{
    for (unsigned degree : {3u, 5u})
        for (double small : {0., 1e-12, 1e-4})
            for (unsigned flags = 0; flags < 16; ++flags)
            {
                V shift = flags & 8 ? V{1e6, -2e6, 3e6} : V{};
                double scale = flags & 8 ? 1000 : 1;
                auto s = Surface(Coordinates(degree, degree, small), flags, scale, shift);
                auto source = Prepare(*s); auto p = Cone(flags, shift);
                auto r = Position(source, p);
                Require(r.status == PairedValidationStatus::Passed, "Near-apex factor failed");
                Require(r.coneFactorAttempts <= 2, "Factor rebuilt on descendants");
                if (small <= 1e-12)
                    Require(r.coneFactorAcceptedNodes > 0, "Near-apex regression did not exercise factor");
                Audit(*s, p, r.motherDistance.upperBound);
            }
    auto s = Surface(Coordinates(5, 3, 1e-12)); auto source = Prepare(*s);
    Require(Position(source, Cone(), false).status == PairedValidationStatus::NotCertified,
        "Near-apex fixture does not reproduce original failure");
    auto good = Position(source, Cone());
    Require(good.status == PairedValidationStatus::Passed && good.surfaceNodes == 1,
        "Near-apex factor should certify the root");
    auto moved = Cone(); moved.origin = sggk::Point3D(.01, 0, 0);
    Require(Position(source, moved).status != PairedValidationStatus::Passed, "Stale candidate factor reused");
    Require(Position(source, Cone()).status == PairedValidationStatus::Passed, "Candidate cache recovery failed");
    for (double noise : {2e-6, 3e-5})
    {
        auto perturbed = Surface(Coordinates(5, 3, 1e-12, noise)); auto cache = Prepare(*perturbed);
        auto r = Position(cache, Cone());
        Require((r.status == PairedValidationStatus::Passed) == (noise < 1e-5), "Perturbation budget lost");
        if (noise < 1e-5)
        {
            Require(r.coneFactorPerturbation > 1e-6, "Removed endpoint perturbation not reported");
            Audit(*perturbed, Cone(), r.motherDistance.upperBound);
        }
        else
            Require(Position(cache, Cone(), true, 12, 1).status == PairedValidationStatus::BudgetExceeded,
                "Factor path bypassed node budget");
    }
    auto negative = Coordinates(3, 3, 0);
    for (double &z : negative[2].c) z = -z;
    auto neg = Surface(negative); auto negCache = Prepare(*neg);
    Require(Position(negCache, Cone()).status != PairedValidationStatus::Passed, "Negative cone sheet accepted");
    auto regular = Surface(Coordinates(5, 3, .1)); auto regularCache = Prepare(*regular);
    auto normal = Position(regularCache, Cone());
    Require(normal.status == PairedValidationStatus::Passed && normal.coneFactorAttempts == 0,
        "Ordinary passing cone triggered factor work");
    auto oscillating = Surface(Coordinates(3, 3, 0, 0, true)); auto osc = Prepare(*oscillating);
    auto split = Position(osc, Cone());
    Require(split.status == PairedValidationStatus::Passed && split.surfaceDepth > 0 &&
        split.coneFactorAttempts <= 2 && split.coneFactorAcceptedNodes > 0, "Factor subdivision regression");
    Audit(*oscillating, Cone(), split.motherDistance.upperBound);
    paired::Products products; MotherDistanceValidationReport diagnostics;
    auto coordinates = Coordinates(5, 3, 1e-4, 2e-6);
    auto factor = TryConeFactor(coordinates, 1e-12, Pi / 4, 1e-5, products, diagnostics, nullptr);
    Require(bool(factor), "Affine factor unavailable");
    auto children = SplitConeFactor(*factor, false);
    Require(children.second.factor.range().lo > .49, "Child lost original factor parameter");
    for (const auto *child : {&children.first, &children.second})
        Require(BoundConeFactor(*child).upper < 1e-5, "Split affine factor lost enclosure");
    diagnostics = {};
    Require(!TryConeFactor(coordinates, 2e-5, Pi / 4, 1e-5, products, diagnostics, nullptr),
        "Source envelope omitted from budget");
    auto tiny = *factor;
    tiny.z = BPoly(0, 0, I(-2e-7)); tiny.radial = BPoly(0, 0, I(1e-14));
    Require(BoundConeFactor(tiny).upper >= std::hypot(1e-7, 2e-7), "Apex bound missed negative height");
    ExecutionDeadline expired(1, ExecutionDeadline::Clock::now() - std::chrono::seconds(1));
    bool cancelled = false;
    try { (void)TryConeFactor(coordinates, 0, Pi / 4, 1e-5, products, diagnostics, &expired); }
    catch (const DeadlineExceeded &) { cancelled = true; }
    Require(cancelled, "Factor construction ignored deadline");
    AnalyticPatch apex; apex.parameters = Cone(); apex.domain = {0, .1, 0, 1};
    Require(ValidateMotherDistanceConversion(*s, apex).status == PairedValidationStatus::Unsupported,
        "Position change silently enabled singular chart contract");
    auto circular = Surface(CircularCoordinates(1e-4));
    auto patch = RecoverApproxPatch(*circular, Cone());
    auto complete = ValidateMotherDistanceConversion(*circular, patch);
    if (!complete.Passed()) std::cerr << "full validator: " << complete.message << '\n';
    Require(complete.Passed(), "Near-apex full position and joint tangent contract failed");
    Require(complete.coneFactorAcceptedNodes > 0, "Full validator did not exercise cone factor");
    auto wrong = patch; wrong.domain[0] += .01; wrong.domain[1] += .01;
    Require(!ValidateMotherDistanceConversion(*circular, wrong).Passed(), "Factor bypassed paired boundary checks");
    std::cout << "Cone factor checks=" << checks << " PASSED\n";
}
void Benchmark(const std::string &path)
{
    std::ofstream out(path); Require(bool(out), "Cannot open benchmark CSV");
    out << std::setprecision(17) << "case,run,enabled,passed,total_ms,position_ms,subdivision_ms,build_ms,factor_ms,nodes,attempts,accepted_nodes,bound\n";
    for (unsigned id = 0; id < 4; ++id)
    {
        auto s = Surface(Coordinates(5, 3, id == 0 ? .1 : 1e-12, id == 2 ? 2e-6 : id == 3 ? 3e-5 : 0));
        auto source = Prepare(*s);
        for (bool enabled : {false, true}) (void)Position(source, Cone(), enabled);
        for (unsigned run = 0; run < 100; ++run)
            for (unsigned order = 0; order < 2; ++order)
            {
                bool enabled = (run + order) % 2 != 0;
                auto r = Position(source, Cone(), enabled);
                out << id << ',' << run << ',' << enabled << ',' << (r.status == PairedValidationStatus::Passed)
                    << ',' << r.timing.totalMs << ',' << r.timing.surfaceBounds.elapsedMs
                    << ',' << r.timing.subdivision.elapsedMs << ',' << r.residualBuildMs << ',' << r.coneFactorBuildMs
                    << ',' << r.surfaceNodes << ',' << r.coneFactorAttempts << ',' << r.coneFactorAcceptedNodes
                    << ',' << r.motherDistance.upperBound << '\n';
            }
    }
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init(); int result = 0;
    try { if (argc == 3 && std::string(argv[1]) == "--benchmark") Benchmark(argv[2]); else Unit(); }
    catch (const std::exception &e) { std::cerr << "FAIL " << e.what() << '\n'; result = 1; }
    sggk::fini(); return result;
}
