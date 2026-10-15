#include "Dataset500.h"
#include "ApproximateInternal.h"
#include <Geometry/3D/Surface/ConeSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <Geometry/3D/Surface/TorusSurface.h>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
namespace bspline2analytic::tests
{
using namespace detail;
sggk::SurfaceType DatasetType(unsigned i)
{
    const sggk::SurfaceType types[] = {sggk::SurfaceType::Plane, sggk::SurfaceType::Cylinder,
                                       sggk::SurfaceType::Cone, sggk::SurfaceType::Sphere,
                                       sggk::SurfaceType::Torus};
    if (i >= 5)
        throw std::runtime_error("Invalid dataset type");
    return types[i];
}
std::string CaseFilename(unsigned id)
{
    std::ostringstream s;
    s << "case-" << std::setfill('0') << std::setw(3) << id << ".b2a";
    return s.str();
}
namespace
{
void Require(bool ok, const std::string &s)
{
    if (!ok)
        throw std::runtime_error(s);
}
double Number(unsigned seed, unsigned channel)
{
    uint32_t x = seed + 0x9e3779b9u * (channel + 1);
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return double(x) / 4294967296.0;
}
void WriteV(std::ostream &o, V v)
{
    o << v.x << ' ' << v.y << ' ' << v.z << '\n';
}
V ReadV(std::istream &in)
{
    V v;
    in >> v.x >> v.y >> v.z;
    return v;
}
void WriteSurface(std::ostream &o, const sggk::BSplineSurface &s)
{
    const auto &cp = s.ControlPoints();
    const auto &w = s.Weights();
    o << s.DegreeU() << ' ' << s.DegreeV() << ' ' << cp.size() << ' ' << cp.front().size() << ' '
      << !w.empty() << '\n';
    auto knots = [&](const auto &k, const auto &m) {
        o << k.size() << '\n';
        for (size_t i = 0; i < k.size(); ++i)
            o << k[i] << ' ' << m[i] << '\n';
    };
    knots(s.KnotsU(), s.MultsU());
    knots(s.KnotsV(), s.MultsV());
    for (size_t i = 0; i < cp.size(); ++i)
        for (size_t j = 0; j < cp[i].size(); ++j)
        {
            V p = Get(cp[i][j]);
            o << p.x << ' ' << p.y << ' ' << p.z;
            if (!w.empty())
                o << ' ' << w[i][j];
            o << '\n';
        }
}
sggk::BSplineSurfacePtr ReadSurface(std::istream &in)
{
    unsigned p = 0, q = 0, nu = 0, nv = 0;
    bool rational = false;
    in >> p >> q >> nu >> nv >> rational;
    Require(in.good() && nu && nv && nu <= 200 && nv <= 200 && nu * nv <= 1000 && p < nu && q < nv,
            "Invalid surface dimensions");
    sggk::RealArray ku, kv;
    sggk::UIntArray mu, mv;
    auto knots = [&](auto &k, auto &m, unsigned degree, unsigned count) {
        unsigned n = 0;
        in >> n;
        Require(n >= 2 && n <= 201, "Invalid knot count");
        k.resize(n);
        m.resize(n);
        unsigned total = 0;
        for (unsigned i = 0; i < n; ++i)
        {
            in >> k[i] >> m[i];
            Require(std::isfinite(k[i]) && m[i] > 0 && m[i] <= degree + 1 &&
                        (!i || k[i] > k[i - 1]),
                    "Invalid knot data");
            total += m[i];
        }
        Require(k.front() == 0 && k.back() == 1 && m.front() == degree + 1 &&
                    m.back() == degree + 1 && total == count + degree + 1,
                "Invalid clamped knot vector");
    };
    knots(ku, mu, p, nu);
    knots(kv, mv, q, nv);
    sggk::Point3DMatrix cp(nu, std::vector<sggk::Point3D>(nv));
    sggk::RealMatrix w;
    if (rational)
        w.assign(nu, std::vector<double>(nv));
    for (unsigned i = 0; i < nu; ++i)
        for (unsigned j = 0; j < nv; ++j)
        {
            V v = ReadV(in);
            Require(Finite(v), "Invalid control point");
            cp[i][j] = Point(v);
            if (rational)
            {
                in >> w[i][j];
                Require(std::isfinite(w[i][j]) && w[i][j] > 0, "Invalid weight");
            }
        }
    Require(!in.fail(), "Truncated surface");
    if (rational)
        return std::make_shared<sggk::BSplineSurface>(p, q, cp, w, ku, kv, mu, mv);
    return std::make_shared<sggk::BSplineSurface>(p, q, cp, ku, kv, mu, mv);
}
void Save(const DatasetCase &c, const std::filesystem::path &file)
{
    std::ofstream o(file);
    Require(bool(o), "Cannot write dataset case");
    o << std::setprecision(17);
    o << "B2A_CASE_V1\n"
      << c.id << ' ' << c.pair << ' ' << c.type << ' ' << c.bucket << ' ' << c.variant << ' '
      << c.seed << ' ' << c.refitted << ' ' << c.collapsedMask << ' ' << c.storedBound << '\n';
    WriteV(o, Get(c.truth.origin));
    WriteV(o, Get(c.truth.axis));
    WriteV(o, Get(c.truth.xDirection));
    o << c.truth.radius << ' ' << c.truth.majorRadius << ' ' << c.truth.halfAngle << '\n';
    for (unsigned i = 0; i < c.analyticDomain.size(); ++i)
        o << c.analyticDomain[i] << (i + 1 == c.analyticDomain.size() ? '\n' : ' ');
    WriteSurface(o, *c.source);
    WriteSurface(o, *c.input);
    Require(bool(o), "Dataset write failed");
}
sggk::BSplineSurfacePtr Normalize(const sggk::BSplineSurface &s)
{
    auto ku = s.KnotsU(), kv = s.KnotsV();
    auto norm = [](auto &k) {
        double a = k.front(), b = k.back();
        for (auto &t : k)
            t = (t - a) / (b - a);
        k.front() = 0;
        k.back() = 1;
    };
    norm(ku);
    norm(kv);
    if (s.Weights().empty())
        return std::make_shared<sggk::BSplineSurface>(s.DegreeU(), s.DegreeV(), s.ControlPoints(),
                                                      ku, kv, s.MultsU(), s.MultsV());
    return std::make_shared<sggk::BSplineSurface>(s.DegreeU(), s.DegreeV(), s.ControlPoints(),
                                                  s.Weights(), ku, kv, s.MultsU(), s.MultsV());
}
double Bernstein(unsigned n, unsigned i, double t)
{
    double c = 1;
    for (unsigned j = 1; j <= i; ++j)
        c *= double(n + 1 - j) / j;
    return c * std::pow(t, int(i)) * std::pow(1 - t, int(n - i));
}
sggk::BSplineSurfacePtr Fit(const sggk::BSplineSurface &s, unsigned p, unsigned q)
{
    std::vector<std::vector<V>> a(p + 1, std::vector<V>(q + 1));
    V origin = Get(s.ControlPoints()[0][0]);
    for (unsigned i = 0; i <= p; ++i)
        for (unsigned j = 0; j <= q; ++j)
            a[i][j] = Get(s.CalcPoint(double(i) / p, double(j) / q)) - origin;
    auto solve = [](const std::vector<V> &v) {
        unsigned n = static_cast<unsigned>(v.size() - 1);
        std::vector<std::vector<double>> m(n + 1, std::vector<double>(n + 1));
        for (unsigned i = 0; i <= n; ++i)
            for (unsigned j = 0; j <= n; ++j)
                m[i][j] = Bernstein(n, j, double(i) / n);
        std::vector<V> out(n + 1);
        for (unsigned d = 0; d < 3; ++d)
        {
            std::vector<double> b(n + 1), x;
            for (unsigned i = 0; i <= n; ++i)
                b[i] = d == 0 ? v[i].x : (d == 1 ? v[i].y : v[i].z);
            Require(SolveLeastSquares(m, b, x), "Polynomial collocation is rank deficient");
            for (unsigned i = 0; i <= n; ++i)
            {
                if (d == 0)
                    out[i].x = x[i];
                else if (d == 1)
                    out[i].y = x[i];
                else
                    out[i].z = x[i];
            }
        }
        return out;
    };
    for (unsigned j = 0; j <= q; ++j)
    {
        std::vector<V> v(p + 1);
        for (unsigned i = 0; i <= p; ++i)
            v[i] = a[i][j];
        auto c = solve(v);
        for (unsigned i = 0; i <= p; ++i)
            a[i][j] = c[i];
    }
    sggk::Point3DMatrix cp(p + 1, std::vector<sggk::Point3D>(q + 1));
    for (unsigned i = 0; i <= p; ++i)
    {
        auto c = solve(a[i]);
        for (unsigned j = 0; j <= q; ++j)
            cp[i][j] = Point(c[j] + origin);
    }
    return std::make_shared<sggk::BSplineSurface>(
        p, q, cp, sggk::RealArray{0, 1}, sggk::RealArray{0, 1}, sggk::UIntArray{p + 1, p + 1},
        sggk::UIntArray{q + 1, q + 1});
}
void Snap(sggk::BSplineSurface &s, unsigned mask, const sggk::BSplineSurface &reference)
{
    unsigned nu = static_cast<unsigned>(s.ControlPoints().size()),
             nv = static_cast<unsigned>(s.ControlPoints()[0].size());
    for (unsigned k = 0; k < 4; ++k)
        if (mask & (1u << k))
        {
            auto point = reference.CalcPoint(k == 1 ? 1 : 0, k == 3 ? 1 : 0);
            for (unsigned j = 0; j < (k < 2 ? nv : nu); ++j)
                s.SetControlPoint(point, k < 2 ? (k == 0 ? 0 : nu - 1) : j,
                                  k < 2 ? j : (k == 2 ? 0 : nv - 1));
        }
}
void Pad(sggk::BSplineSurface &s, unsigned nu, unsigned nv, unsigned variant)
{
    for (unsigned axis = 0; axis < 2; ++axis)
    {
        auto size = [&] { return axis ? s.ControlPoints()[0].size() : s.ControlPoints().size(); };
        unsigned target = axis ? nv : nu, k = 1;
        Require(size() <= target, "Fitting exceeds control budget");
        while (size() < target)
        {
            double t = std::fmod(double(k++) * 0.6180339887498949, 1.0);
            if (variant % 2)
                t = std::pow(t, 1.4);
            unsigned p = axis ? s.DegreeV() : s.DegreeU();
            unsigned mult = variant % 3 == 0 && target - size() >= 2 ? std::min(2u, p) : 1;
            if (axis)
                s.InsertVKnots(t, mult);
            else
                s.InsertUKnots(t, mult);
        }
    }
}
struct CollocationAxis
{
    sggk::RealArray knots, nodes;
    sggk::UIntArray mults;
    std::vector<std::vector<double>> matrix;
    CollocationAxis(unsigned p, unsigned n)
    {
        Require((p == 3 || p == 5) && n > p, "Invalid low-degree fit dimensions");
        std::vector<double> full;
        for (unsigned i = 0; i <= n - p; ++i)
        {
            double t = double(i) / (n - p);
            unsigned mult = i == 0 || i == n - p ? p + 1 : 1;
            knots.push_back(t);
            mults.push_back(mult);
            full.insert(full.end(), mult, t);
        }
        // Greville interpolation on a clamped spline space. Extra accuracy comes
        // from simple internal knots, never from increasing the degree.
        for (unsigned i = 0; i < n; ++i)
        {
            double t = 0;
            for (unsigned j = 1; j <= p; ++j)
                t += full[i + j];
            t /= p;
            nodes.push_back(t);
            unsigned span =
                t == 1 ? n - 1
                       : static_cast<unsigned>(std::upper_bound(full.begin(), full.end(), t) -
                                               full.begin() - 1);
            std::vector<double> b(p + 1), left(p + 1), right(p + 1), row(n);
            b[0] = 1;
            for (unsigned j = 1; j <= p; ++j)
            {
                left[j] = t - full[span + 1 - j];
                right[j] = full[span + j] - t;
                double saved = 0;
                for (unsigned k = 0; k < j; ++k)
                {
                    double temp = b[k] / (right[k + 1] + left[j - k]);
                    b[k] = saved + right[k + 1] * temp;
                    saved = left[j - k] * temp;
                }
                b[j] = saved;
            }
            for (unsigned j = 0; j <= p; ++j)
                row[span - p + j] = b[j];
            matrix.push_back(std::move(row));
        }
    }
    std::vector<V> Solve(const std::vector<V> &values) const
    {
        const size_t n = nodes.size();
        std::vector<V> out(n);
        for (unsigned d = 0; d < 3; ++d)
        {
            std::vector<double> b(n), x;
            for (size_t i = 0; i < n; ++i)
                b[i] = d == 0 ? values[i].x : (d == 1 ? values[i].y : values[i].z);
            Require(SolveLeastSquares(matrix, b, x), "Spline collocation is rank deficient");
            for (size_t i = 0; i < n; ++i)
            {
                if (d == 0)
                    out[i].x = x[i];
                else if (d == 1)
                    out[i].y = x[i];
                else
                    out[i].z = x[i];
            }
        }
        return out;
    }
};
sggk::BSplineSurfacePtr FitSpline(const sggk::BSplineSurface &s, unsigned p, unsigned q,
                                  unsigned nu, unsigned nv)
{
    CollocationAxis u(p, nu), v(q, nv);
    V origin = Get(s.ControlPoints()[0][0]);
    std::vector<std::vector<V>> a(nu, std::vector<V>(nv));
    for (unsigned i = 0; i < nu; ++i)
        for (unsigned j = 0; j < nv; ++j)
            a[i][j] = Get(s.CalcPoint(u.nodes[i], v.nodes[j])) - origin;
    for (unsigned j = 0; j < nv; ++j)
    {
        std::vector<V> column(nu);
        for (unsigned i = 0; i < nu; ++i)
            column[i] = a[i][j];
        auto solved = u.Solve(column);
        for (unsigned i = 0; i < nu; ++i)
            a[i][j] = solved[i];
    }
    sggk::Point3DMatrix cp(nu, std::vector<sggk::Point3D>(nv));
    for (unsigned i = 0; i < nu; ++i)
    {
        auto solved = v.Solve(a[i]);
        for (unsigned j = 0; j < nv; ++j)
            cp[i][j] = Point(solved[j] + origin);
    }
    return std::make_shared<sggk::BSplineSurface>(p, q, cp, u.knots, v.knots, u.mults, v.mults);
}
sggk::BSplineSurfacePtr TryFitLowDegree(const DatasetCase &c, unsigned targetU, unsigned targetV)
{
    // The stored source may already have exchanged UV. Detect linear directions
    // from that source instead of assuming that V is the cylinder/cone axis.
    unsigned p = c.source->DegreeU() == 1 ? 3 : 5;
    unsigned q = c.source->DegreeV() == 1 ? 3 : 5;
    const unsigned maxU = p == 3 ? p + 1 : targetU, maxV = q == 3 ? q + 1 : targetV;
    Require(p < targetU && q < targetV && targetU * targetV < 200,
            "Low-degree fit exceeds control budget");
    for (unsigned extra = 0; extra <= maxU + maxV - p - q - 2; ++extra)
        for (unsigned nu = p + 1; nu <= maxU; ++nu)
        {
            const int nvValue = int(p + q + 2 + extra) - int(nu);
            if (nvValue < int(q + 1) || nvValue > int(maxV))
                continue;
            auto poly = FitSpline(*c.source, p, q, nu, static_cast<unsigned>(nvValue));
            Snap(*poly, c.collapsedMask, *c.source);
            // Sampling only rejects poor trials. Acceptance always requires the
            // independent continuous interval certificate after padding/snapping.
            double sampled = 0;
            for (unsigned i = 0; i <= 24; ++i)
                for (unsigned j = 0; j <= 24; ++j)
                    sampled = std::max(
                        sampled, Norm(Get(c.source->CalcPoint(double(i) / 24, double(j) / 24)) -
                                      Get(poly->CalcPoint(double(i) / 24, double(j) / 24))));
            if (!(sampled < 8e-7))
                continue;
            Pad(*poly, targetU, targetV, c.variant);
            Snap(*poly, c.collapsedMask, *c.source);
            if (BoundPair(*c.source, *poly, 8e-7).passed)
                return poly;
        }
    return {};
}
sggk::BSplineSurfacePtr FitLowDegree(const DatasetCase &c, unsigned targetU, unsigned targetV)
{
    if (auto fit = TryFitLowDegree(c, targetU, targetV))
        return fit;
    // Keep the original grid whenever possible. If it is insufficient, allocate
    // more controls within the SAME workload band, preferring smaller products.
    const unsigned hi[] = {79, 109, 139, 169, 199};
    unsigned minU = c.source->DegreeU() == 1 ? 4 : 6;
    unsigned minV = c.source->DegreeV() == 1 ? 4 : 6;
    std::vector<std::pair<unsigned, unsigned>> grids;
    for (unsigned nu = minU; nu * minV <= hi[c.bucket]; ++nu)
        for (unsigned nv = minV; nu * nv <= hi[c.bucket]; ++nv)
            if (nu * nv >= targetU * targetV && (nu != targetU || nv != targetV))
                grids.emplace_back(nu, nv);
    std::sort(grids.begin(), grids.end(), [](auto a, auto b) {
        if (a.first * a.second != b.first * b.second)
            return a.first * a.second < b.first * b.second;
        int da = std::abs(int(a.first) - int(a.second)),
            db = std::abs(int(b.first) - int(b.second));
        return da == db ? a < b : da < db;
    });
    for (const auto &[nu, nv] : grids)
        if (auto fit = TryFitLowDegree(c, nu, nv))
            return fit;
    throw std::runtime_error("No certified degree-3/5 spline within control band for pair " +
                             std::to_string(c.pair));
}
void ManifestHeader(std::ostream &o)
{
    o << "id,pair,type,source_kind,bucket,nu,nv,degree_u,degree_v,collapsed_mask,distance_upper\n"
      << std::setprecision(17);
}
void ManifestRow(std::ostream &o, const DatasetCase &c)
{
    o << c.id << ',' << c.pair << ',' << TypeName(c.truth.type) << ','
      << (c.refitted ? "refitted" : "analytic") << ',' << c.bucket << ','
      << c.input->ControlPoints().size() << ',' << c.input->ControlPoints()[0].size() << ','
      << c.input->DegreeU() << ',' << c.input->DegreeV() << ',' << c.collapsedMask << ','
      << c.storedBound << '\n';
}
void CheckIdentity(const DatasetCase &c, unsigned id)
{
    Require(id < 500 && c.id == id && c.pair == id / 2 && c.type == id / 100 &&
                c.bucket == (id % 100) / 20 && c.variant == (id % 20) / 2 &&
                c.refitted == bool(id % 2) && c.seed == 0x500b2a17u + c.pair,
            "Dataset500 metadata mismatch: " + std::to_string(id));
}
void Orient(sggk::BSplineSurface &s, unsigned variant)
{
    if (variant & 1)
        s.UReverse();
    if (variant & 2)
        s.VReverse();
    if (variant & 4)
        s.ExchangeUV();
}
DatasetCase Reference(unsigned type, unsigned bucket, unsigned variant)
{
    DatasetCase c;
    c.type = type;
    c.bucket = bucket;
    c.variant = variant;
    c.pair = type * 50 + bucket * 10 + variant;
    c.seed = 0x500b2a17u + c.pair;
    auto f = [&](unsigned k) { return Number(c.seed, k); };
    auto &a = c.truth;
    a.type = DatasetType(type);
    V z = Unit(V{f(0) - .5, f(1) - .5, .3 + f(2)}), x = Perpendicular(z);
    a.origin = Point({40 * (f(3) - .5), 40 * (f(4) - .5), 40 * (f(5) - .5)});
    a.axis = Direction(z);
    a.xDirection = Direction(x);
    a.radius = .7 + 5 * f(6);
    a.halfAngle = .25 + .65 * f(7);
    a.majorRadius = a.radius * (1.7 + 2 * f(8));
    sggk::Ucs3D frame(a.origin, a.xDirection, a.axis);
    sggk::SurfacePtr surface;
    switch (type)
    {
    case 0:
        surface = std::make_shared<sggk::Plane>(frame);
        break;
    case 1:
        surface = std::make_shared<sggk::CylinderSurface>(frame, a.radius);
        break;
    case 2:
        surface = std::make_shared<sggk::ConeSurface>(
            sggk::ConeSurface::CreateFromHalfAngle(frame, a.halfAngle));
        break;
    case 3:
        surface = std::make_shared<sggk::SphereSurface>(frame, a.radius);
        break;
    default:
        surface = std::make_shared<sggk::TorusSurface>(frame, a.majorRadius, a.radius);
    }
    double u0 = variant % 5 == 0 ? 6.1 : .1 + 4 * f(9), u1 = u0 + .36 + .25 * f(10), v0 = 0, v1 = 1;
    if (type == 0)
    {
        u0 = -2;
        u1 = u0 + 1 + 6 * f(10);
        v0 = -1;
        v1 = 1 + 5 * f(11);
    }
    if (type == 1)
    {
        v0 = -2;
        v1 = 1 + 4 * f(11);
    }
    if (type == 2)
    {
        v0 = variant >= 8 ? 0 : 1 + f(12);
        v1 = v0 + 1 + 4 * f(11);
    }
    if (type == 3)
    {
        v0 = .45 + 1.3 * f(12);
        v1 = v0 + .36 + .25 * f(11);
        if (variant == 8)
        {
            v0 = 0;
            v1 = .5;
        }
        if (variant == 9)
        {
            v0 = Pi - .5;
            v1 = Pi;
        }
    }
    if (type == 4)
    {
        v0 = variant % 5 == 0 ? 6.1 : .2 + 4 * f(12);
        v1 = v0 + .36 + .25 * f(11);
    }
    c.analyticDomain = {u0, u1, v0, v1};
    c.source = Normalize(
        *surface->ToBSpline(sggk::UVRange(sggk::Interval(u0, u1), sggk::Interval(v0, v1))));
    // Snap only mathematically collapsed analytic boundaries (roundoff at sin(pi)).
    if ((type == 2 || type == 3) && variant >= 8)
    {
        unsigned j = type == 3 && variant == 9
                         ? static_cast<unsigned>(c.source->ControlPoints()[0].size() - 1)
                         : 0;
        V point = Get(a.origin);
        if (type == 3)
            point = point + z * (variant == 9 ? a.radius : -a.radius);
        for (unsigned i = 0; i < c.source->ControlPoints().size(); ++i)
            c.source->SetControlPoint(Point(point), i, j);
    }
    c.collapsedMask = CollapsedMask(*c.source);
    return c;
}
} // namespace
void SaveDatasetCase(const DatasetCase &c, const std::filesystem::path &file)
{
    Save(c, file);
}
sggk::BSplineSurfacePtr NormalizeDatasetSurface(const sggk::BSplineSurface &s)
{
    return Normalize(s);
}
sggk::BSplineSurfacePtr FitDatasetPolynomial(const sggk::BSplineSurface &s, unsigned p, unsigned q)
{
    return Fit(s, p, q);
}
void PadDatasetSurface(sggk::BSplineSurface &s, unsigned nu, unsigned nv, unsigned variant)
{
    Pad(s, nu, nv, variant);
}
unsigned CollapsedMask(const sggk::BSplineSurface &s)
{
    const auto &cp = s.ControlPoints();
    unsigned nu = static_cast<unsigned>(cp.size()), nv = static_cast<unsigned>(cp[0].size()),
             mask = 0;
    for (unsigned k = 0; k < 4; ++k)
    {
        V origin = Get(cp[k == 1 ? nu - 1 : 0][k == 3 ? nv - 1 : 0]);
        double error = 0;
        for (unsigned j = 0; j < (k < 2 ? nv : nu); ++j)
            error = std::max(
                error,
                Norm(Get(cp[k < 2 ? (k == 0 ? 0 : nu - 1) : j][k < 2 ? j : (k == 2 ? 0 : nv - 1)]) -
                     origin));
        if (error < 1e-11)
            mask |= 1u << k;
    }
    return mask;
}
sggk::BSplineSurfacePtr RebuildDatasetSource(unsigned type, unsigned bucket, unsigned variant)
{
    auto c = Reference(type, bucket, variant);
    Orient(*c.source, variant);
    return Normalize(*c.source);
}
DatasetCase ReadDatasetCase(const std::filesystem::path &file)
{
    std::ifstream in(file);
    std::string magic;
    in >> magic;
    Require(magic == "B2A_CASE_V1", "Invalid dataset header");
    DatasetCase c;
    in >> c.id >> c.pair >> c.type >> c.bucket >> c.variant >> c.seed >> c.refitted >>
        c.collapsedMask >> c.storedBound;
    c.truth.type = DatasetType(c.type);
    c.truth.origin = Point(ReadV(in));
    c.truth.axis = Direction(ReadV(in));
    c.truth.xDirection = Direction(ReadV(in));
    in >> c.truth.radius >> c.truth.majorRadius >> c.truth.halfAngle;
    for (auto &x : c.analyticDomain)
        in >> x;
    c.source = ReadSurface(in);
    c.input = ReadSurface(in);
    in >> std::ws;
    Require(in.eof(), "Trailing dataset data");
    return c;
}
void CheckDataset500Contract(const DatasetCase &c, unsigned id)
{
    CheckIdentity(c, id);
    Require(c.source && c.input, "Missing Dataset500 surface");
    Require(c.source->DegreeU() <= 5 && c.source->DegreeV() <= 5 && c.input->DegreeU() <= 5 &&
                c.input->DegreeV() <= 5,
            "Dataset500 degree exceeds 5: " + std::to_string(id));
    const auto &net = c.input->ControlPoints();
    const size_t nu = net.size(), nv = net[0].size(), controls = nu * nv;
    const unsigned lo[] = {50, 80, 110, 140, 170}, hi[] = {79, 109, 139, 169, 199};
    Require(controls >= lo[c.bucket] && controls <= hi[c.bucket],
            "Dataset500 control count outside allocated band: " + std::to_string(id));
    Require(c.input->IsRational() == (!c.refitted && c.type != 0),
            "Dataset500 rationality mismatch");
    Require(std::isfinite(c.storedBound) && c.storedBound >= 0 &&
                c.storedBound < (c.refitted ? 1e-6 : 1e-9),
            "Invalid stored fixture bound");
    bool expectedCollapsed = (c.type == 2 || c.type == 3) && c.variant >= 8;
    Require(bool(c.collapsedMask) == expectedCollapsed &&
                CollapsedMask(*c.source) == c.collapsedMask &&
                CollapsedMask(*c.input) == c.collapsedMask,
            "Collapsed boundary mismatch");
    for (unsigned edge = 0; edge < 4; ++edge)
        if (c.collapsedMask & (1u << edge))
        {
            auto origin = Get(net[edge == 1 ? nu - 1 : 0][edge == 3 ? nv - 1 : 0]);
            for (size_t j = 0; j < (edge < 2 ? nv : nu); ++j)
            {
                auto point = Get(net[edge < 2 ? (edge == 0 ? 0 : nu - 1) : j]
                                    [edge < 2 ? j : (edge == 2 ? 0 : nv - 1)]);
                Require(point.x == origin.x && point.y == origin.y && point.z == origin.z,
                        "Degenerate control row is not exactly collapsed");
            }
        }
}
void AuditDataset500(const std::filesystem::path &dir)
{
    unsigned files = 0, rational = 0, collapsed = 0;
    for (const auto &f : std::filesystem::directory_iterator(dir))
        if (f.path().extension() == ".b2a" && f.path().filename().string().rfind("case-", 0) == 0)
            ++files;
    Require(files == 500, "Dataset must contain exactly 500 cases");
    for (unsigned id = 0; id < 500; ++id)
    {
        auto c = ReadDatasetCase(dir / CaseFilename(id));
        CheckDataset500Contract(c, id);
        auto source = RebuildDatasetSource(c.type, c.bucket, c.variant);
        Require(BoundPair(*source, *c.source, 1e-9).passed, "Stored source differs from generator");
        Require(BoundPair(*c.source, *c.input, c.refitted ? 1e-6 : 1e-9).passed,
                "Reloaded fixture certificate failed: " + std::to_string(id));
        rational += c.input->IsRational();
        collapsed += bool(c.collapsedMask);
        if (id % 100 == 99)
            std::cout << "Audited " << id + 1 << "/500" << std::endl;
    }
    Require(rational == 200 && collapsed == 40, "Dataset500 quota mismatch");
}
void RebuildLowDegreeDataset500(const std::filesystem::path &inputDir,
                                const std::filesystem::path &outputDir)
{
    Require(!std::filesystem::exists(outputDir) || std::filesystem::is_empty(outputDir),
            "Migration output directory must be empty; original data is never overwritten");
    std::filesystem::create_directories(outputDir);
    std::ofstream manifest(outputDir / "manifest.csv"), changes(outputDir / "migration.csv");
    Require(bool(manifest) && bool(changes), "Cannot write migration evidence");
    ManifestHeader(manifest);
    changes << std::setprecision(17)
            << "id,changed,old_degree_u,old_degree_v,new_degree_u,new_degree_v,old_nu,old_nv,"
               "new_nu,new_nv,old_spans_u,old_spans_v,new_spans_u,new_spans_v,certified_bound\n";
    unsigned changed = 0;
    for (unsigned id = 0; id < 500; ++id)
    {
        auto c = ReadDatasetCase(inputDir / CaseFilename(id));
        CheckIdentity(c, id);
        const auto old = c.input;
        unsigned nu = static_cast<unsigned>(old->ControlPoints().size());
        unsigned nv = static_cast<unsigned>(old->ControlPoints()[0].size());
        bool rebuild = old->DegreeU() > 5 || old->DegreeV() > 5 || nu * nv >= 200;
        if (rebuild)
        {
            if (nu * nv >= 200)
            {
                Require(nu * nv == 200 && c.bucket == 4, "Unexpected legacy control budget");
                nu = nv = 14;
            }
            if (c.refitted)
            {
                c.input = FitLowDegree(c, nu, nv);
                nu = static_cast<unsigned>(c.input->ControlPoints().size());
                nv = static_cast<unsigned>(c.input->ControlPoints()[0].size());
            }
            else
            {
                // Rebuild an exact representation from the already oriented source.
                c.input = std::make_shared<sggk::BSplineSurface>(*c.source);
                Require(old->DegreeU() <= 5 && old->DegreeV() <= 5 &&
                            old->DegreeU() >= c.source->DegreeU() &&
                            old->DegreeV() >= c.source->DegreeV(),
                        "Unexpected high-degree exact input");
                c.input->UDegreeElevation(old->DegreeU() - c.source->DegreeU());
                c.input->VDegreeElevation(old->DegreeV() - c.source->DegreeV());
                Pad(*c.input, nu, nv, c.variant);
                Snap(*c.input, c.collapsedMask, *c.source);
            }
            auto bound = BoundPair(*c.source, *c.input, c.refitted ? 1e-6 : 1e-9);
            Require(bound.passed, "Migrated fixture is not certified: " + std::to_string(id));
            c.storedBound = bound.upper;
            CheckDataset500Contract(c, id);
            Save(c, outputDir / CaseFilename(id));
            ++changed;
        }
        else
        {
            CheckDataset500Contract(c, id);
            std::filesystem::copy_file(inputDir / CaseFilename(id), outputDir / CaseFilename(id));
        }
        ManifestRow(manifest, c);
        changes << id << ',' << rebuild << ',' << old->DegreeU() << ',' << old->DegreeV() << ','
                << c.input->DegreeU() << ',' << c.input->DegreeV() << ','
                << old->ControlPoints().size() << ',' << old->ControlPoints()[0].size() << ',' << nu
                << ',' << nv << ',' << old->KnotsU().size() - 1 << ',' << old->KnotsV().size() - 1
                << ',' << c.input->KnotsU().size() - 1 << ',' << c.input->KnotsV().size() - 1 << ','
                << c.storedBound << '\n';
        if (id % 20 == 19)
            std::cout << "Migrated " << id + 1 << "/500; changed=" << changed << std::endl;
    }
    Require(bool(manifest) && bool(changes), "Migration evidence write failed");
    manifest.close();
    changes.close();
    // Read every file back and independently certify the complete batch before use.
    AuditDataset500(outputDir);
}
void GenerateDataset500(const std::filesystem::path &dir)
{
    std::filesystem::create_directories(dir);
    std::ofstream manifest(dir / "manifest.csv");
    ManifestHeader(manifest);
    for (unsigned type = 0; type < 5; ++type)
        for (unsigned bucket = 0; bucket < 5; ++bucket)
            for (unsigned variant = 0; variant < 10; ++variant)
            {
                auto base = Reference(type, bucket, variant);
                unsigned sizes[] = {8, 10, 11, 12, 14}, nu = sizes[bucket], nv = nu;
                if (type == 0 && bucket == 0 && variant == 0)
                {
                    nu = 5;
                    nv = 10;
                }
                for (unsigned mode = 0; mode < 2; ++mode)
                {
                    auto c = base;
                    c.id = 2 * c.pair + mode;
                    c.refitted = mode != 0;
                    if (!mode)
                    {
                        c.input = std::make_shared<sggk::BSplineSurface>(*base.source);
                        unsigned elevate = variant % 3;
                        if (elevate)
                        {
                            c.input->UDegreeElevation(elevate);
                            c.input->VDegreeElevation(elevate);
                        }
                        Pad(*c.input, nu, nv, variant);
                        Snap(*c.input, c.collapsedMask, *base.source);
                    }
                    else
                        c.input = FitLowDegree(c, nu, nv);
                    c.source = std::make_shared<sggk::BSplineSurface>(*base.source);
                    Orient(*c.source, variant);
                    Orient(*c.input, variant);
                    // Rebuild the clamped [0,1] representation to remove SDK endpoint ulps.
                    c.source = Normalize(*c.source);
                    c.input = Normalize(*c.input);
                    c.collapsedMask = CollapsedMask(*c.source);
                    Snap(*c.input, c.collapsedMask, *c.source);
                    auto bound = BoundPair(*c.source, *c.input, mode ? 1e-6 : 1e-9);
                    Require(bound.passed, "Final pair bound failed: " + std::to_string(c.id));
                    c.storedBound = bound.upper;
                    CheckDataset500Contract(c, c.id);
                    Save(c, dir / CaseFilename(c.id));
                    ManifestRow(manifest, c);
                }
                if (variant == 9)
                    std::cout << "Generated " << (2 * base.pair + 2) << "/500" << std::endl;
            }
    Require(bool(manifest), "Manifest write failed");
    manifest.close();
    AuditDataset500(dir);
}
} // namespace bspline2analytic::tests
