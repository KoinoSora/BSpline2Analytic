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
void GenerateDataset500(const std::filesystem::path &dir)
{
    std::filesystem::create_directories(dir);
    std::ofstream manifest(dir / "manifest.csv");
    manifest
        << "id,pair,type,source_kind,bucket,nu,nv,degree_u,degree_v,collapsed_mask,distance_upper\n"
        << std::setprecision(17);
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
                if (type == 4 && bucket == 4 && variant == 9)
                {
                    nu = 10;
                    nv = 20;
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
                    {
                        bool passed = false;
                        for (unsigned degree :
                             type == 0 ? std::vector<unsigned>{3} : std::vector<unsigned>{5, 7, 9})
                        {
                            unsigned p = degree, q = type >= 3 ? degree : 3;
                            if (p >= nu || q >= nv)
                                continue;
                            auto poly = Fit(*base.source, p, q);
                            Snap(*poly, c.collapsedMask, *base.source);
                            Pad(*poly, nu, nv, variant);
                            Snap(*poly, c.collapsedMask, *base.source);
                            auto check = BoundPair(*base.source, *poly, 8e-7);
                            if (check.passed)
                            {
                                c.input = poly;
                                passed = true;
                                break;
                            }
                        }
                        Require(passed, "No certified polynomial within control budget for pair " +
                                            std::to_string(c.pair));
                    }
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
                    const auto &cp = c.input->ControlPoints();
                    Require(cp.size() * cp[0].size() >= 50 && cp.size() * cp[0].size() <= 200,
                            "Control count outside [50,200]");
                    Save(c, dir / CaseFilename(c.id));
                    manifest << c.id << ',' << c.pair << ',' << TypeName(c.truth.type) << ','
                             << (mode ? "refitted" : "analytic") << ',' << bucket << ','
                             << cp.size() << ',' << cp[0].size() << ',' << c.input->DegreeU() << ','
                             << c.input->DegreeV() << ',' << c.collapsedMask << ',' << bound.upper
                             << '\n';
                }
                if (variant == 9)
                    std::cout << "Generated " << (2 * base.pair + 2) << "/500" << std::endl;
            }
    Require(bool(manifest), "Manifest write failed");
}
} // namespace bspline2analytic::tests
