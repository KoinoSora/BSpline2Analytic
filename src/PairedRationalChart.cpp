#include "PairedValidationInternal.h"

namespace bspline2analytic::detail::paired
{
namespace
{
BPoly C(I x)
{
    return BPoly(0, 0, x);
}
I Choose(unsigned n, unsigned k)
{
    I r(1);
    for (unsigned i = 1; i <= k; ++i)
        r = r * I(n - k + i) / I(i);
    return r;
}
std::pair<BPoly, BPoly> Cut(const BPoly &p, bool u, I t)
{
    BPoly a(p.u, p.v), b(p.u, p.v);
    unsigned n = u ? p.u : p.v, m = u ? p.v : p.u;
    for (unsigned j = 0; j <= m; ++j)
    {
        std::vector<I> d(n + 1);
        for (unsigned i = 0; i <= n; ++i)
            d[i] = u ? p.at(i, j) : p.at(j, i);
        auto put = [&](BPoly &q, unsigned i, I x) {
            if (u)
                q.at(i, j) = x;
            else
                q.at(j, i) = x;
        };
        put(a, 0, d[0]);
        put(b, n, d[n]);
        for (unsigned r = 1; r <= n; ++r)
        {
            for (unsigned i = 0; i <= n - r; ++i)
                d[i] = d[i] * (I(1) - t) + d[i + 1] * t;
            put(a, r, d[0]);
            put(b, n - r, d[n - r]);
        }
    }
    return {std::move(a), std::move(b)};
}
BPoly Orient(const BPoly &p, const ParameterCorrespondence &m)
{
    BPoly q(m.swappedUV ? p.v : p.u, m.swappedUV ? p.u : p.v);
    for (unsigned i = 0; i <= q.u; ++i)
        for (unsigned j = 0; j <= q.v; ++j)
        {
            unsigned a = m.reversedU ? q.u - i : i, b = m.reversedV ? q.v - j : j;
            q.at(i, j) = m.swappedUV ? p.at(b, a) : p.at(a, b);
        }
    return q;
}
BV Orient(const BV &p, const ParameterCorrespondence &m)
{
    return {Orient(p.x, m), Orient(p.y, m), Orient(p.z, m)};
}
struct Axis
{
    BPoly c, s, w, t;
};
Axis MakeAxis(double lo, double hi, bool angular, bool u)
{
    if (!angular)
    {
        BPoly t(u ? 1 : 0, u ? 0 : 1);
        t.c = {I(lo), I(hi)};
        return {C(I(0)), C(I(0)), C(I(1)), t};
    }
    I a(lo), b(hi), mid = (a + b) / I(2), half = (b - a) / I(2);
    BPoly c(u ? 2 : 0, u ? 0 : 2), s = c, w = c;
    c.c = {Cos(a), Cos(mid), Cos(b)};
    s.c = {Sin(a), Sin(mid), Sin(b)};
    w.c = {I(1), Cos(half), I(1)};
    if (w.range().lo <= 0)
        throw Failure(Status::Unsupported, "Nonpositive rational arc weight");
    return {c, s, w, C(I(0))};
}
} // namespace
const std::vector<I> &Products::Weights(unsigned p, unsigned q)
{
    auto key = std::make_pair(p, q);
    auto it = factors.find(key);
    if (it != factors.end())
        return it->second;
    std::vector<I> w((p + 1) * (q + 1));
    for (unsigned i = 0; i <= p; ++i)
        for (unsigned j = 0; j <= q; ++j)
            w[i * (q + 1) + j] = Choose(p, i) * Choose(q, j) / Choose(p + q, i + j);
    return factors.emplace(key, std::move(w)).first->second;
}
BPoly Products::Mul(const BPoly &a, const BPoly &b)
{
    if (!a.u && !a.v)
        return b * a.c[0];
    if (!b.u && !b.v)
        return a * b.c[0];
    const auto &wu = Weights(a.u, b.u), &wv = Weights(a.v, b.v);
    BPoly r(a.u + b.u, a.v + b.v);
    for (unsigned i = 0; i <= a.u; ++i)
        for (unsigned j = 0; j <= a.v; ++j)
            for (unsigned k = 0; k <= b.u; ++k)
                for (unsigned l = 0; l <= b.v; ++l)
                    r.at(i + k, j + l) = r.at(i + k, j + l) + a.at(i, j) * b.at(k, l) *
                                                                  wu[i * (b.u + 1) + k] *
                                                                  wv[j * (b.v + 1) + l];
    return r;
}
BV Products::Mul(const BV &a, const BPoly &b)
{
    return {Mul(a.x, b), Mul(a.y, b), Mul(a.z, b)};
}
BPoly Products::Dot(const BV &a, const BV &b)
{
    return Mul(a.x, b.x) + Mul(a.y, b.y) + Mul(a.z, b.z);
}
BV Products::Cross(const BV &a, const BV &b)
{
    return {Mul(a.y, b.z) - Mul(a.z, b.y), Mul(a.z, b.x) - Mul(a.x, b.z),
            Mul(a.x, b.y) - Mul(a.y, b.x)};
}
BPoly Restrict(const BPoly &p, bool u, I a, I b)
{
    BPoly q = (b.lo == 1 && b.hi == 1) ? p : Cut(p, u, b).first;
    if (a.lo == 0 && a.hi == 0)
        return q;
    if (b.lo <= 0)
        throw Failure(Status::VerificationFailed, "Unresolved restriction divisor");
    return Cut(q, u, a / b).second;
}
BV Restrict(const BV &p, bool u, I a, I b)
{
    return {Restrict(p.x, u, a, b), Restrict(p.y, u, a, b), Restrict(p.z, u, a, b)};
}
BezierCell Crop(const BezierCell &cell, double u0, double u1, double v0, double v1)
{
    BezierCell c = cell;
    auto axis = [&](bool u, double lo, double hi, double a, double b) {
        I x = lo == a ? I(0) : (I(lo) - I(a)) / (I(b) - I(a));
        I y = hi == b ? I(1) : (I(hi) - I(a)) / (I(b) - I(a));
        c.p = Restrict(c.p, u, x, y);
        c.w = Restrict(c.w, u, x, y);
    };
    axis(true, u0, u1, cell.u0, cell.u1);
    axis(false, v0, v1, cell.v0, cell.v1);
    c.u0 = u0;
    c.u1 = u1;
    c.v0 = v0;
    c.v1 = v1;
    return c;
}
BV LocalSource(const BV &p, const AnalyticParameters &a)
{
    IV z = Normalize(Exact(Get(a.axis))), x = Exact(Get(a.xDirection));
    x = Normalize(x - z * detail::Dot(x, z));
    IV y = detail::Cross(z, x), o = Exact(Get(a.origin));
    BV d = {p.x - C(o.x), p.y - C(o.y), p.z - C(o.z)};
    return {detail::Dot(d, x), detail::Dot(d, y), detail::Dot(d, z)};
}
void CheckPatch(const AnalyticPatch &p)
{
    using T = sggk::SurfaceType;
    const auto &a = p.parameters;
    if (a.type != T::Plane && a.type != T::Cylinder && a.type != T::Cone && a.type != T::Sphere &&
        a.type != T::Torus)
        throw Failure(Status::Unsupported, "Unsupported mother surface");
    for (double x : p.domain)
        if (!std::isfinite(x))
            throw Failure(Status::InvalidInput, "Nonfinite domain");
    if (!(p.domain[0] < p.domain[1] && p.domain[2] < p.domain[3]) ||
        !(p.normalSign == 1 || p.normalSign == -1))
        throw Failure(Status::InvalidInput, "Invalid domain or normal orientation");
    auto f = FrameOf(a);
    if (!Finite(f.o) || !Finite(f.x) || !Finite(f.z))
        throw Failure(Status::InvalidInput, "Invalid frame");
    if (a.type != T::Plane && a.type != T::Cone && !(a.radius > 0 && std::isfinite(a.radius)))
        throw Failure(Status::InvalidInput, "Invalid radius");
    if (a.type != T::Plane && p.domain[1] - p.domain[0] > 2 * Pi + 1e-12)
        throw Failure(Status::Unsupported, "Multiple-turn chart");
    if (a.type == T::Sphere && !(p.domain[2] > 0 && p.domain[3] < Pi))
        throw Failure(Status::Unsupported, "Sphere pole chart");
    if (a.type == T::Cone && !(a.halfAngle > 0 && a.halfAngle < Pi / 2 && p.domain[2] > 0))
        throw Failure(Status::Unsupported, "Cone apex or invalid angle");
    if (a.type == T::Torus && (!(a.majorRadius > a.radius && std::isfinite(a.majorRadius)) ||
                               p.domain[3] - p.domain[2] > 2 * Pi + 1e-12))
        throw Failure(Status::Unsupported, "Only regular ring torus charts supported");
}
Chart BuildChart(const AnalyticPatch &p)
{
    using T = sggk::SurfaceType;
    Chart chart;
    const auto &a = p.parameters;
    chart.angular = {a.type != T::Plane, a.type == T::Sphere || a.type == T::Torus};
    for (unsigned d = 0; d < 2; ++d)
    {
        double lo = p.domain[2 * d], hi = p.domain[2 * d + 1];
        unsigned n = chart.angular[d] ? static_cast<unsigned>(std::ceil((hi - lo) / (Pi / 2))) : 1;
        for (unsigned i = 0; i <= n; ++i)
            chart.angles[d].push_back(i == n ? hi : lo + (hi - lo) * i / n);
    }
    const unsigned nu = static_cast<unsigned>(chart.angles[0].size() - 1),
                   nv = static_cast<unsigned>(chart.angles[1].size() - 1);
    Products products;
    for (unsigned i = 0; i < nu; ++i)
        for (unsigned j = 0; j < nv; ++j)
        {
            auto u = MakeAxis(chart.angles[0][i], chart.angles[0][i + 1], chart.angular[0], true);
            auto v = MakeAxis(chart.angles[1][j], chart.angles[1][j + 1], chart.angular[1], false);
            BV q, n;
            BPoly w = products.Mul(u.w, v.w), zero = C(I(0));
            switch (a.type)
            {
            case T::Plane:
                q = {u.t, v.t, zero};
                n = {zero, zero, C(I(1))};
                break;
            case T::Cylinder:
                q = {u.c * I(a.radius), u.s * I(a.radius), products.Mul(v.t, u.w)};
                n = {u.c, u.s, zero};
                break;
            case T::Cone: {
                I sn = Sin(I(a.halfAngle)), cs = Cos(I(a.halfAngle));
                q = {products.Mul(u.c, v.t) * (sn / cs), products.Mul(u.s, v.t) * (sn / cs),
                     products.Mul(v.t, u.w)};
                n = {u.c * cs, u.s * cs, u.w * (-sn)};
                break;
            }
            case T::Sphere:
                n = {products.Mul(u.c, v.s), products.Mul(u.s, v.s),
                     products.Mul(u.w, v.c) * I(-1)};
                q = n * I(a.radius);
                break;
            default: {
                auto radial = v.w * I(a.majorRadius) + v.c * I(a.radius);
                q = {products.Mul(u.c, radial), products.Mul(u.s, radial),
                     products.Mul(u.w, v.s) * I(a.radius)};
                n = {products.Mul(u.c, v.c), products.Mul(u.s, v.c), products.Mul(u.w, v.s)};
                break;
            }
            }
            double low[2] = {double(i) / nu, double(j) / nv},
                   high[2] = {double(i + 1) / nu, double(j + 1) / nv};
            unsigned du = p.mapping.swappedUV ? 1 : 0, dv = 1 - du;
            TargetCell cell;
            cell.surface = {Orient(q, p.mapping),
                            p.mapping.reversedU ? 1 - high[du] : low[du],
                            p.mapping.reversedU ? 1 - low[du] : high[du],
                            p.mapping.reversedV ? 1 - high[dv] : low[dv],
                            p.mapping.reversedV ? 1 - low[dv] : high[dv],
                            Orient(w, p.mapping)};
            cell.normal = Orient(n, p.mapping) * I(p.normalSign);
            chart.cells.push_back(std::move(cell));
        }
    return chart;
}
double RationalCoordinate(const Chart &c, unsigned d, double angle)
{
    const auto &a = c.angles[d];
    if (angle <= a.front())
        return 0;
    if (angle >= a.back())
        return 1;
    unsigned i = static_cast<unsigned>(std::upper_bound(a.begin(), a.end(), angle) - a.begin() - 1);
    double t = (angle - a[i]) / (a[i + 1] - a[i]);
    if (c.angular[d])
    {
        double h = (a[i + 1] - a[i]) * .5, k = std::tan((angle - a[i]) * .5);
        t = k / (std::sin(h) + k * (1 - std::cos(h)));
    }
    return (i + std::clamp(t, 0.0, 1.0)) / (a.size() - 1);
}
} // namespace bspline2analytic::detail::paired
