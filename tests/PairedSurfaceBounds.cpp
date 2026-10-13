#include "PairedSurfaceBounds.h"
#include "IntervalBounds.h"
#include <functional>
namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
struct Cell
{
    BV x;
    BPoly w;
    double u0, u1, v0, v1;
};
std::vector<double> Expand(const sggk::RealArray &k, const sggk::UIntArray &m)
{
    std::vector<double> r;
    for (size_t i = 0; i < k.size(); ++i)
        r.insert(r.end(), m[i], k[i]);
    return r;
}
// Interval blossom extraction uses the input control net directly, avoiding an
// unaccounted floating-point SDK knot-insertion error in the certificate.
std::vector<std::vector<I>> Extract(const std::vector<double> &k, unsigned p, unsigned span)
{
    std::vector<std::vector<I>> e(p + 1, std::vector<I>(p + 1));
    for (unsigned j = 0; j <= p; ++j)
    {
        std::vector<std::vector<I>> d(p + 1, std::vector<I>(p + 1));
        for (unsigned i = 0; i <= p; ++i)
            d[i][i] = I(1);
        for (unsigned r = 1; r <= p; ++r)
        {
            double t = r <= p - j ? k[span] : k[span + 1];
            for (unsigned i = p; i >= r; --i)
            {
                I alpha =
                    (I(t) - I(k[span - p + i])) / (I(k[span + 1 + i - r]) - I(k[span - p + i]));
                for (unsigned l = 0; l <= p; ++l)
                    d[i][l] = d[i - 1][l] * (I(1) - alpha) + d[i][l] * alpha;
            }
        }
        e[j] = d[p];
    }
    return e;
}
std::vector<Cell> Cells(const sggk::BSplineSurface &s, V origin)
{
    auto ku = Expand(s.KnotsU(), s.MultsU()), kv = Expand(s.KnotsV(), s.MultsV());
    unsigned p = s.DegreeU(), q = s.DegreeV();
    const auto &cp = s.ControlPoints();
    const auto &w = s.Weights();
    std::vector<Cell> cells;
    for (unsigned a = p; a < cp.size(); ++a)
        if (ku[a] < ku[a + 1])
        {
            auto eu = Extract(ku, p, a);
            for (unsigned b = q; b < cp.front().size(); ++b)
                if (kv[b] < kv[b + 1])
                {
                    auto ev = Extract(kv, q, b);
                    Cell c{{BPoly(p, q), BPoly(p, q), BPoly(p, q)},
                           BPoly(p, q),
                           ku[a],
                           ku[a + 1],
                           kv[b],
                           kv[b + 1]};
                    for (unsigned i = 0; i <= p; ++i)
                        for (unsigned j = 0; j <= q; ++j)
                        {
                            IV x{I(0), I(0), I(0)};
                            I sum(0);
                            for (unsigned k = 0; k <= p; ++k)
                                for (unsigned l = 0; l <= q; ++l)
                                {
                                    unsigned u = a - p + k, v = b - q + l;
                                    I wt(w.empty() ? 1 : w[u][v]);
                                    if (!(wt.lo > 0))
                                        throw std::runtime_error(
                                            "Pair certificate requires positive weights");
                                    I f = eu[i][k] * ev[j][l] * wt;
                                    x = x + (Exact(Get(cp[u][v])) - Exact(origin)) * f;
                                    sum = sum + f;
                                }
                            c.x.x.at(i, j) = x.x;
                            c.x.y.at(i, j) = x.y;
                            c.x.z.at(i, j) = x.z;
                            c.w.at(i, j) = sum;
                        }
                    if (w.empty())
                        c.w = BPoly(0, 0, I(1));
                    cells.push_back(std::move(c));
                }
        }
    return cells;
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
        auto put = [&](BPoly &o, unsigned i, I x) {
            if (u)
                o.at(i, j) = x;
            else
                o.at(j, i) = x;
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
    return {a, b};
}
Cell Crop(Cell c, double u0, double u1, double v0, double v1)
{
    auto axis = [&](bool u, double lo, double hi) {
        double &a = u ? c.u0 : c.v0, &b = u ? c.u1 : c.v1;
        auto cut = [&](double t, bool right) {
            I f = (I(t) - I(a)) / (I(b) - I(a));
            auto apply = [&](BPoly &p) {
                auto q = Cut(p, u, f);
                p = right ? q.second : q.first;
            };
            apply(c.x.x);
            apply(c.x.y);
            apply(c.x.z);
            apply(c.w);
            if (right)
                a = t;
            else
                b = t;
        };
        if (lo > a)
            cut(lo, true);
        if (hi < b)
            cut(hi, false);
    };
    axis(true, u0, u1);
    axis(false, v0, v1);
    return c;
}
} // namespace
PairBound BoundPair(const sggk::BSplineSurface &a, const sggk::BSplineSurface &b, double tolerance,
                    unsigned maxPieces)
{
    PairBound out;
    if (!IntervalEnvironmentSupported() || !(tolerance > 0) || !std::isfinite(tolerance) ||
        !maxPieces)
        return out;
    if (a.MinParamU() != b.MinParamU() || a.MaxParamU() != b.MaxParamU() ||
        a.MinParamV() != b.MinParamV() || a.MaxParamV() != b.MaxParamV())
        return out;
    V origin = Get(a.ControlPoints().front().front());
    auto aa = Cells(a, origin), bb = Cells(b, origin);
    if (aa.empty() || bb.empty())
        return out;
    std::function<double(const BV &, const BPoly &, unsigned)> visit;
    visit = [&](const BV &e, const BPoly &w, unsigned depth) -> double {
        if (++out.pieces > maxPieces)
            return INFINITY;
        auto z = Range(e);
        auto wr = w.range();
        double bound = wr.lo > 0 ? (Sqrt(Sqr(z.x) + Sqr(z.y) + Sqr(z.z)) / I(wr.lo)).hi : INFINITY;
        if (bound < tolerance)
            return bound;
        if (depth >= 18)
            return INFINITY;
        bool u = depth % 2 == 0;
        auto x = Split(e.x, u), y = Split(e.y, u), v = Split(e.z, u), d = Split(w, u);
        double left = visit({x.first, y.first, v.first}, d.first, depth + 1);
        if (!std::isfinite(left))
            return INFINITY;
        return std::max(left, visit({x.second, y.second, v.second}, d.second, depth + 1));
    };
    double bound = 0;
    for (const auto &x : aa)
        for (const auto &y : bb)
        {
            double u0 = std::max(x.u0, y.u0), u1 = std::min(x.u1, y.u1), v0 = std::max(x.v0, y.v0),
                   v1 = std::min(x.v1, y.v1);
            if (u0 >= u1 || v0 >= v1)
                continue;
            auto p = Crop(x, u0, u1, v0, v1), q = Crop(y, u0, u1, v0, v1);
            auto e = p.x * q.w - q.x * p.w;
            auto w = p.w * q.w;
            bound = std::max(bound, visit(e, w, 0));
            if (!std::isfinite(bound))
                return out;
        }
    out.upper = bound;
    out.passed = out.pieces > 0 && bound < tolerance;
    return out;
}
} // namespace bspline2analytic::tests
