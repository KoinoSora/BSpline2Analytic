#include "IntervalBounds.h"
#include <cfenv>
#include <xmmintrin.h>
namespace bspline2analytic::detail
{
bool IntervalEnvironmentSupported()
{
    return std::fegetround() == FE_TONEAREST && (_mm_getcsr() & 0x8040u) == 0;
}
static double Down(double x)
{
    return std::nextafter(x, -INFINITY);
}
static double Up(double x)
{
    return std::nextafter(x, INFINITY);
}
I IPi()
{
    return {Down(Pi), Up(Pi)};
}
I operator+(I a, I b)
{
    if (b.lo == 0 && b.hi == 0)
        return a;
    if (a.lo == 0 && a.hi == 0)
        return b;
    return {Down(a.lo + b.lo), Up(a.hi + b.hi)};
}
I operator-(I a, I b)
{
    return {Down(a.lo - b.hi), Up(a.hi - b.lo)};
}
I operator-(I a)
{
    return {-a.hi, -a.lo};
}
I operator*(I a, I b)
{
    if (a.lo == 1 && a.hi == 1)
        return b;
    if (b.lo == 1 && b.hi == 1)
        return a;
    if (a.finite() && b.finite() && ((a.lo == 0 && a.hi == 0) || (b.lo == 0 && b.hi == 0)))
        return I(0);
    double q[] = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    for (double x : q)
        if (std::isnan(x))
            return I::Whole();
    return {Down(*std::min_element(q, q + 4)), Up(*std::max_element(q, q + 4))};
}
I operator/(I a, I b)
{
    if (b.lo <= 0 && b.hi >= 0)
        return I::Whole();
    return a * I(Down(1 / b.hi), Up(1 / b.lo));
}
I Sqr(I a)
{
    if (a.lo <= 0 && a.hi >= 0)
        return {0, Up(a.absMax() * a.absMax())};
    return {Down(std::min(a.lo * a.lo, a.hi * a.hi)), Up(std::max(a.lo * a.lo, a.hi * a.hi))};
}
I Sqrt(I a)
{
    if (a.hi < 0)
        return I::Whole();
    return {a.lo <= 0 ? 0 : std::max(0.0, Down(std::sqrt(a.lo))),
            Up(std::sqrt(std::max(0.0, a.hi)))};
}
static I AtanPoint(double x)
{
    if (!std::isfinite(x))
        return x < 0 ? -IPi() / I(2) : IPi() / I(2);
    if (x < 0)
        return -AtanPoint(-x);
    if (x > 1)
        return IPi() / I(2) - Atan(I(1) / I(x));
    I z(x), offset(0);
    if (x > 0.4)
    {
        z = (z - I(1)) / (z + I(1));
        offset = IPi() / I(4);
    }
    I zz = z * z, term = z, sum = z;
    for (unsigned k = 1; k < 48; ++k)
    {
        term = term * zz;
        sum = sum + term * I(k % 2 ? -1 : 1) / I(2 * k + 1);
    }
    I tail(1);
    for (unsigned k = 0; k < 97; ++k)
        tail = tail * I(z.absMax());
    tail = tail / I(97);
    return offset + sum + I(-tail.hi, tail.hi);
}
I Atan(I a)
{
    auto l = AtanPoint(a.lo), h = AtanPoint(a.hi);
    return {l.lo, h.hi};
}
static I Trig(I x, bool cosine)
{
    if (!x.finite() || x.hi - x.lo > 2 * Pi)
        return {-1, 1};
    double middle = x.lo + (x.hi - x.lo) * 0.5;
    double k = std::round(middle / (2 * Pi));
    if (std::abs(k) > 1e9)
        return {-1, 1};
    I z = x - IPi() * I(2 * k), zz = z * z;
    I term = cosine ? I(1) : z, sum = term;
    for (unsigned n = 1; n <= 28; ++n)
    {
        double a = cosine ? 2 * n - 1 : 2 * n, b = cosine ? 2 * n : 2 * n + 1;
        term = -term * zz / I(a * b);
        sum = sum + term;
    }
    // Reduced |z| <= 2*pi: the order-57 remainder is below 1e-30.
    // Use a deliberately larger, rigorously safe remainder.
    if (z.absMax() > 7)
        return {-1, 1};
    I rem(1);
    for (unsigned j = 1; j <= 57; ++j)
        rem = rem * I(7) / I(j);
    return {std::max(-1.0, (sum - I(rem.hi)).lo), std::min(1.0, (sum + I(rem.hi)).hi)};
}
I Sin(I a)
{
    return Trig(a, false);
}
I Cos(I a)
{
    return Trig(a, true);
}
static I AcosPoint(double a)
{
    a = std::clamp(a, -1.0, 1.0);
    if (a == 1)
        return I(0);
    if (a == -1)
        return IPi();
    if (a == 0)
        return IPi() / I(2);
    I r = Sqrt(I(1) - Sqr(I(a))) / I(std::abs(a));
    return a > 0 ? Atan(r) : IPi() - Atan(r);
}
I Acos(I a)
{
    auto l = AcosPoint(a.hi), h = AcosPoint(a.lo);
    return {l.lo, h.hi};
}
IV Exact(V p)
{
    return {I(p.x), I(p.y), I(p.z)};
}
IV operator+(IV a, IV b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
IV operator-(IV a, IV b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
IV operator*(IV a, I b)
{
    return {a.x * b, a.y * b, a.z * b};
}
I Dot(IV a, IV b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
IV Cross(IV a, IV b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
IV Normalize(IV a)
{
    return a * (I(1) / Sqrt(Sqr(a.x) + Sqr(a.y) + Sqr(a.z)));
}
I BPoly::range() const
{
    if (c.empty())
        return I::Whole();
    I r = c.front();
    for (auto a : c)
    {
        r.lo = std::min(r.lo, a.lo);
        r.hi = std::max(r.hi, a.hi);
    }
    return r;
}
static BPoly Elevate(BPoly p, unsigned u, unsigned v)
{
    while (p.u < u || p.v < v)
    {
        bool along = p.u < u;
        unsigned n = along ? p.u : p.v, lines = along ? p.v : p.u;
        BPoly q(p.u + (along ? 1 : 0), p.v + (along ? 0 : 1));
        for (unsigned j = 0; j <= lines; ++j)
            for (unsigned i = 0; i <= n + 1; ++i)
            {
                I value(0);
                if (i > 0)
                    value = value + (along ? p.at(i - 1, j) : p.at(j, i - 1)) * (I(i) / I(n + 1));
                if (i <= n)
                    value = value + (along ? p.at(i, j) : p.at(j, i)) * (I(n + 1 - i) / I(n + 1));
                if (along)
                    q.at(i, j) = value;
                else
                    q.at(j, i) = value;
            }
        p = std::move(q);
    }
    return p;
}
BPoly operator+(const BPoly &x, const BPoly &y)
{
    unsigned u = std::max(x.u, y.u), v = std::max(x.v, y.v);
    BPoly a = Elevate(x, u, v), b = Elevate(y, u, v);
    for (size_t i = 0; i < a.c.size(); ++i)
        a.c[i] = a.c[i] + b.c[i];
    return a;
}
BPoly operator-(const BPoly &a, const BPoly &b)
{
    return a + b * I(-1);
}
BPoly operator*(const BPoly &a, I b)
{
    BPoly r = a;
    for (auto &x : r.c)
        x = x * b;
    return r;
}
static I Choose(unsigned n, unsigned k)
{
    I r(1);
    for (unsigned i = 1; i <= k; ++i)
        r = r * I(n - k + i) / I(i);
    return r;
}
BPoly operator*(const BPoly &a, const BPoly &b)
{
    BPoly r(a.u + b.u, a.v + b.v);
    std::vector<I> wu((a.u + 1) * (b.u + 1)), wv((a.v + 1) * (b.v + 1));
    for (unsigned i = 0; i <= a.u; ++i)
        for (unsigned k = 0; k <= b.u; ++k)
            wu[i * (b.u + 1) + k] = Choose(a.u, i) * Choose(b.u, k) / Choose(r.u, i + k);
    for (unsigned j = 0; j <= a.v; ++j)
        for (unsigned l = 0; l <= b.v; ++l)
            wv[j * (b.v + 1) + l] = Choose(a.v, j) * Choose(b.v, l) / Choose(r.v, j + l);
    for (unsigned i = 0; i <= a.u; ++i)
        for (unsigned k = 0; k <= b.u; ++k)
            for (unsigned j = 0; j <= a.v; ++j)
                for (unsigned l = 0; l <= b.v; ++l)
                    r.at(i + k, j + l) = r.at(i + k, j + l) + a.at(i, j) * b.at(k, l) *
                                                                  wu[i * (b.u + 1) + k] *
                                                                  wv[j * (b.v + 1) + l];
    return r;
}
BPoly Constant(const BPoly &p, I a)
{
    return BPoly(p.u, p.v, a);
}
BPoly Derivative(const BPoly &p, bool along)
{
    unsigned n = along ? p.u : p.v;
    if (!n)
        return BPoly(0, 0, I(0));
    BPoly r(p.u - (along ? 1 : 0), p.v - (along ? 0 : 1));
    for (unsigned i = 0; i <= r.u; ++i)
        for (unsigned j = 0; j <= r.v; ++j)
            r.at(i, j) = (p.at(i + (along ? 1 : 0), j + (along ? 0 : 1)) - p.at(i, j)) * I(n);
    return r;
}
static I Eval1(std::vector<I> a, double t)
{
    for (size_t k = a.size() - 1; k > 0; --k)
        for (size_t i = 0; i < k; ++i)
            a[i] = a[i] * (I(1) - I(t)) + a[i + 1] * I(t);
    return a[0];
}
BPoly Evaluate(const BPoly &p, double u, double v)
{
    std::vector<I> rows;
    for (unsigned i = 0; i <= p.u; ++i)
    {
        std::vector<I> a;
        for (unsigned j = 0; j <= p.v; ++j)
            a.push_back(p.at(i, j));
        rows.push_back(Eval1(a, v));
    }
    return BPoly(0, 0, Eval1(rows, u));
}
static std::pair<BPoly, BPoly> SplitAt(const BPoly &p, bool along, I fraction)
{
    BPoly a = p, b = p;
    unsigned n = along ? p.u : p.v, lines = along ? p.v : p.u;
    for (unsigned j = 0; j <= lines; ++j)
    {
        std::vector<I> t(n + 1);
        for (unsigned i = 0; i <= n; ++i)
            t[i] = along ? p.at(i, j) : p.at(j, i);
        for (unsigned k = 0; k <= n; ++k)
        {
            if (along)
            {
                a.at(k, j) = t[0];
                b.at(n - k, j) = t[n - k];
            }
            else
            {
                a.at(j, k) = t[0];
                b.at(j, n - k) = t[n - k];
            }
            for (unsigned i = 0; i < n - k; ++i)
                t[i] = t[i] * (I(1) - fraction) + t[i + 1] * fraction;
        }
    }
    return {a, b};
}
std::pair<BPoly, BPoly> Split(const BPoly &p, bool along)
{
    return SplitAt(p, along, I(.5));
}
BV operator+(const BV &a, const BV &b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
BV operator-(const BV &a, const BV &b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
BV operator*(const BV &a, I b)
{
    return {a.x * b, a.y * b, a.z * b};
}
BV operator*(const BV &a, const BPoly &b)
{
    return {a.x * b, a.y * b, a.z * b};
}
BV Constant(const BV &a, IV b)
{
    return {Constant(a.x, b.x), Constant(a.y, b.y), Constant(a.z, b.z)};
}
BPoly Dot(const BV &a, const BV &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
BPoly Dot(const BV &a, IV b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
BV Cross(const BV &a, const BV &b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
BV Derivative(const BV &a, bool u)
{
    return {Derivative(a.x, u), Derivative(a.y, u), Derivative(a.z, u)};
}
BV Evaluate(const BV &a, double u, double v)
{
    return {Evaluate(a.x, u, v), Evaluate(a.y, u, v), Evaluate(a.z, u, v)};
}
IV Range(const BV &a)
{
    return {a.x.range(), a.y.range(), a.z.range()};
}
static std::vector<double> Expanded(const sggk::RealArray &k, const sggk::UIntArray &m)
{
    std::vector<double> out;
    if (k.size() != m.size())
        throw Failure(Status::InvalidInput, "Knot/multiplicity size mismatch");
    for (size_t i = 0; i < k.size(); ++i)
        out.insert(out.end(), m[i], k[i]);
    return out;
}
static std::vector<std::vector<I>> Extraction(const std::vector<double> &k, unsigned degree,
                                              unsigned span)
{
    std::vector<std::vector<I>> out(degree + 1, std::vector<I>(degree + 1));
    if (k[span - degree + 1] == k[span] && k[span + degree] == k[span + 1])
    {
        for (unsigned i = 0; i <= degree; ++i)
            out[i][i] = I(1);
        return out;
    }
    for (unsigned j = 0; j <= degree; ++j)
    {
        std::vector<std::vector<I>> d(degree + 1, std::vector<I>(degree + 1));
        for (unsigned i = 0; i <= degree; ++i)
            d[i][i] = I(1);
        for (unsigned r = 1; r <= degree; ++r)
        {
            double t = r <= degree - j ? k[span] : k[span + 1];
            for (unsigned i = degree; i >= r; --i)
            {
                I alpha = (I(t) - I(k[span - degree + i])) /
                          (I(k[span + 1 + i - r]) - I(k[span - degree + i]));
                for (unsigned b = 0; b <= degree; ++b)
                    d[i][b] = d[i - 1][b] * (I(1) - alpha) + d[i][b] * alpha;
            }
        }
        out[j] = d[degree];
    }
    return out;
}
std::vector<BezierCell> ExtractCells(const sggk::BSplineSurface &s)
{
    auto ku = Expanded(s.KnotsU(), s.MultsU()), kv = Expanded(s.KnotsV(), s.MultsV());
    unsigned p = s.DegreeU(), q = s.DegreeV();
    const auto &cp = s.ControlPoints();
    if (ku.size() != cp.size() + p + 1 || kv.size() != cp.front().size() + q + 1)
        throw Failure(Status::InvalidInput, "Invalid full knot count");
    std::vector<BezierCell> out;
    for (unsigned a = p; a < cp.size(); ++a)
        if (ku[a] < ku[a + 1])
        {
            auto eu = Extraction(ku, p, a);
            for (unsigned b = q; b < cp.front().size(); ++b)
                if (kv[b] < kv[b + 1])
                {
                    if (ku[a + 1] <= s.MinParamU() || ku[a] >= s.MaxParamU() ||
                        kv[b + 1] <= s.MinParamV() || kv[b] >= s.MaxParamV())
                        continue;
                    auto ev = Extraction(kv, q, b);
                    BezierCell cell{{BPoly(p, q), BPoly(p, q), BPoly(p, q)},
                                    ku[a],
                                    ku[a + 1],
                                    kv[b],
                                    kv[b + 1]};
                    if (s.IsRational()) cell.w = BPoly(p, q);
                    for (unsigned i = 0; i <= p; ++i)
                        for (unsigned j = 0; j <= q; ++j)
                        {
                            IV sum{I(0), I(0), I(0)};
                            I denominator(0);
                            for (unsigned k = 0; k <= p; ++k)
                                for (unsigned l = 0; l <= q; ++l)
                                {
                                    I weight = eu[i][k] * ev[j][l];
                                    if (weight.lo == 0 && weight.hi == 0)
                                        continue;
                                    if (s.IsRational())
                                    {
                                        weight = weight * I(s.Weights()[a - p + k][b - q + l]);
                                        denominator = denominator + weight;
                                    }
                                    sum = sum + Exact(Get(cp[a - p + k][b - q + l])) * weight;
                                }
                            cell.p.x.at(i, j) = sum.x;
                            cell.p.y.at(i, j) = sum.y;
                            cell.p.z.at(i, j) = sum.z;
                            if (s.IsRational()) cell.w.at(i, j) = denominator;
                        }
                    // SGK may retain an endpoint one ulp inside the reversed knot
                    // domain. Restrict the polynomial itself, not only its UV labels.
                    auto crop = [&](bool along, double lo, double hi) {
                        double &start = along ? cell.u0 : cell.v0;
                        double &end = along ? cell.u1 : cell.v1;
                        auto cut = [&](I t, bool keepRight) {
                            auto x = SplitAt(cell.p.x, along, t), y = SplitAt(cell.p.y, along, t),
                                 z = SplitAt(cell.p.z, along, t);
                            cell.p = keepRight ? BV{x.second, y.second, z.second}
                                               : BV{x.first, y.first, z.first};
                            auto w = SplitAt(cell.w, along, t);
                            cell.w = keepRight ? w.second : w.first;
                        };
                        if (start < lo)
                        {
                            cut((I(lo) - I(start)) / (I(end) - I(start)), true);
                            start = lo;
                        }
                        if (end > hi)
                        {
                            cut((I(hi) - I(start)) / (I(end) - I(start)), false);
                            end = hi;
                        }
                    };
                    crop(true, s.MinParamU(), s.MaxParamU());
                    crop(false, s.MinParamV(), s.MaxParamV());
                    out.push_back(std::move(cell));
                }
        }
    return out;
}
std::pair<BezierCell, BezierCell> Split(const BezierCell &p, bool u)
{
    auto x = Split(p.p.x, u), y = Split(p.p.y, u), z = Split(p.p.z, u);
    BezierCell a = p, b = p;
    a.p = {x.first, y.first, z.first};
    b.p = {x.second, y.second, z.second};
    auto w = Split(p.w, u);
    a.w = w.first;
    b.w = w.second;
    if (u)
    {
        double t = (p.u0 + p.u1) * 0.5;
        a.u1 = t;
        b.u0 = t;
    }
    else
    {
        double t = (p.v0 + p.v1) * 0.5;
        a.v1 = t;
        b.v0 = t;
    }
    return {a, b};
}
BV Boundary(const BV &a, BoundaryId edge)
{
    auto get = [&](const BPoly &p) {
        bool fixedU = edge == BoundaryId::Umin || edge == BoundaryId::Umax;
        unsigned n = fixedU ? p.v : p.u;
        BPoly r(n, 0);
        for (unsigned i = 0; i <= n; ++i)
            r.at(i, 0) = fixedU ? p.at(edge == BoundaryId::Umin ? 0 : p.u, i)
                                : p.at(i, edge == BoundaryId::Vmin ? 0 : p.v);
        return r;
    };
    return {get(a.x), get(a.y), get(a.z)};
}
I AngleRange(const BPoly &x, const BPoly &y, double ref)
{
    I c = Cos(I(ref)), s = Sin(I(ref));
    BPoly a = x * c + y * s, b = y * c - x * s;
    if (a.range().lo <= 0)
        return I::Whole();
    // With positive denominator Bernstein coefficients, b(t)/a(t) is a
    // convex combination of b_i/a_i. Keep shared radial/height dependence:
    // dividing two independent coordinate ranges was extremely loose on cones.
    double lo = INFINITY, hi = -INFINITY;
    for (size_t i = 0; i < a.c.size(); ++i)
    {
        I r = b.c[i] / a.c[i];
        lo = std::min(lo, r.lo);
        hi = std::max(hi, r.hi);
    }
    return I(ref) + Atan(I(lo, hi));
}
} // namespace bspline2analytic::detail
