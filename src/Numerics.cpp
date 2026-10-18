/**
 * @file Numerics.cpp
 * @brief 精确识别共用的向量、线性求解与 Bernstein 数值工具。
 * 负责基础几何计算与浮点误差界支持，区别于 IntervalBounds 的向外舍入区间算术。
 */
#include "Numerics.h"
#include <Geometry/3D/Curve/BezierCurve3D.h>
#include <Geometry/3D/Surface/BezierSurface.h>
#include <numeric>

namespace bspline2analytic::detail
{
double Roundoff(double scale) noexcept
{
    return 256 * std::numeric_limits<double>::epsilon() * std::max(1.0, scale);
}
double Poly::maxAbs() const
{
    double r = 0;
    for (double a : c)
        r = std::max(r, std::abs(a));
    return r;
}
Poly operator+(const Poly &a, const Poly &b)
{
    if (a.u != b.u || a.v != b.v)
        throw std::logic_error("Bernstein degree mismatch");
    Poly r = a;
    for (std::size_t i = 0; i < r.c.size(); ++i)
        r.c[i] += b.c[i];
    return r;
}
Poly operator*(const Poly &a, double s)
{
    Poly r = a;
    for (auto &c : r.c)
        c *= s;
    return r;
}
Poly operator-(const Poly &a, const Poly &b)
{
    return a + b * (-1);
}
static double Binomial(unsigned n, unsigned k)
{
    double r = 1;
    k = std::min(k, n - k);
    for (unsigned i = 1; i <= k; ++i)
        r *= double(n - k + i) / i;
    return r;
}
Poly operator*(const Poly &a, const Poly &b)
{
    Poly r(a.u + b.u, a.v + b.v);
    for (unsigned i = 0; i <= a.u; ++i)
        for (unsigned k = 0; k <= b.u; ++k)
        {
            const double wu = Binomial(a.u, i) * Binomial(b.u, k) / Binomial(r.u, i + k);
            for (unsigned j = 0; j <= a.v; ++j)
                for (unsigned l = 0; l <= b.v; ++l)
                    r.at(i + k, j + l) += a.at(i, j) * b.at(k, l) * wu * Binomial(a.v, j) *
                                          Binomial(b.v, l) / Binomial(r.v, j + l);
        }
    return r;
}
std::pair<Poly, Poly> Split(const Poly &p, bool alongU)
{
    Poly a(p.u, p.v), b(p.u, p.v);
    unsigned n = alongU ? p.u : p.v, lines = alongU ? p.v : p.u;
    for (unsigned l = 0; l <= lines; ++l)
    {
        std::vector<double> t(n + 1);
        for (unsigned i = 0; i <= n; ++i)
            t[i] = alongU ? p.at(i, l) : p.at(l, i);
        for (unsigned k = 0; k <= n; ++k)
        {
            if (alongU)
            {
                a.at(k, l) = t[0];
                b.at(n - k, l) = t[n - k];
            }
            else
            {
                a.at(l, k) = t[0];
                b.at(l, n - k) = t[n - k];
            }
            for (unsigned i = 0; i < n - k; ++i)
                t[i] = (t[i] + t[i + 1]) * 0.5;
        }
    }
    return {a, b};
}
std::pair<Homogeneous, Homogeneous> Split(const Homogeneous &h, bool u)
{
    auto x = Split(h.x, u), y = Split(h.y, u), z = Split(h.z, u), w = Split(h.w, u);
    return {{x.first, y.first, z.first, w.first}, {x.second, y.second, z.second, w.second}};
}
static Homogeneous Make(unsigned u, unsigned v)
{
    return {Poly(u, v), Poly(u, v), Poly(u, v), Poly(u, v)};
}
static void Set(Homogeneous &h, unsigned i, unsigned j, V p, double w, V origin, double scale,
                double ws)
{
    p = (p - origin) / scale;
    w /= ws;
    h.x.at(i, j) = p.x * w;
    h.y.at(i, j) = p.y * w;
    h.z.at(i, j) = p.z * w;
    h.w.at(i, j) = w;
}
std::vector<Homogeneous> CurvePieces(const sggk::BSplineCurve3D &c, V origin, double scale)
{
    std::vector<Homogeneous> out;
    for (const auto &b : c.ToBezierCurves())
    {
        const auto &pts = b->ControlPoints();
        const auto &weights = b->Weights();
        auto h = Make(static_cast<unsigned>(pts.size() - 1), 0);
        double ws = 1;
        for (double w : weights)
            ws = std::max(ws, std::abs(w));
        if (!weights.empty() && weights.front() < 0)
            ws = -ws;
        for (unsigned i = 0; i < pts.size(); ++i)
            Set(h, i, 0, Get(pts[i]), weights.empty() ? 1 : weights[i], origin, scale, ws);
        out.push_back(std::move(h));
    }
    return out;
}
std::vector<Homogeneous> SurfacePieces(const sggk::BSplineSurface &s, V origin, double scale)
{
    std::vector<Homogeneous> out;
    for (const auto &row : s.ToBezierSurfaces())
        for (const auto &b : row)
        {
            const auto &pts = b->ControlPoints();
            const auto &weights = b->Weights();
            auto h = Make(static_cast<unsigned>(pts.size() - 1),
                          static_cast<unsigned>(pts.front().size() - 1));
            double ws = 1;
            for (const auto &rw : weights)
                for (double w : rw)
                    ws = std::max(ws, std::abs(w));
            if (!weights.empty() && weights.front().front() < 0)
                ws = -ws;
            for (unsigned i = 0; i < pts.size(); ++i)
                for (unsigned j = 0; j < pts[i].size(); ++j)
                    Set(h, i, j, Get(pts[i][j]), weights.empty() ? 1 : weights[i][j], origin, scale,
                        ws);
            out.push_back(std::move(h));
        }
    return out;
}
Poly DotNumerator(const Homogeneous &h, V d)
{
    return h.x * d.x + h.y * d.y + h.z * d.z;
}
Poly SquaredNorm(const Homogeneous &h)
{
    return h.x * h.x + h.y * h.y + h.z * h.z;
}

// 只有三个未知量的列主元 Householder QR，避免平方放大条件数。
bool LeastSquares3(std::vector<std::array<double, 3>> a, std::vector<double> b, double tol,
                   std::array<double, 3> &result)
{
    if (a.size() < 3 || a.size() != b.size())
        return false;
    std::array<unsigned, 3> perm{0, 1, 2};
    double first = 0;
    for (unsigned k = 0; k < 3; ++k)
    {
        unsigned pivot = k;
        double best = -1;
        for (unsigned j = k; j < 3; ++j)
        {
            double s = 0;
            for (std::size_t i = k; i < a.size(); ++i)
                s = std::hypot(s, a[i][j]);
            if (s > best)
            {
                best = s;
                pivot = j;
            }
        }
        if (k == 0)
            first = best;
        if (!(best > tol * first) || !std::isfinite(best))
            return false;
        for (auto &row : a)
            std::swap(row[k], row[pivot]);
        std::swap(perm[k], perm[pivot]);
        double alpha = -std::copysign(best, a[k][k]);
        std::vector<double> v(a.size() - k);
        for (std::size_t i = k; i < a.size(); ++i)
            v[i - k] = a[i][k];
        v[0] -= alpha;
        double vn = 0;
        for (double x : v)
            vn = std::hypot(vn, x);
        for (auto &x : v)
            x /= vn;
        for (unsigned j = k; j < 3; ++j)
        {
            double d = 0;
            for (std::size_t i = k; i < a.size(); ++i)
                d += v[i - k] * a[i][j];
            for (std::size_t i = k; i < a.size(); ++i)
                a[i][j] -= 2 * v[i - k] * d;
        }
        double d = 0;
        for (std::size_t i = k; i < a.size(); ++i)
            d += v[i - k] * b[i];
        for (std::size_t i = k; i < a.size(); ++i)
            b[i] -= 2 * v[i - k] * d;
    }
    std::array<double, 3> x{};
    for (int k = 2; k >= 0; --k)
    {
        double d = b[k];
        for (unsigned j = k + 1; j < 3; ++j)
            d -= a[k][j] * x[j];
        x[k] = d / a[k][k];
    }
    for (unsigned k = 0; k < 3; ++k)
        result[perm[k]] = x[k];
    return true;
}
std::vector<double> SampleParameters(const sggk::RealArray &knots, double lo, double hi, unsigned n)
{
    std::vector<double> out{lo};
    double left = lo;
    for (double k : knots)
    {
        if (k <= left || k >= hi)
            continue;
        for (unsigned j = 1; j <= n; ++j)
            out.push_back(left + (k - left) * double(j) / n);
        left = k;
    }
    for (unsigned j = 1; j <= n; ++j)
        out.push_back(left + (hi - left) * double(j) / n);
    return out;
}
std::vector<double> SampleParameters(const sggk::BSplineCurve3D &c, unsigned n)
{
    return SampleParameters(c.Knots(), c.MinParam(), c.MaxParam(), n);
}
void ValidateOptions(const RecognitionOptions &o)
{
    if (!(o.distanceTolerance > 0) || !std::isfinite(o.distanceTolerance) ||
        !(o.rankTolerance > 0 && o.rankTolerance < 1) ||
        !(o.angularTolerance > 0 && o.angularTolerance < 1) || o.samplesPerSpan < 2 ||
        o.samplesPerSpan > 64 || o.maxSubdivisionDepth > 24 || o.maxValidationPieces == 0)
        throw Failure(Status::InvalidInput, "Invalid recognition tolerances or validation limits");
}
void ValidateSurface(const sggk::BSplineSurface &s)
{
    if (!std::isfinite(s.MinParamU()) || !std::isfinite(s.MaxParamU()) ||
        !std::isfinite(s.MinParamV()) || !std::isfinite(s.MaxParamV()) ||
        !(s.MinParamU() < s.MaxParamU() && s.MinParamV() < s.MaxParamV()))
        throw Failure(Status::InvalidInput, "Invalid parameter domain");
    const auto &p = s.ControlPoints();
    const auto &w = s.Weights();
    if (p.size() < 2 || p.front().size() < 2)
        throw Failure(Status::Degenerate, "Insufficient surface control points");
    if (!w.empty() && w.size() != p.size())
        throw Failure(Status::InvalidInput, "Weight dimensions differ from control net");
    for (std::size_t i = 0; i < p.size(); ++i)
    {
        if (p[i].size() != p.front().size() || (!w.empty() && w[i].size() != p[i].size()))
            throw Failure(Status::InvalidInput, "Nonrectangular control net");
        for (std::size_t j = 0; j < p[i].size(); ++j)
            if (!Finite(Get(p[i][j])) || (!w.empty() && !std::isfinite(w[i][j])))
                throw Failure(Status::InvalidInput, "Nonfinite surface data");
    }
}
} // namespace bspline2analytic::detail
