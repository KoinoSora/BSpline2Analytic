#include "MotherDistanceInternal.h"

namespace bspline2analytic::detail::mother
{
namespace
{
double Up(double x)
{
    return std::nextafter(x, INFINITY);
}
V Abs(V a)
{
    return {std::abs(a.x), std::abs(a.y), std::abs(a.z)};
}
V Max(V a, V b)
{
    return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}
V Round(V a)
{
    return {RoundedError(a.x), RoundedError(a.y), RoundedError(a.z)};
}
double DotError(double A, double e, double M, double r, double gamma)
{
    // Coefficient uncertainty, input uncertainty, then ordinary dot-product rounding.
    return ErrorAdd(
        ErrorAdd(ErrorMul(ErrorAdd(A, e), r), ErrorMul(e, M)),
        ErrorAdd(ErrorMul(ErrorMul(gamma, A), M), 64 * std::numeric_limits<double>::denorm_min()));
}
V DotError(const ExtractionSpan &s, unsigned row, V m, V r)
{
    double A = s.rowNorm[row], e = s.rowError[row], g = Gamma(2 * (s.degree + 1));
    return {DotError(A, e, m.x, r.x, g), DotError(A, e, m.y, r.y, g), DotError(A, e, m.z, r.z, g)};
}
std::array<double, 4> Domain(const sggk::BSplineSurface &s)
{
    return {s.MinParamU(), s.MaxParamU(), s.MinParamV(), s.MaxParamV()};
}
std::vector<ExtractionSpan> AxisSpans(const sggk::RealArray &knots, const sggk::UIntArray &mults,
                                      unsigned degree, size_t count, double lo, double hi)
{
    std::vector<double> full;
    for (size_t i = 0; i < knots.size(); ++i)
        full.insert(full.end(), mults[i], knots[i]);
    unsigned n = degree + 1;
    std::vector<ExtractionSpan> spans;
    for (unsigned a = degree; a < count; ++a)
    {
        if (full[a] >= full[a + 1] || full[a + 1] <= lo || full[a] >= hi)
            continue;
        ExtractionSpan s;
        s.first = a - degree;
        s.degree = degree;
        s.lo = std::max(full[a], lo);
        s.hi = std::min(full[a + 1], hi);
        s.enclosure.resize(n * n);
        if (full[a - degree + 1] == full[a] && full[a + degree] == full[a + 1])
            for (unsigned i = 0; i < n; ++i)
                s.enclosure[i * n + i] = I(1);
        else
        {
            std::vector<I> work(n * n);
            for (unsigned j = 0; j < n; ++j)
            {
                std::fill(work.begin(), work.end(), I(0));
                for (unsigned i = 0; i < n; ++i)
                    work[i * n + i] = I(1);
                for (unsigned r = 1; r <= degree; ++r)
                {
                    double t = r <= degree - j ? full[a] : full[a + 1];
                    for (unsigned i = degree; i >= r; --i)
                    {
                        I alpha = (I(t) - I(full[a - degree + i])) /
                                  (I(full[a + 1 + i - r]) - I(full[a - degree + i]));
                        for (unsigned k = 0; k < n; ++k)
                            work[i * n + k] =
                                work[(i - 1) * n + k] * (I(1) - alpha) + work[i * n + k] * alpha;
                    }
                }
                for (unsigned k = 0; k < n; ++k)
                    s.enclosure[j * n + k] = work[degree * n + k];
            }
        }
        // Restrict the small matrix, not every 3D cell. Do not snap one-ulp trims.
        if (s.lo != full[a] || s.hi != full[a + 1])
        {
            I l = s.lo == full[a] ? I(0) : (I(s.lo) - I(full[a])) / (I(full[a + 1]) - I(full[a]));
            I h =
                s.hi == full[a + 1] ? I(1) : (I(s.hi) - I(full[a])) / (I(full[a + 1]) - I(full[a]));
            for (unsigned k = 0; k < n; ++k)
            {
                BPoly col(degree, 0);
                for (unsigned i = 0; i < n; ++i)
                    col.c[i] = s.enclosure[i * n + k];
                col = paired::Restrict(col, true, l, h);
                for (unsigned i = 0; i < n; ++i)
                    s.enclosure[i * n + k] = col.c[i];
            }
        }
        s.weights.resize(n * n);
        s.rowNorm.resize(n);
        s.rowError.resize(n);
        for (unsigned i = 0; i < n; ++i)
            for (unsigned k = 0; k < n; ++k)
            {
                I v = s.enclosure[i * n + k];
                double mid = v.lo * .5 + v.hi * .5;
                if (!v.finite())
                    throw Failure(Status::VerificationFailed, "Nonfinite extraction matrix");
                s.weights[i * n + k] = mid;
                s.rowNorm[i] = ErrorAdd(s.rowNorm[i], std::abs(mid));
                s.rowError[i] = ErrorAdd(s.rowError[i], (v - I(mid)).absMax());
            }
        spans.push_back(std::move(s));
    }
    return spans;
}
BezierCell Enclose(const SourceCell &s)
{
    if (s.tight)
        return *s.tight;
    BezierCell b{{BPoly(s.x.u, s.x.v), BPoly(s.x.u, s.x.v), BPoly(s.x.u, s.x.v)},
                 s.domain[0],
                 s.domain[1],
                 s.domain[2],
                 s.domain[3]};
    for (size_t i = 0; i < s.x.c.size(); ++i)
    {
        b.p.x.c[i] = I(s.x.c[i]) + I(-s.error.x, s.error.x);
        b.p.y.c[i] = I(s.y.c[i]) + I(-s.error.y, s.error.y);
        b.p.z.c[i] = I(s.z.c[i]) + I(-s.error.z, s.error.z);
    }
    return b;
}
} // namespace
double ErrorAdd(double a, double b)
{
    return a == 0 ? b : b == 0 ? a : Up(a + b);
}
double ErrorMul(double a, double b)
{
    return a == 0 || b == 0 ? 0 : Up(a * b);
}
double Gamma(unsigned n)
{
    double nu = n * 0x1p-53;
    return Up(nu / std::nextafter(1 - nu, 0.));
}
double RoundedError(double x)
{
    return ErrorAdd(ErrorMul(Gamma(1), std::abs(x)), std::numeric_limits<double>::denorm_min());
}
double ErrorNorm(V a)
{
    return Up(
        std::sqrt(ErrorAdd(ErrorAdd(ErrorMul(a.x, a.x), ErrorMul(a.y, a.y)), ErrorMul(a.z, a.z))));
}

bool SourceCache::Matches(const sggk::BSplineSurface &s) const
{
    if (cells.empty() || s.IsRational() || p != s.DegreeU() || q != s.DegreeV() ||
        domain != Domain(s) || knotsU != s.KnotsU() || knotsV != s.KnotsV() ||
        multsU != s.MultsU() || multsV != s.MultsV())
        return false;
    const auto &cp = s.ControlPoints();
    if (controls.size() != cp.size())
        return false;
    for (size_t i = 0; i < cp.size(); ++i)
    {
        if (controls[i].size() != cp[i].size())
            return false;
        for (size_t j = 0; j < cp[i].size(); ++j)
        {
            V a = detail::Get(cp[i][j]), b = controls[i][j];
            if (a.x != b.x || a.y != b.y || a.z != b.z)
                return false;
        }
    }
    return true;
}
void SourceCache::Initialize(const sggk::BSplineSurface &s, MotherDistanceValidationReport &out)
{
    SourceCache next;
    next.p = s.DegreeU();
    next.q = s.DegreeV();
    next.domain = Domain(s);
    next.knotsU = s.KnotsU();
    next.knotsV = s.KnotsV();
    next.multsU = s.MultsU();
    next.multsV = s.MultsV();
    V lo = detail::Get(s.ControlPoints().front().front()), hi = lo;
    for (const auto &row : s.ControlPoints())
    {
        std::vector<V> points;
        for (const auto &point : row)
        {
            V v = detail::Get(point);
            points.push_back(v);
            lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
            hi = Max(hi, v);
        }
        next.controls.push_back(std::move(points));
    }
    next.origin = lo * .5 + hi * .5;
    {
        Measure time(out.sourceEnvelopeMs);
        for (const auto &row : next.controls)
            for (V v : row)
            {
                V c = v - next.origin;
                next.centered.push_back(c);
                next.initialError = Max(next.initialError, Round(c));
                next.magnitude = Max(next.magnitude, Abs(c));
            }
    }
    {
        Measure time(out.extractionPlanMs);
        next.us = AxisSpans(next.knotsU, next.multsU, next.p, next.controls.size(), next.domain[0],
                            next.domain[1]);
        next.vs = AxisSpans(next.knotsV, next.multsV, next.q, next.controls.front().size(),
                            next.domain[2], next.domain[3]);
        ++out.extractionPlanBuilds;
    }
    next.strips.resize(next.us.size());
    next.cells.resize(next.us.size() * next.vs.size());
    if (next.cells.empty())
        throw Failure(Status::InvalidInput, "Empty source domain");
    *this = std::move(next);
}
const SourceCell &SourceCache::Get(size_t id, MotherDistanceValidationReport &out)
{
    if (cells.at(id))
    {
        ++out.cachedCellHits;
        return *cells[id];
    }
    size_t ui = id / vs.size(), vi = id % vs.size(), columns = controls.front().size();
    const auto &u = us[ui], &v = vs[vi];
    if (!strips[ui])
    {
        SourceStrip strip;
        strip.points.resize((p + 1) * columns);
        strip.error.resize(p + 1);
        {
            Measure time(out.sourceContractionMs);
            for (unsigned i = 0; i <= p; ++i)
                for (size_t l = 0; l < columns; ++l)
                {
                    V sum{};
                    for (unsigned k = 0; k <= p; ++k)
                        sum = sum +
                              centered[(u.first + k) * columns + l] * u.weights[i * (p + 1) + k];
                    strip.points[i * columns + l] = sum;
                }
        }
        {
            Measure time(out.sourceEnvelopeMs);
            for (unsigned i = 0; i <= p; ++i)
                strip.error[i] = DotError(u, i, magnitude, initialError);
        }
        strips[ui] = std::move(strip);
    }
    const auto &strip = *strips[ui];
    SourceCell cell{Poly(p, q), Poly(p, q), Poly(p, q), {}, {u.lo, u.hi, v.lo, v.hi}, {}};
    {
        Measure time(out.sourceContractionMs);
        for (unsigned i = 0; i <= p; ++i)
            for (unsigned j = 0; j <= q; ++j)
            {
                V sum{};
                for (unsigned l = 0; l <= q; ++l)
                    sum =
                        sum + strip.points[i * columns + v.first + l] * v.weights[j * (q + 1) + l];
                cell.x.at(i, j) = sum.x;
                cell.y.at(i, j) = sum.y;
                cell.z.at(i, j) = sum.z;
            }
    }
    {
        Measure time(out.sourceEnvelopeMs);
        for (unsigned i = 0; i <= p; ++i)
        {
            V m{};
            for (unsigned l = 0; l <= q; ++l)
                m = Max(m, Abs(strip.points[i * columns + v.first + l]));
            for (unsigned j = 0; j <= q; ++j)
                cell.error = Max(cell.error, DotError(v, j, m, strip.error[i]));
        }
    }
    double error = ErrorNorm(cell.error);
    if (!std::isfinite(error))
        throw Failure(Status::VerificationFailed, "Nonfinite extraction envelope");
    out.maxExtractionError = std::max(out.maxExtractionError, error);
    ++out.extractedCells;
    cells[id] = std::move(cell);
    return *cells[id];
}
const SourceCell &SourceCache::Tighten(size_t id, MotherDistanceValidationReport &out)
{
    Get(id, out);
    auto &cell = *cells[id];
    if (cell.tight)
        return cell;
    Measure time(out.sourceContractionMs);
    const auto &u = us[id / vs.size()], &v = vs[id % vs.size()];
    std::vector<IV> strip((p + 1) * (q + 1));
    for (unsigned i = 0; i <= p; ++i)
        for (unsigned l = 0; l <= q; ++l)
            for (unsigned k = 0; k <= p; ++k)
                strip[i * (q + 1) + l] =
                    strip[i * (q + 1) + l] +
                    (Exact(controls[u.first + k][v.first + l]) - Exact(origin)) *
                        u.enclosure[i * (p + 1) + k];
    auto b = Enclose(cell);
    cell.error = {};
    for (unsigned i = 0; i <= p; ++i)
        for (unsigned j = 0; j <= q; ++j)
        {
            IV sum{};
            for (unsigned l = 0; l <= q; ++l)
                sum = sum + strip[i * (q + 1) + l] * v.enclosure[j * (q + 1) + l];
            b.p.x.at(i, j) = sum.x;
            b.p.y.at(i, j) = sum.y;
            b.p.z.at(i, j) = sum.z;
            V mid{sum.x.lo * .5 + sum.x.hi * .5, sum.y.lo * .5 + sum.y.hi * .5,
                  sum.z.lo * .5 + sum.z.hi * .5};
            cell.x.at(i, j) = mid.x;
            cell.y.at(i, j) = mid.y;
            cell.z.at(i, j) = mid.z;
            cell.error = Max(cell.error, {(sum.x - I(mid.x)).absMax(), (sum.y - I(mid.y)).absMax(),
                                          (sum.z - I(mid.z)).absMax()});
        }
    cell.tight = std::move(b);
    boundary.reset();
    ++out.tightenedCells;
    return cell;
}
const paired::PreparedSource &SourceCache::PrepareBoundary(MotherDistanceValidationReport &out,
                                                           bool tight)
{
    if (boundary && !tight)
        return *boundary;
    paired::PreparedSource next;
    next.origin = origin;
    for (size_t id = 0; id < cells.size(); ++id)
    {
        size_t i = id / vs.size(), j = id % vs.size();
        if (i && i + 1 < us.size() && j && j + 1 < vs.size())
            continue;
        next.cells.push_back(Enclose(tight ? Tighten(id, out) : Get(id, out)));
    }
    for (unsigned edge = 0; edge < 4; ++edge)
        for (const auto &cell : next.cells)
        {
            bool fixedU = edge < 2;
            if (!(edge == 0   ? cell.u0 == domain[0]
                  : edge == 1 ? cell.u1 == domain[1]
                  : edge == 2 ? cell.v0 == domain[2]
                              : cell.v1 == domain[3]))
                continue;
            double begin = fixedU ? cell.v0 : cell.u0, end = fixedU ? cell.v1 : cell.u1;
            BV c = Boundary(cell.p, static_cast<BoundaryId>(edge));
            BV tangent = Derivative(c, true) * (I(1) / (I(end) - I(begin)));
            next.boundaries[edge].push_back({begin, end, std::move(tangent), {}});
        }
    boundary = std::move(next);
    return *boundary;
}
} // namespace bspline2analytic::detail::mother
