#include "PairedValidationInternal.h"

namespace bspline2analytic::detail::paired
{
namespace
{
struct Span
{
    unsigned first;
    double lo, hi;
    std::vector<I> matrix;
};
std::vector<Span> AxisSpans(const sggk::RealArray &knots, const sggk::UIntArray &mults,
                            unsigned degree, size_t count, double lo, double hi)
{
    std::vector<double> full;
    full.reserve(count + degree + 1);
    for (size_t i = 0; i < knots.size(); ++i)
        full.insert(full.end(), mults[i], knots[i]);
    const unsigned n = degree + 1;
    std::vector<Span> spans;
    for (unsigned a = degree; a < count; ++a)
    {
        if (full[a] >= full[a + 1] || full[a + 1] <= lo || full[a] >= hi)
            continue;
        Span span{a - degree, full[a], full[a + 1], std::vector<I>(n * n)};
        if (full[a - degree + 1] == full[a] && full[a + degree] == full[a + 1])
        {
            for (unsigned i = 0; i < n; ++i)
                span.matrix[i * n + i] = I(1);
        }
        else
        {
            // The same interval blossom recurrence as the shared extractor,
            // computed once per active one-dimensional span.
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
                        I complement = I(1) - alpha;
                        for (unsigned b = 0; b < n; ++b)
                            work[i * n + b] =
                                work[(i - 1) * n + b] * complement + work[i * n + b] * alpha;
                    }
                }
                for (unsigned k = 0; k < n; ++k)
                    span.matrix[j * n + k] = work[degree * n + k];
            }
        }
        spans.push_back(std::move(span));
    }
    return spans;
}
bool Zero(I a)
{
    return a.lo == 0 && a.hi == 0;
}
std::array<double, 4> Domain(const sggk::BSplineSurface &s)
{
    return {s.MinParamU(), s.MaxParamU(), s.MinParamV(), s.MaxParamV()};
}
} // namespace

std::vector<BezierCell> ExtractPolynomialCells(const sggk::BSplineSurface &s,
                                               PairedValidationReport &report)
{
    // Caller validates dimensions, continuity, degree and cell budgets first.
    if (s.IsRational())
        throw Failure(Status::Unsupported, "Paired source extraction requires polynomial input");
    const auto &cp = s.ControlPoints();
    unsigned p = s.DegreeU(), q = s.DegreeV();
    auto us = AxisSpans(s.KnotsU(), s.MultsU(), p, cp.size(), s.MinParamU(), s.MaxParamU());
    auto vs = AxisSpans(s.KnotsV(), s.MultsV(), q, cp.front().size(), s.MinParamV(), s.MaxParamV());
    report.uExtractionMatrices += static_cast<unsigned>(us.size());
    report.vExtractionMatrices += static_cast<unsigned>(vs.size());
    std::vector<BezierCell> out;
    out.reserve(us.size() * vs.size());
    // Transform every source column along U once and reuse it for all V spans.
    const size_t columns = cp.front().size();
    std::vector<IV> alongU((p + 1) * columns);
    for (const auto &u : us)
    {
        for (unsigned i = 0; i <= p; ++i)
            for (size_t l = 0; l < columns; ++l)
            {
                IV sum{};
                for (unsigned k = 0; k <= p; ++k)
                {
                    I weight = u.matrix[i * (p + 1) + k];
                    if (!Zero(weight))
                        sum = sum + Exact(Get(cp[u.first + k][l])) * weight;
                }
                alongU[i * columns + l] = sum;
            }
        for (const auto &v : vs)
        {
            BezierCell cell{{BPoly(p, q), BPoly(p, q), BPoly(p, q)}, u.lo, u.hi, v.lo, v.hi};
            for (unsigned i = 0; i <= p; ++i)
                for (unsigned j = 0; j <= q; ++j)
                {
                    IV sum{};
                    for (unsigned l = 0; l <= q; ++l)
                    {
                        I weight = v.matrix[j * (q + 1) + l];
                        if (!Zero(weight))
                            sum = sum + alongU[i * columns + v.first + l] * weight;
                    }
                    cell.p.x.at(i, j) = sum.x;
                    cell.p.y.at(i, j) = sum.y;
                    cell.p.z.at(i, j) = sum.z;
                }
            // Preserve partial end spans, including SDK reversal's one-ulp trims.
            if (u.lo < s.MinParamU() || u.hi > s.MaxParamU() || v.lo < s.MinParamV() ||
                v.hi > s.MaxParamV())
                cell = Crop(cell, std::max(u.lo, s.MinParamU()), std::min(u.hi, s.MaxParamU()),
                            std::max(v.lo, s.MinParamV()), std::min(v.hi, s.MaxParamV()));
            out.push_back(std::move(cell));
        }
    }
    return out;
}

bool SourceCache::Matches(const sggk::BSplineSurface &s) const
{
    if (cells.empty() || s.IsRational() || degreeU != s.DegreeU() || degreeV != s.DegreeV() ||
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
            V point = Get(cp[i][j]), old = controls[i][j];
            if (point.x != old.x || point.y != old.y || point.z != old.z)
                return false;
        }
    }
    return true;
}
void SourceCache::Store(const sggk::BSplineSurface &s, const std::vector<BezierCell> &extracted)
{
    SourceCache next;
    next.degreeU = s.DegreeU();
    next.degreeV = s.DegreeV();
    next.domain = Domain(s);
    next.knotsU = s.KnotsU();
    next.knotsV = s.KnotsV();
    next.multsU = s.MultsU();
    next.multsV = s.MultsV();
    for (const auto &row : s.ControlPoints())
    {
        std::vector<V> points;
        points.reserve(row.size());
        for (const auto &p : row)
            points.push_back(Get(p));
        next.controls.push_back(std::move(points));
    }
    next.cells = extracted;
    *this = std::move(next);
}
} // namespace bspline2analytic::detail::paired
