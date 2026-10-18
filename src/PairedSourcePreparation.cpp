/**
 * @file PairedSourcePreparation.cpp
 * @brief 配对验证器的源面区间提取及边界缓存。
 * 按 U/V 分离提取复用矩阵和条带，源内容精确匹配才可复用，保留端点微小裁剪。
 */
#include "PairedValidationInternal.h"
#include <chrono>

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
            // 与共享提取器相同的区间多仿射递推，每个活动一维结点区间只计算一次。
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

static std::vector<BezierCell> ExtractPrepared(const sggk::BSplineSurface &s,
                                               PairedValidationReport &report, V origin,
                                               Boundaries *boundaries)
{
    // 调用方必须先验证维数、连续性、次数和单元预算。
    if (s.IsRational())
        throw Failure(Status::Unsupported, "Paired source extraction requires polynomial input");
    const auto &cp = s.ControlPoints();
    unsigned p = s.DegreeU(), q = s.DegreeV();
    using Clock = std::chrono::steady_clock;
    auto elapsed = [](Clock::time_point start) {
        return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    };
    auto start = Clock::now();
    std::vector<std::vector<IV>> points(cp.size(), std::vector<IV>(cp.front().size()));
    for (size_t i = 0; i < cp.size(); ++i)
        for (size_t j = 0; j < cp[i].size(); ++j)
            points[i][j] = Exact(Get(cp[i][j])) - Exact(origin);
    report.sourceCenterMs += elapsed(start);
    start = Clock::now();
    auto us = AxisSpans(s.KnotsU(), s.MultsU(), p, cp.size(), s.MinParamU(), s.MaxParamU());
    auto vs = AxisSpans(s.KnotsV(), s.MultsV(), q, cp.front().size(), s.MinParamV(), s.MaxParamV());
    report.uExtractionMatrices += static_cast<unsigned>(us.size());
    report.vExtractionMatrices += static_cast<unsigned>(vs.size());
    std::vector<BezierCell> out;
    out.reserve(us.size() * vs.size());
    // 每列源控制点沿 U 只变换一次，结果复用于全部 V 结点区间。
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
                        sum = sum + points[u.first + k][l] * weight;
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
            // 保留端部不完整结点区间，包括 SDK 反向操作产生的一个 ulp 级裁剪。
            if (u.lo < s.MinParamU() || u.hi > s.MaxParamU() || v.lo < s.MinParamV() ||
                v.hi > s.MaxParamV())
                cell = Crop(cell, std::max(u.lo, s.MinParamU()), std::min(u.hi, s.MaxParamU()),
                            std::max(v.lo, s.MinParamV()), std::min(v.hi, s.MaxParamV()));
            out.push_back(std::move(cell));
        }
    }
    report.sourceExtractionMs += elapsed(start);
    if (boundaries)
    {
        start = Clock::now();
        Products products;
        // 曲面认证已需要张量提取，因此复用最外层控制行；不重复张量收缩，也不构造二维导数。
        for (unsigned edge = 0; edge < 4; ++edge)
            for (const auto &cell : out)
            {
                bool fixedU = edge < 2, high = (edge % 2) != 0;
                bool outer = edge == 0   ? cell.u0 == out.front().u0
                             : edge == 1 ? cell.u1 == out.back().u1
                             : edge == 2 ? cell.v0 == out.front().v0
                                         : cell.v1 == out.back().v1;
                if (!outer)
                    continue;
                unsigned fixedDegree = fixedU ? cell.p.x.u : cell.p.x.v;
                unsigned varyingDegree = fixedU ? cell.p.x.v : cell.p.x.u;
                double fixedLo = fixedU ? cell.u0 : cell.v0, fixedHi = fixedU ? cell.u1 : cell.v1;
                double value = fixedU ? (high ? s.MaxParamU() : s.MinParamU())
                                      : (high ? s.MaxParamV() : s.MinParamV());
                double begin = fixedU ? cell.v0 : cell.u0, end = fixedU ? cell.v1 : cell.u1;
                I t = value == fixedLo   ? I(0)
                      : value == fixedHi ? I(1)
                                         : (I(value) - I(fixedLo)) / (I(fixedHi) - I(fixedLo));
                BV c{BPoly(varyingDegree, 0), BPoly(varyingDegree, 0), BPoly(varyingDegree, 0)},
                    d = c;
                auto get = [&](unsigned i, unsigned j) {
                    unsigned u = fixedU ? i : j, v = fixedU ? j : i;
                    return IV{cell.p.x.at(u, v), cell.p.y.at(u, v), cell.p.z.at(u, v)};
                };
                auto eval = [t](std::vector<IV> coefficients) {
                    for (size_t n = coefficients.size(); n > 1; --n)
                        for (size_t i = 0; i + 1 < n; ++i)
                            coefficients[i] =
                                coefficients[i] * (I(1) - t) + coefficients[i + 1] * t;
                    return coefficients.front();
                };
                for (unsigned j = 0; j <= varyingDegree; ++j)
                {
                    IV point{}, derivative{};
                    if ((t.lo == 0 && t.hi == 0) || (t.lo == 1 && t.hi == 1))
                    {
                        bool last = t.lo == 1;
                        point = get(last ? fixedDegree : 0, j);
                        if (fixedDegree)
                            derivative = last ? get(fixedDegree, j) - get(fixedDegree - 1, j)
                                              : get(1, j) - get(0, j);
                    }
                    else
                    {
                        // SDK 的一个 ulp 级端点偏移按精确参数保留，不能以此为由吸附边界。
                        std::vector<IV> row(fixedDegree + 1), differences(fixedDegree);
                        for (unsigned i = 0; i <= fixedDegree; ++i)
                            row[i] = get(i, j);
                        for (unsigned i = 0; i < fixedDegree; ++i)
                            differences[i] = row[i + 1] - row[i];
                        point = eval(std::move(row));
                        if (fixedDegree)
                            derivative = eval(std::move(differences));
                    }
                    derivative = derivative * (I(fixedDegree) / (I(fixedHi) - I(fixedLo)));
                    c.x.c[j] = point.x;
                    c.y.c[j] = point.y;
                    c.z.c[j] = point.z;
                    d.x.c[j] = derivative.x;
                    d.y.c[j] = derivative.y;
                    d.z.c[j] = derivative.z;
                }
                BV tangent = Derivative(c, true) * (I(1) / (I(end) - I(begin)));
                BV normal = fixedU ? products.Cross(d, tangent) : products.Cross(tangent, d);
                (*boundaries)[edge].push_back({begin, end, std::move(tangent), std::move(normal)});
            }
        report.sourceBoundaryMs += elapsed(start);
        ++report.sourceBoundaryBuilds;
    }
    return out;
}

std::vector<BezierCell> ExtractPolynomialCells(const sggk::BSplineSurface &s,
                                               PairedValidationReport &report)
{
    return ExtractPrepared(s, report, {}, nullptr);
}
PreparedSource PrepareSource(const sggk::BSplineSurface &s, PairedValidationReport &report)
{
    PreparedSource out;
    V lo = Get(s.ControlPoints().front().front()), hi = lo;
    for (const auto &row : s.ControlPoints())
        for (const auto &point : row)
        {
            V p = Get(point);
            lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
            hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
        }
    out.origin = lo * .5 + hi * .5;
    out.cells = ExtractPrepared(s, report, out.origin, &out.boundaries);
    return out;
}

bool SourceCache::Matches(const sggk::BSplineSurface &s) const
{
    if (prepared.cells.empty() || s.IsRational() || degreeU != s.DegreeU() ||
        degreeV != s.DegreeV() || domain != Domain(s) || knotsU != s.KnotsU() ||
        knotsV != s.KnotsV() || multsU != s.MultsU() || multsV != s.MultsV())
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
void SourceCache::Store(const sggk::BSplineSurface &s, PreparedSource &&extracted)
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
    next.prepared = std::move(extracted);
    *this = std::move(next);
}
} // namespace bspline2analytic::detail::paired
