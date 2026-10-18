/**
 * @file PairedBoundaryValidation.cpp
 * @brief 配对验证器的连续边界位置、切向和法向认证。
 * 在共同参数对应下组合源/目标导数并求区间界；该法向合同与母面验证器的联合切向合同不同。
 */
#include "PairedValidationInternal.h"
#include <chrono>
#include <functional>

namespace bspline2analytic::detail::paired
{
namespace
{
struct Measure
{
    double &value;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    explicit Measure(double &v) : value(v) {}
    ~Measure()
    {
        value += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                     .count();
    }
};
I Length(IV a)
{
    return Sqrt(Sqr(a.x) + Sqr(a.y) + Sqr(a.z));
}
I Value(const PairedParameterMap &m, double t)
{
    auto it = std::lower_bound(m.source.begin(), m.source.end(), t);
    if (it != m.source.end() && *it == t)
        return I(m.target[it - m.source.begin()]);
    if (it == m.source.begin() || it == m.source.end())
        throw Failure(Status::VerificationFailed, "Boundary map coverage failure");
    size_t i = static_cast<size_t>(it - m.source.begin() - 1);
    return I(m.target[i]) + (I(t) - I(m.source[i])) / (I(m.source[i + 1]) - I(m.source[i])) *
                                (I(m.target[i + 1]) - I(m.target[i]));
}
BV CropEdge(const BV &p, double begin, double end, I a, I b)
{
    I x = a.lo == begin && a.hi == begin ? I(0) : (a - I(begin)) / (I(end) - I(begin));
    I y = b.lo == end && b.hi == end ? I(1) : (b - I(begin)) / (I(end) - I(begin));
    return Restrict(p, true, x, y);
}
struct AngleData
{
    BPoly dot;
    BV cross;
};
struct AngleResult
{
    double ratio;
    bool exceeded;
};
AngleResult Bound(const AngleData &p, I limit)
{
    I dot = p.dot.range();
    double ratio = dot.lo > 0 ? (Length(Range(p.cross)) / I(dot.lo)).hi : INFINITY;
    if (ratio < limit.lo)
        return {ratio, false};
    I d = Evaluate(p.dot, .5, 0).range(), c = Length(Range(Evaluate(p.cross, .5, 0)));
    return {ratio, d.hi < 0 || (d.lo > 0 && (c / d).lo >= limit.hi)};
}
std::pair<AngleData, AngleData> SplitAngle(const AngleData &p)
{
    auto d = Split(p.dot, true), x = Split(p.cross.x, true), y = Split(p.cross.y, true),
         z = Split(p.cross.z, true);
    return {{std::move(d.first), {std::move(x.first), std::move(y.first), std::move(z.first)}},
            {std::move(d.second), {std::move(x.second), std::move(y.second), std::move(z.second)}}};
}
} // namespace

Boundaries PrepareTargetBoundaries(const Chart &chart, Products &products)
{
    Boundaries out;
    for (const auto &cell : chart.cells)
        for (unsigned edge = 0; edge < 4; ++edge)
        {
            const auto &s = cell.surface;
            bool outer = edge == 0   ? s.u0 == 0
                         : edge == 1 ? s.u1 == 1
                         : edge == 2 ? s.v0 == 0
                                     : s.v1 == 1;
            if (!outer)
                continue;
            auto id = static_cast<BoundaryId>(edge);
            BV p = Boundary(s.p, id);
            BPoly w = Boundary(BV{s.w, s.w, s.w}, id).x;
            if (w.range().lo <= 0)
                throw Failure(Status::VerificationFailed, "Nonpositive boundary weight");
            BV tangent =
                products.Mul(Derivative(p, true), w) - products.Mul(p, Derivative(w, true));
            out[edge].push_back({edge < 2 ? s.v0 : s.u0, edge < 2 ? s.v1 : s.u1, std::move(tangent),
                                 chart.hasNormals ? Boundary(cell.normal, id) : BV{}});
        }
    for (auto &edge : out)
    {
        std::sort(edge.begin(), edge.end(),
                  [](const auto &a, const auto &b) { return a.begin < b.begin; });
        double end = 0;
        for (const auto &segment : edge)
        {
            if (segment.begin != end || !(segment.begin < segment.end))
                throw Failure(Status::VerificationFailed, "Target boundary partition gap");
            end = segment.end;
        }
        if (end != 1)
            throw Failure(Status::VerificationFailed, "Incomplete target boundary");
    }
    return out;
}

BoundaryResult ValidateBoundaries(const Boundaries &source, const Boundaries &target,
                                  const PairedParameterMap &uMap, const PairedParameterMap &vMap,
                                  const PairedValidationOptions &options,
                                  PairedValidationReport &report, Products &products,
                                  bool checkNormals)
{
    BoundaryResult result;
    I t = I(options.validation.boundaryTangentToleranceDeg) * IPi() / I(180);
    I tLimit = Sin(t) / Cos(t), nLimit;
    if (checkNormals)
    {
        I n = I(options.validation.boundaryNormalToleranceDeg) * IPi() / I(180);
        nLimit = Sin(n) / Cos(n);
    }
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        const auto &map = edge < 2 ? vMap : uMap;
        const auto &src = source[edge], &dst = target[edge];
        if (src.empty() || dst.empty() || src.front().begin != map.source.front() ||
            src.back().end != map.source.back())
            throw Failure(Status::VerificationFailed, "Source boundary coverage failure");
        std::vector<double> breaks = map.source;
        for (const auto &s : src)
        {
            breaks.push_back(s.begin);
            breaks.push_back(s.end);
        }
        std::sort(breaks.begin(), breaks.end());
        breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
        size_t si = 0, ti = 0;
        for (size_t k = 0; k + 1 < breaks.size(); ++k)
        {
            double a = breaks[k], b = breaks[k + 1];
            BV st, sn, tt, tn;
            {
                Measure time(report.edgeComposeMs);
                while (si + 1 < src.size() && src[si].end <= a)
                    ++si;
                I mid = Value(map, a + (b - a) * .5);
                while (ti + 1 < dst.size() && dst[ti].end <= mid.lo)
                    ++ti;
                if (!(src[si].begin <= a && b <= src[si].end && dst[ti].begin < mid.lo &&
                      mid.hi < dst[ti].end))
                    throw Failure(Status::VerificationFailed,
                                  "Boundary partition coverage failure");
                st = CropEdge(src[si].tangent, src[si].begin, src[si].end, I(a), I(b));
                if (checkNormals)
                    sn = CropEdge(src[si].normal, src[si].begin, src[si].end, I(a), I(b));
                I ma = Value(map, a), mb = Value(map, b);
                tt = CropEdge(dst[ti].tangent, dst[ti].begin, dst[ti].end, ma, mb);
                if (checkNormals)
                    tn = CropEdge(dst[ti].normal, dst[ti].begin, dst[ti].end, ma, mb);
            }
            AngleData tangent, normal;
            {
                Measure time(report.edgeAlgebraMs);
                tangent = {products.Dot(st, tt), products.Cross(st, tt)};
                if (checkNormals)
                    normal = {products.Dot(sn, tn), products.Cross(sn, tn)};
            }
            ++report.boundarySegments;
            std::function<bool(const AngleData &, const AngleData &, unsigned)> certify;
            certify = [&](const AngleData &ta, const AngleData &no, unsigned depth) {
                if (++report.certificate.piecesVisited > options.validation.maxValidationPieces)
                    throw Budget("Validation node budget exhausted");
                ++report.boundaryNodes;
                AngleResult tr, nr{0, false};
                {
                    Measure time(report.edgeBoundsMs);
                    tr = Bound(ta, tLimit);
                    if (checkNormals)
                        nr = Bound(no, nLimit);
                }
                result.tangentExceeded |= tr.exceeded;
                result.normalExceeded |= nr.exceeded;
                if (tr.ratio < tLimit.lo && (!checkNormals || nr.ratio < nLimit.lo))
                {
                    result.tangentRatio = std::max(result.tangentRatio, tr.ratio);
                    result.normalRatio = std::max(result.normalRatio, nr.ratio);
                    return true;
                }
                if (tr.exceeded || nr.exceeded || depth >= options.validation.maxSubdivisionDepth)
                    return false;
                std::pair<AngleData, AngleData> ts, ns;
                {
                    Measure time(report.edgeSubdivisionMs);
                    ts = SplitAngle(ta);
                    if (checkNormals)
                        ns = SplitAngle(no);
                }
                return certify(ts.first, ns.first, depth + 1) &&
                       certify(ts.second, ns.second, depth + 1);
            };
            if (!certify(tangent, normal, 0))
                return result;
        }
    }
    result.passed = true;
    return result;
}
} // namespace bspline2analytic::detail::paired
