/**
 * @file MotherDistanceBoundary.cpp
 * @brief 母面验证器的四边区间对应位置检查。
 * 源段与目标有理段按共同参数映射限制和细分；与 U/V 联合切向共享验收流程但分别报告误差。
 */
#include "MotherDistanceInternal.h"
#include <functional>

namespace bspline2analytic::detail::mother
{
namespace
{
struct Edge { double begin, end; BV p; BPoly w; };
using Edges = std::array<std::vector<Edge>, 4>;
void Add(Edges &out, const BezierCell &cell, const std::array<double, 4> &domain)
{
    for (unsigned e = 0; e < 4; ++e)
    {
        bool outer = e == 0 ? cell.u0 == domain[0] : e == 1 ? cell.u1 == domain[1] :
                     e == 2 ? cell.v0 == domain[2] : cell.v1 == domain[3];
        if (outer)
        {
            auto id = static_cast<BoundaryId>(e);
            out[e].push_back({e < 2 ? cell.v0 : cell.u0, e < 2 ? cell.v1 : cell.u1,
                Boundary(cell.p, id), Boundary(BV{cell.w, cell.w, cell.w}, id).x});
        }
    }
}
void Sort(Edges &edges, const PairedParameterMap &u, const PairedParameterMap &v)
{
    for (unsigned e = 0; e < 4; ++e)
    {
        auto &edge = edges[e];
        const auto &domain = e < 2 ? v.source : u.source;
        std::sort(edge.begin(), edge.end(), [](const Edge &a, const Edge &b) { return a.begin < b.begin; });
        double end = domain.front();
        for (const auto &s : edge)
        {
            if (s.begin != end || !(s.end > s.begin))
                throw Failure(Status::VerificationFailed, "Boundary position partition gap");
            end = s.end;
        }
        if (end != domain.back())
            throw Failure(Status::VerificationFailed, "Boundary position coverage incomplete");
    }
}
I Value(const PairedParameterMap &m, double t)
{
    auto it = std::lower_bound(m.source.begin(), m.source.end(), t);
    if (it != m.source.end() && *it == t) return I(m.target[it - m.source.begin()]);
    if (it == m.source.begin() || it == m.source.end())
        throw Failure(Status::VerificationFailed, "Boundary position map coverage failure");
    size_t k = it - m.source.begin() - 1;
    return I(m.target[k]) + (I(t) - I(m.source[k])) / (I(m.source[k + 1]) - I(m.source[k])) *
        (I(m.target[k + 1]) - I(m.target[k]));
}
Edge Crop(const Edge &edge, I lo, I hi)
{
    I a = lo.lo == edge.begin && lo.hi == edge.begin ? I(0) :
        (lo - I(edge.begin)) / (I(edge.end) - I(edge.begin));
    I b = hi.lo == edge.end && hi.hi == edge.end ? I(1) :
        (hi - I(edge.begin)) / (I(edge.end) - I(edge.begin));
    return {edge.begin, edge.end, paired::Restrict(edge.p, true, a, b),
        paired::Restrict(edge.w, true, a, b)};
}
} // namespace
bool ValidateEdgePositions(const paired::PreparedSource &source, const paired::Chart &chart,
                           const PairedValidationOptions &options, PairedValidationReport &aux,
                           paired::Products &products, MotherDistanceValidationReport &out, const ExecutionDeadline *deadline)
{
    Measure time(out.boundaryPositionMs);
    out.boundaryPosition = {};
    Edges src, dst;
    std::array<double, 4> domain{aux.uMap.source.front(), aux.uMap.source.back(),
        aux.vMap.source.front(), aux.vMap.source.back()};
    for (const auto &cell : source.cells) Add(src, cell, domain);
    for (const auto &cell : chart.cells) Add(dst, cell.surface, {0, 1, 0, 1});
    Sort(src, aux.uMap, aux.vMap);
    PairedParameterMap unit{{0, 1}, {0, 1}};
    Sort(dst, unit, unit);
    double maximum = 0;
    std::function<bool(const BV &, const BPoly &, unsigned)> visit;
    visit = [&](const BV &error, const BPoly &weight, unsigned depth) {
        if (deadline) deadline->Check("iso_boundary_position");
        if (++aux.certificate.piecesVisited > options.validation.maxValidationPieces)
            throw paired::Budget("Boundary node budget exhausted");
        ++out.boundaryPositionNodes;
        double upper = paired::PositionBoxUpper(error, weight);
        if (!(upper < options.validation.distanceTolerance))
            upper = paired::PositionControlUpper(error, weight, products);
        if (upper < options.validation.distanceTolerance)
        {
            maximum = std::max(maximum, upper);
            return true;
        }
        IV e = Range(Evaluate(error, .5, 0));
        I w = Evaluate(weight, .5, 0).range();
        I distance = Sqrt(Sqr(e.x) + Sqr(e.y) + Sqr(e.z)) / w;
        if (w.lo > 0 && distance.lo >= options.validation.distanceTolerance)
        {
            out.boundaryPosition.state = CheckState::Exceeded;
            out.boundaryPosition.hasWitness = true;
            return false;
        }
        if (depth >= options.validation.maxSubdivisionDepth) return false;
        auto x = Split(error.x, true), y = Split(error.y, true), z = Split(error.z, true);
        auto wsplit = Split(weight, true);
        return visit(BV{x.first, y.first, z.first}, wsplit.first, depth + 1) &&
            visit(BV{x.second, y.second, z.second}, wsplit.second, depth + 1);
    };
    for (unsigned e = 0; e < 4; ++e)
    {
        const auto &map = e < 2 ? aux.vMap : aux.uMap;
        std::vector<double> breaks = map.source;
        for (const auto &s : src[e]) { breaks.push_back(s.begin); breaks.push_back(s.end); }
        std::sort(breaks.begin(), breaks.end());
        breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
        size_t si = 0, ti = 0;
        for (size_t k = 0; k + 1 < breaks.size(); ++k)
        {
            double a = breaks[k], b = breaks[k + 1];
            while (si + 1 < src[e].size() && src[e][si].end <= a) ++si;
            I mid = Value(map, a + (b - a) * .5);
            while (ti + 1 < dst[e].size() && dst[e][ti].end <= mid.lo) ++ti;
            if (!(src[e][si].begin <= a && b <= src[e][si].end &&
                  dst[e][ti].begin < mid.lo && mid.hi < dst[e][ti].end))
                throw Failure(Status::VerificationFailed, "Boundary position segment coverage failure");
            Edge s = Crop(src[e][si], I(a), I(b));
            Edge t = Crop(dst[e][ti], Value(map, a), Value(map, b));
            BV error = products.Subtract(products.Mul(s.p, t.w), t.p);
            if (!visit(error, t.w, 0))
            {
                out.boundaryPosition.boundary = static_cast<BoundaryId>(e);
                out.boundaryPosition.parameterRegion = e < 2 ?
                    std::array<double, 4>{domain[e], domain[e], a, b} :
                    std::array<double, 4>{a, b, domain[e], domain[e]};
                if (out.boundaryPosition.state != CheckState::Exceeded)
                    out.boundaryPosition.state = CheckState::Unresolved;
                return false;
            }
        }
    }
    out.boundaryPosition.state = CheckState::Passed;
    out.boundaryPosition.upperBound = maximum;
    return true;
}
} // namespace bspline2analytic::detail::mother
