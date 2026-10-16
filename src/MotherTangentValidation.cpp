#include "MotherDistanceInternal.h"
#include <functional>

namespace bspline2analytic::detail::mother
{
namespace
{
BPoly Partial(const BPoly &p, BoundaryId id, bool alongU)
{
    const unsigned edge = static_cast<unsigned>(id);
    const bool fixedU = edge < 2, high = (edge % 2) != 0;
    const unsigned fixedDegree = fixedU ? p.u : p.v, varyingDegree = fixedU ? p.v : p.u;
    auto get = [&](unsigned fixed, unsigned varying) {
        return fixedU ? p.at(fixed, varying) : p.at(varying, fixed);
    };
    BPoly row(varyingDegree, 0);
    for (unsigned j = 0; j <= varyingDegree; ++j)
        if (alongU == fixedU)
            row.c[j] = fixedDegree ?
                (high ? get(fixedDegree, j) - get(fixedDegree - 1, j)
                      : get(1, j) - get(0, j)) * I(fixedDegree) : I(0);
        else
            row.c[j] = get(high ? fixedDegree : 0, j);
    return alongU == fixedU ? row : Derivative(row, true);
}
I Value(const PairedParameterMap &m, double t)
{
    auto it = std::lower_bound(m.source.begin(), m.source.end(), t);
    if (it != m.source.end() && *it == t)
        return I(m.target[it - m.source.begin()]);
    if (it == m.source.begin() || it == m.source.end())
        throw Failure(Status::VerificationFailed, "Joint tangent map coverage failure");
    const size_t i = static_cast<size_t>(it - m.source.begin() - 1);
    return I(m.target[i]) + (I(t) - I(m.source[i])) /
        (I(m.source[i + 1]) - I(m.source[i])) * (I(m.target[i + 1]) - I(m.target[i]));
}
BV Crop(const BV &p, double begin, double end, I a, I b)
{
    I lo = a.lo == begin && a.hi == begin ? I(0) : (a - I(begin)) / (I(end) - I(begin));
    I hi = b.lo == end && b.hi == end ? I(1) : (b - I(begin)) / (I(end) - I(begin));
    return paired::Restrict(p, true, lo, hi);
}
struct AngleData
{
    BPoly dot;
    BV cross;
};
std::pair<AngleData, AngleData> SplitAngle(const AngleData &p)
{
    auto d = Split(p.dot, true), x = Split(p.cross.x, true),
         y = Split(p.cross.y, true), z = Split(p.cross.z, true);
    return {{std::move(d.first), {std::move(x.first), std::move(y.first), std::move(z.first)}},
            {std::move(d.second), {std::move(x.second), std::move(y.second), std::move(z.second)}}};
}
I LengthSquared(IV v)
{
    return Sqr(v.x) + Sqr(v.y) + Sqr(v.z);
}
// For positive dot, sin^2(angle)=c^2/(dot^2+c^2). This function is
// increasing in c^2 and decreasing in dot^2; use the corresponding endpoints.
double SineSquaredUpper(const AngleData &p)
{
    I d = p.dot.range(), c = LengthSquared(Range(p.cross));
    if (!(d.lo > 0) || !d.finite() || !c.finite())
        return INFINITY;
    I q = I(c.hi) / (Sqr(I(d.lo)) + I(c.hi));
    return q.finite() ? std::min(1., q.hi) : INFINITY;
}
double Degrees(double sineSquared)
{
    if (!(sineSquared >= 0 && sineSquared < 1))
        return INFINITY;
    // asin(sqrt(q)), expressed with existing outward-rounded operations.
    return (Atan(Sqrt(I(sineSquared) / (I(1) - I(sineSquared)))) * I(180) / IPi()).hi;
}
bool Witness(const std::array<AngleData, 2> &p, I limit)
{
    I sum(0);
    for (const auto &a : p)
    {
        I d = Evaluate(a.dot, .5, 0).range();
        if (d.hi < 0)
            return true; // Oriented parameter directions disagree.
        if (!(d.lo > 0))
            return false; // Zero/uncertain derivatives cannot certify a violation.
        I c = LengthSquared(Range(Evaluate(a.cross, .5, 0)));
        I q = c / (Sqr(d) + c);
        if (!q.finite())
            return false;
        sum = sum + q;
    }
    return sum.lo >= limit.hi;
}
void Coverage(const std::vector<TangentSegment> &segments, double lo, double hi)
{
    double end = lo;
    for (const auto &s : segments)
    {
        if (s.begin != end || !(s.begin < s.end))
            throw Failure(Status::VerificationFailed, "Joint tangent partition gap");
        end = s.end;
    }
    if (segments.empty() || end != hi)
        throw Failure(Status::VerificationFailed, "Joint tangent boundary coverage incomplete");
}
} // namespace

std::array<BV, 2> BoundaryPartials(const BV &p, BoundaryId id)
{
    return {BV{Partial(p.x, id, true), Partial(p.y, id, true), Partial(p.z, id, true)},
            BV{Partial(p.x, id, false), Partial(p.y, id, false), Partial(p.z, id, false)}};
}
BoundaryTangents PrepareTargetTangents(const paired::Chart &chart, paired::Products &products)
{
    BoundaryTangents out;
    // BuildChart supplies orthogonal canonical parameter directions for the five
    // supported mothers. Separate monotone axis maps, swaps and reversals keep
    // them orthogonal. Do not infer orthogonality from a floating-point dot test.
    for (const auto &cell : chart.cells)
        for (unsigned edge = 0; edge < 4; ++edge)
        {
            const auto &s = cell.surface;
            bool outer = edge == 0 ? s.u0 == 0 : edge == 1 ? s.u1 == 1 :
                         edge == 2 ? s.v0 == 0 : s.v1 == 1;
            if (!outer)
                continue;
            auto id = static_cast<BoundaryId>(edge);
            BV h = Boundary(s.p, id);
            BPoly w = Boundary(BV{s.w, s.w, s.w}, id).x;
            if (!(w.range().lo > 0))
                throw Failure(Status::VerificationFailed, "Nonpositive joint tangent weight");
            auto derivatives = BoundaryPartials(s.p, id);
            for (unsigned axis = 0; axis < 2; ++axis)
                derivatives[axis] = products.Subtract(products.Mul(derivatives[axis], w),
                    products.Mul(h, Partial(s.w, id, axis == 0)));
            // H_axis*W-H*W_axis has the direction of Q_axis. The omitted W^2
            // and chart span are positive; no normal or co-normal is needed.
            out[edge].push_back({edge < 2 ? s.v0 : s.u0, edge < 2 ? s.v1 : s.u1,
                                 std::move(derivatives)});
        }
    for (auto &edge : out)
    {
        std::sort(edge.begin(), edge.end(), [](const auto &a, const auto &b) {
            return a.begin < b.begin;
        });
        Coverage(edge, 0, 1);
    }
    return out;
}

bool ValidateJointTangents(const BoundaryTangents &source, const BoundaryTangents &target,
                          const MotherDistanceValidationOptions &options,
                          PairedValidationReport &aux, paired::Products &products,
                          MotherDistanceValidationReport &out)
{
    out.boundaryTangentU = {};
    out.boundaryTangentV = {};
    out.boundaryJointTangent = {};
    out.boundaryTangentU.state = out.boundaryTangentV.state =
        out.boundaryJointTangent.state = CheckState::Unresolved;
    out.jointTangentSineSquaredUpperBound = INFINITY;
    out.jointTangentFailure = "JointTangentConditionUnresolved";
    out.message = "Joint tangent condition not certified; no direct normal fallback";
    const I limit = Sqr(Sin(I(options.jointTangentToleranceDeg) * IPi() / I(180)));
    double maxU = 0, maxV = 0, maxJoint = 0;
    const auto &budget = options.boundary.validation;
    auto region = [&](unsigned edge, double a, double b) {
        const double fixed = edge < 2 ? (edge ? aux.uMap.source.back() : aux.uMap.source.front()) :
            (edge == 3 ? aux.vMap.source.back() : aux.vMap.source.front());
        out.boundaryJointTangent.boundary = static_cast<BoundaryId>(edge);
        out.boundaryJointTangent.parameterRegion = edge < 2 ?
            std::array<double, 4>{fixed, fixed, a, b} : std::array<double, 4>{a, b, fixed, fixed};
    };
    std::function<bool(const std::array<AngleData, 2> &, unsigned, unsigned, double, double)> visit;
    visit = [&](const std::array<AngleData, 2> &p, unsigned depth, unsigned edge, double a, double b) {
        region(edge, a, b);
        if (++aux.certificate.piecesVisited > budget.maxValidationPieces)
        {
            out.jointTangentFailure = "JointTangentBudgetExceeded";
            out.message = "Joint tangent condition not certified: shared boundary node budget "
                          "exhausted; no direct normal fallback";
            throw paired::Budget(out.message);
        }
        ++aux.boundaryNodes;
        double u, v, joint;
        {
            Measure time(aux.edgeBoundsMs);
            u = SineSquaredUpper(p[0]);
            v = SineSquaredUpper(p[1]);
            joint = (I(u) + I(v)).hi;
            if (joint < limit.lo && Degrees(joint) < options.jointTangentToleranceDeg)
            {
                maxU = std::max(maxU, u);
                maxV = std::max(maxV, v);
                maxJoint = std::max(maxJoint, joint);
                return true;
            }
            if (Witness(p, limit))
            {
                out.boundaryJointTangent.state = CheckState::Exceeded;
                out.boundaryJointTangent.hasWitness = true;
                out.jointTangentFailure = "JointTangentConditionExceeded";
                out.message = "Joint tangent condition exceeded at a certified witness: "
                    "sin^2(angle(Su,Qu))+sin^2(angle(Sv,Qv)) must be strictly below "
                    "sin^2(jointTangentToleranceDeg), with positive direction dot products. "
                    "This does not prove the actual normal angle exceeds tolerance; "
                    "no direct normal fallback";
                return false;
            }
        }
        const double middle = a + (b - a) * .5;
        if (depth >= budget.maxSubdivisionDepth || !(a < middle && middle < b))
        {
            out.message = "Joint tangent condition not certified within subdivision/parameter "
                "resolution limit; the joint upper bound or positive nonzero tangent directions "
                "remain unresolved. This is not a proven normal-angle violation; "
                "no direct normal fallback";
            return false;
        }
        std::pair<AngleData, AngleData> uSplit, vSplit;
        {
            Measure time(aux.edgeSubdivisionMs);
            uSplit = SplitAngle(p[0]);
            vSplit = SplitAngle(p[1]);
        }
        return visit({uSplit.first, vSplit.first}, depth + 1, edge, a, middle) &&
               visit({uSplit.second, vSplit.second}, depth + 1, edge, middle, b);
    };
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        const auto &map = edge < 2 ? aux.vMap : aux.uMap;
        const auto &src = source[edge], &dst = target[edge];
        Coverage(src, map.source.front(), map.source.back());
        Coverage(dst, map.target.front(), map.target.back());
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
            const double a = breaks[k], b = breaks[k + 1];
            std::array<BV, 2> s, t;
            {
                Measure time(aux.edgeComposeMs);
                while (si + 1 < src.size() && src[si].end <= a)
                    ++si;
                I mid = Value(map, a + (b - a) * .5);
                while (ti + 1 < dst.size() && dst[ti].end <= mid.lo)
                    ++ti;
                if (!(src[si].begin <= a && b <= src[si].end &&
                      dst[ti].begin < mid.lo && mid.hi < dst[ti].end))
                    throw Failure(Status::VerificationFailed, "Joint tangent segment coverage failure");
                const I ma = Value(map, a), mb = Value(map, b);
                for (unsigned axis = 0; axis < 2; ++axis)
                {
                    s[axis] = Crop(src[si].directions[axis], src[si].begin, src[si].end, I(a), I(b));
                    t[axis] = Crop(dst[ti].directions[axis], dst[ti].begin, dst[ti].end, ma, mb);
                }
            }
            std::array<AngleData, 2> p;
            {
                Measure time(aux.edgeAlgebraMs);
                for (unsigned axis = 0; axis < 2; ++axis)
                    p[axis] = {products.Dot(s[axis], t[axis]), products.Cross(s[axis], t[axis])};
            }
            ++aux.boundarySegments;
            if (!visit(p, 0, edge, a, b))
                return false;
        }
    }
    out.boundaryTangentU.state = out.boundaryTangentV.state =
        out.boundaryJointTangent.state = CheckState::Passed;
    out.boundaryTangentU.upperBound = Degrees(maxU);
    out.boundaryTangentV.upperBound = Degrees(maxV);
    out.boundaryJointTangent.upperBound = Degrees(maxJoint);
    out.jointTangentSineSquaredUpperBound = maxJoint;
    out.jointTangentFailure.clear();
    out.message.clear();
    return true;
}
} // namespace bspline2analytic::detail::mother
