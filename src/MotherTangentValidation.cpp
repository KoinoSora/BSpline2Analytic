/**
 * @file MotherTangentValidation.cpp
 * @brief 同一边界子区间上的 U/V 参数切向联合认证。
 * 用正比例残差给出两方向正弦平方上界并相加；无结论则同步细分，有见证才报告条件超限。
 */
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
    BV source, target;
};
std::pair<AngleData, AngleData> SplitAngle(const AngleData &p)
{
    auto split = [](const BV &v) {
        auto x = Split(v.x, true), y = Split(v.y, true), z = Split(v.z, true);
        return std::make_pair(BV{std::move(x.first), std::move(y.first), std::move(z.first)},
                              BV{std::move(x.second), std::move(y.second), std::move(z.second)});
    };
    auto s = split(p.source), t = split(p.target);
    return {{std::move(s.first), std::move(t.first)},
            {std::move(s.second), std::move(t.second)}};
}
I LengthSquared(IV v)
{
    return Sqr(v.x) + Sqr(v.y) + Sqr(v.z);
}
double Estimate(const BPoly &p, double t)
{
    // 这里只作启发式估计；拟合比例和参考方向无需精确插值原向量场，后续独立构造完整证书。
    auto midpoint = [](I x) { return x.lo * .5 + x.hi * .5; };
    if (t == 0) return midpoint(p.c.front());
    if (t == 1) return midpoint(p.c.back());
    std::vector<double> values;
    values.reserve(p.c.size());
    for (I c : p.c) values.push_back(midpoint(c));
    for (size_t n = values.size(); n > 1; --n)
        for (size_t i = 0; i + 1 < n; ++i)
            values[i] = (1 - t) * values[i] + t * values[i + 1];
    return values.front();
}
V Estimate(const BV &p, double t)
{
    return {Estimate(p.x, t), Estimate(p.y, t), Estimate(p.z, t)};
}
struct ResidualData
{
    BPoly scale;
    BV residual;
};
// 在端点和中点选择 lambda 的二次 Bernstein 候选；近常数或近线性时降低次数。
// 必须证明 lambda 在整段严格为正；正性无保证时改用正的中点常数。
// 最终以区间运算构造 A-lambda*B，三点拟合残差小本身不构成连续角度证书。
ResidualData BuildResidual(const AngleData &p, paired::Products &products)
{
    // 点值只用于选择辅助正多项式；验收始终包围完整残差，包括全部系数误差。
    double ratio[3];
    for (unsigned i = 0; i < 3; ++i)
    {
        V a = Estimate(p.source, .5 * i);
        V b = Estimate(p.target, .5 * i);
        ratio[i] = Dot(a, b) / Dot(b, b);
        if (!(ratio[i] > 0) || !std::isfinite(ratio[i]))
            return {BPoly(0, 0, I(0)), {}}; // 无法判定，不是角度超限见证。
    }
    double magnitude = std::max({ratio[0], ratio[1], ratio[2]});
    BPoly scale(2, 0);
    scale.c = {I(ratio[0]), I(2 * ratio[1] - .5 * ratio[0] - .5 * ratio[2]), I(ratio[2])};
    if (std::max(std::abs(ratio[0] - ratio[1]), std::abs(ratio[2] - ratio[1])) <=
        64 * std::numeric_limits<double>::epsilon() * magnitude || !(scale.range().lo > 0))
        scale = BPoly(0, 0, I(ratio[1]));
    else if (std::abs(ratio[1] - .5 * ratio[0] - .5 * ratio[2]) <=
             64 * std::numeric_limits<double>::epsilon() * magnitude)
    {
        scale = BPoly(1, 0);
        scale.c = {I(ratio[0]), I(ratio[2])};
    }
    if (!scale.range().finite())
        return {BPoly(0, 0, I(0)), {}};
    return {scale, products.Subtract(p.source, products.Mul(p.target, scale))};
}
// 用固定参考方向 e 的正投影下界得到 ||A||² >= lower(A·e)² / ||e||²。
// 再利用 sin²(theta) <= ||A-lambda*B||² / ||A||² 得到保守上界。
// 比值严格小于 1 才同时保证目标非零且同向；退化或过宽区间返回无穷，交给细分处理。
double SineSquaredUpper(const AngleData &p, const ResidualData &r)
{
    if (!(r.scale.range().lo > 0))
        return INFINITY;
    V e = Estimate(p.source, .5);
    double size = std::max({std::abs(e.x), std::abs(e.y), std::abs(e.z)});
    if (!(size > 0) || !std::isfinite(size))
        return INFINITY;
    e = e * (1 / size); // 固定参考方向，无需精确单位长度。
    I projection = Dot(p.source, Exact(e)).range();
    if (!(projection.lo > 0) || !projection.finite())
        return INFINITY;
    I length2 = Sqr(I(projection.lo)) / LengthSquared(Exact(e));
    if (!(length2.lo > 0) || !length2.finite())
        return INFINITY;
    const BV &v = r.residual;
    double error2 = 0;
    if (v.x.u == v.y.u && v.x.u == v.z.u && v.x.v == v.y.v && v.x.v == v.z.v)
    {
        // 利用范数凸性直接对向量控制系数求界，不新增范数平方多项式；系数已包含向外舍入保护。
        for (size_t i = 0; i < v.x.c.size(); ++i)
        {
            I norm2 = LengthSquared({v.x.c[i], v.y.c[i], v.z.c[i]});
            if (!norm2.finite()) return INFINITY;
            error2 = std::max(error2, norm2.hi);
        }
    }
    else
    {
        I norm2 = LengthSquared(Range(v));
        if (!norm2.finite()) return INFINITY;
        error2 = norm2.hi;
    }
    I q = I(error2) / I(length2.lo);
    // 当 ||A-lambda*B||<||A|| 且 lambda>0 时，还可证明 A.B>0 及 B 非零。
    // 残差上界未通过不等于真实夹角超限。
    return q.finite() && q.hi < 1 ? std::max(0., q.hi) : INFINITY;
}
double Degrees(double sineSquared)
{
    if (!(sineSquared >= 0 && sineSquared < 1))
        return INFINITY;
    // 用已有向外舍入运算表示 asin(sqrt(q))，保持报告角度的保守性。
    return (Atan(Sqrt(I(sineSquared) / (I(1) - I(sineSquared)))) * I(180) / IPi()).hi;
}
bool Witness(const std::array<AngleData, 2> &p, I limit)
{
    I sum(0);
    for (const auto &a : p)
    {
        IV s = Range(Evaluate(a.source, .5, 0)), t = Range(Evaluate(a.target, .5, 0));
        I d = Dot(s, t);
        if (d.hi < 0)
            return true; // 有向参数切向不一致。
        if (!(d.lo > 0))
            return false; // 零导数或不确定导数不足以证明条件超限。
        I denominator = LengthSquared(s) * LengthSquared(t);
        if (!(denominator.lo > 0) || !denominator.finite())
            return false;
        // 这里只使用点值处向外舍入的点积/范数恒等式；可证明失败，不能据此认证整个区间通过。
        I q = I(1) - Sqr(d) / denominator;
        if (!q.finite())
            return false;
        sum = sum + I(std::max(0., q.lo), std::max(0., q.hi));
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
    // BuildChart 为五类支持母面构造解析正交的规范参数方向；分离的单调轴映射、U/V 交换和反向保持该正交性，不用浮点点积近似测试来推断。
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
            // H_axis*W-H*W_axis 与 Q_axis 同向；省略的 W² 与参数图跨度均为正，无需构造法向或余法向。
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
                          MotherDistanceValidationReport &out, const ExecutionDeadline *deadline)
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
    // 两方向始终在同一个 [a,b] 上联合验收和细分，不能把不同区间的独立最大角相加。
    // 上界未通过先检查独立超限见证；没有见证时只能细分或报告无法认证。
    visit = [&](const std::array<AngleData, 2> &p, unsigned depth, unsigned edge, double a, double b) {
        if (deadline) deadline->Check("iso_joint_tangent");
        region(edge, a, b);
        if (++aux.certificate.piecesVisited > budget.maxValidationPieces)
        {
            out.jointTangentFailure = "JointTangentBudgetExceeded";
            out.message = "Joint tangent condition not certified: shared boundary node budget "
                          "exhausted; no direct normal fallback";
            throw paired::Budget(out.message);
        }
        ++aux.boundaryNodes;
        std::array<ResidualData, 2> residuals;
        {
            Measure time(aux.edgeAlgebraMs);
            for (unsigned axis = 0; axis < 2; ++axis)
                residuals[axis] = BuildResidual(p[axis], products);
            out.tangentResidualFits += 2;
        }
        double u, v, joint;
        {
            Measure time(aux.edgeBoundsMs);
            u = SineSquaredUpper(p[0], residuals[0]);
            v = SineSquaredUpper(p[1], residuals[1]);
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
                "resolution limit; the positive-scale residual estimate or nonzero tangent directions "
                "remain unresolved. This is not a proven normal-angle violation; "
                "no direct normal fallback";
            return false;
        }
        std::pair<AngleData, AngleData> uSplit, vSplit;
        {
            Measure time(aux.edgeSubdivisionMs);
            uSplit = SplitAngle(p[0]);
            vSplit = SplitAngle(p[1]);
            ++out.tangentResidualSplits;
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
            std::array<AngleData, 2> p{{{std::move(s[0]), std::move(t[0])},
                                       {std::move(s[1]), std::move(t[1])}}};
            ++aux.boundarySegments;
            if (!visit(p, 0, edge, a, b))
                return false;
        }
    }
    // 只有四边所有公共子区间均通过，才能发布全边界上界；部分遍历的最大值不够。
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
