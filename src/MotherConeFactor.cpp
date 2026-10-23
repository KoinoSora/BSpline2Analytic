/**
 * @file MotherConeFactor.cpp
 * @brief 用完整重构残差认证近锥顶端点因子，保留闭参数域。
 */
#include "MotherConeFactor.h"

namespace bspline2analytic::detail::mother
{
namespace
{
double VectorUpper(const BV &v)
{
    double upper = 0;
    for (size_t i = 0; i < v.x.c.size(); ++i)
    {
        I a = Sqrt(Sqr(v.x.c[i]) + Sqr(v.y.c[i]) + Sqr(v.z.c[i]));
        if (!a.finite()) return INFINITY;
        upper = std::max(upper, a.hi);
    }
    return upper;
}
// 双精度递推只提议一个存储系数精确确定的 Q；R=P-gQ 的区间重构覆盖全部舍入。
bool Divide(const std::array<Poly, 3> &p, bool along, bool reverse, double small,
            BV &q, BV &r)
{
    unsigned n = along ? p[0].u : p[0].v, m = along ? p[0].v : p[0].u;
    std::array<BPoly *, 3> qs{&q.x, &q.y, &q.z}, rs{&r.x, &r.y, &r.z};
    for (unsigned d = 0; d < 3; ++d)
    {
        *qs[d] = BPoly(p[d].u - unsigned(along), p[d].v - unsigned(!along));
        *rs[d] = BPoly(p[d].u, p[d].v);
        for (unsigned line = 0; line <= m; ++line)
        {
            auto coefficient = [&](unsigned i) {
                unsigned k = reverse ? n - i : i;
                return along ? p[d].at(k, line) : p[d].at(line, k);
            };
            // 源面每方向至多五次；此处不分配逐行动态工作区。
            std::array<double, 5> a{};
            a[n - 1] = coefficient(n); // 大端因子固定为 1。
            for (unsigned i = n - 1; i > 0; --i)
                a[i - 1] = (n * coefficient(i) - (n - i) * small * a[i]) / i;
            for (unsigned i = 0; i < n; ++i)
            {
                if (!std::isfinite(a[i])) return false;
                unsigned k = reverse ? n - 1 - i : i;
                (along ? qs[d]->at(k, line) : qs[d]->at(line, k)) = I(a[i]);
            }
            for (unsigned i = 0; i <= n; ++i)
            {
                I product(0);
                if (i < n) product = I(n - i) / I(n) * I(small) * I(a[i]);
                if (i) product = product + I(i) / I(n) * I(a[i - 1]);
                unsigned k = reverse ? n - i : i;
                I residual = I(coefficient(i)) - product;
                if (!residual.finite()) return false;
                (along ? rs[d]->at(k, line) : rs[d]->at(line, k)) = residual;
            }
        }
    }
    return true;
}
} // namespace

ConeFactorEstimate BoundConeFactor(const ConeFactorCell &p)
{
    ConeFactorEstimate out;
    I g = p.factor.range(), z = p.z.range(), radial = p.radial.range();
    if (!g.finite() || !z.finite() || !radial.finite() || !(radial.hi >= 0)) return out;
    // 非负性来自构造，子片区间舍入不能把这个解析事实丢掉。
    double gmax = std::max(0., g.hi);
    out.perturbation = ErrorAdd(p.enclosure, VectorUpper(p.remainder));
    I rho = Sqrt(I(std::max(0., radial.lo), radial.hi));
    // 锥顶属于母面，因此该界对负高度也成立，不会把负锥瓣当作正锥瓣接受。
    double distance = Sqrt(I(0, radial.hi) + Sqr(z)).hi;
    if (z.lo >= 0)
    {
        I f = p.residual.range(), denominator = rho * p.cosine + z * p.sine;
        if (!f.finite() || !denominator.finite()) return out;
        double numerator = f.absMax();
        distance = std::min(distance, Sqrt(I(numerator)).hi);
        if (denominator.lo > 0)
        {
            out.denominator = denominator.lo;
            distance = std::min(distance, (I(numerator) / I(denominator.lo)).hi);
        }
    }
    if (std::isfinite(distance) && std::isfinite(out.perturbation))
        out.upper = ErrorAdd(out.perturbation, ErrorMul(gmax, distance));
    return out;
}

std::optional<ConeFactorCell> TryConeFactor(const std::array<Poly, 3> &p, double enclosure,
    double angle, double tolerance, paired::Products &products,
    MotherDistanceValidationReport &report, const ExecutionDeadline *deadline)
{
    if (!(tolerance > enclosure) || p[0].u > 5 || p[0].v > 5) return {};
    bool along = false, reverse = false;
    double small = 1;
    // 仅选择高度变化最明显的一个方向。这个启发式只影响是否尝试，不能用来接受。
    for (bool u : {true, false})
    {
        unsigned n = u ? p[2].u : p[2].v, m = u ? p[2].v : p[2].u;
        if (!n) continue;
        double lo = 0, hi = 0;
        for (unsigned j = 0; j <= m; ++j)
        {
            lo += (u ? p[2].at(0, j) : p[2].at(j, 0)) / (m + 1);
            hi += (u ? p[2].at(n, j) : p[2].at(j, n)) / (m + 1);
        }
        double high = std::max(lo, hi);
        if (!(high > 0) || !std::isfinite(high)) continue;
        double ratio = std::max(0., std::min(lo, hi) / high);
        if (std::isfinite(ratio) && ratio < small)
        { small = ratio; along = u; reverse = lo > hi; }
    }
    // 同尺度高度的单元不值得尝试这一后备；不改变原有细分或拒绝语义。
    if (small > .25) return {};
    std::optional<ConeFactorCell> best;
    double bestUpper = INFINITY;
    std::optional<std::pair<I, I>> trig;
    for (unsigned attempt = 0; attempt < (small > 0 ? 2u : 1u); ++attempt)
    {
        if (deadline) deadline->Check("iso_cone_factor_build");
        ++report.coneFactorAttempts;
        double g0 = attempt ? 0 : small;
        BV q, remainder;
        if (!Divide(p, along, reverse, g0, q, remainder)) continue;
        double perturbation = ErrorAdd(enclosure, VectorUpper(remainder));
        if (!(perturbation < tolerance)) { ++report.coneFactorRemainderSkips; continue; }
        if (!trig) trig = std::make_pair(Cos(I(angle)), Sin(I(angle)));
        if (!(trig->first.lo > 0 && trig->second.lo > 0)) return {};
        ConeFactorCell cell;
        cell.cosine = trig->first; cell.sine = trig->second; cell.enclosure = enclosure;
        cell.factor = BPoly(along ? 1 : 0, along ? 0 : 1);
        cell.factor.c = reverse ? std::vector<I>{I(1), I(g0)} : std::vector<I>{I(g0), I(1)};
        cell.remainder = std::move(remainder);
        cell.radial = products.Mul(q.x, q.x) + products.Mul(q.y, q.y);
        if (deadline) deadline->Check("iso_cone_factor_product");
        cell.residual = cell.radial * Sqr(cell.cosine) - products.Mul(q.z, q.z) * Sqr(cell.sine);
        cell.z = std::move(q.z);
        double upper = BoundConeFactor(cell).upper;
        if (upper < bestUpper) { bestUpper = upper; best = std::move(cell); }
        if (bestUpper < tolerance) break;
    }
    if (best) ++report.coneFactorCells;
    return best;
}

std::pair<ConeFactorCell, ConeFactorCell> SplitConeFactor(const ConeFactorCell &p, bool along)
{
    ConeFactorCell a, b;
    auto split = [&](const BPoly &x, BPoly &l, BPoly &h) {
        auto parts = Split(x, along); l = std::move(parts.first); h = std::move(parts.second);
    };
    split(p.residual, a.residual, b.residual); split(p.radial, a.radial, b.radial);
    split(p.z, a.z, b.z); split(p.factor, a.factor, b.factor);
    split(p.remainder.x, a.remainder.x, b.remainder.x);
    split(p.remainder.y, a.remainder.y, b.remainder.y);
    split(p.remainder.z, a.remainder.z, b.remainder.z);
    a.cosine = b.cosine = p.cosine; a.sine = b.sine = p.sine;
    a.enclosure = b.enclosure = p.enclosure;
    return {std::move(a), std::move(b)};
}
} // namespace bspline2analytic::detail::mother
