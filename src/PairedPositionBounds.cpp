/**
 * @file PairedPositionBounds.cpp
 * @brief 配对位置残差的快速与紧致连续上界。
 * 先用系数包围盒筛查，必要时升阶对齐并采用控制向量界；缓存只优化运算量，不改变容差。
 */
#include "PairedValidationInternal.h"

namespace bspline2analytic::detail::paired
{
BV Products::MulShared(const BV &a, const BPoly &b)
{
    if (a.x.u != a.y.u || a.x.u != a.z.u || a.x.v != a.y.v || a.x.v != a.z.v)
        return Mul(a, b);
    if (!b.u && !b.v)
        return a * b.c[0];
    unsigned p = a.x.u, q = a.x.v;
    const auto *wu = p && b.u ? &Weights(p, b.u) : nullptr;
    const auto *wv = q && b.v ? &Weights(q, b.v) : nullptr;
    BV out{BPoly(p + b.u, q + b.v), BPoly(p + b.u, q + b.v), BPoly(p + b.u, q + b.v)};
    for (unsigned i = 0; i <= p; ++i)
        for (unsigned j = 0; j <= q; ++j)
            for (unsigned k = 0; k <= b.u; ++k)
                for (unsigned l = 0; l <= b.v; ++l)
                {
                    I factor = b.at(k, l);
                    if (wu)
                        factor = factor * (*wu)[i * (b.u + 1) + k];
                    if (wv)
                        factor = factor * (*wv)[j * (b.v + 1) + l];
                    size_t from = i * (q + 1) + j, to = (i + k) * (out.x.v + 1) + j + l;
                    out.x.c[to] = out.x.c[to] + a.x.c[from] * factor;
                    out.y.c[to] = out.y.c[to] + a.y.c[from] * factor;
                    out.z.c[to] = out.z.c[to] + a.z.c[from] * factor;
                }
    return out;
}
BV Products::MulWeight(const BV &a, const TargetCell &target)
{
    if (!target.separableWeight)
        return MulShared(a, target.surface.w);
    if (target.weightFactors[0].v || target.weightFactors[1].u)
        throw Failure(Status::VerificationFailed, "Invalid oriented weight factor axes");
    return MulShared(MulShared(a, target.weightFactors[0]), target.weightFactors[1]);
}
BPoly Products::Elevate(const BPoly &p, unsigned u, unsigned v)
{
    if (u < p.u || v < p.v)
        throw Failure(Status::InvalidInput, "Degree elevation cannot reduce degree");
    BPoly current = p;
    for (bool along : {true, false})
    {
        unsigned n = along ? current.u : current.v, target = along ? u : v;
        if (n == target)
            continue;
        unsigned delta = target - n, lines = along ? current.v : current.u;
        BPoly next(along ? target : current.u, along ? current.v : target);
        if (n == 0)
        {
            // 利用单位分解，某方向上为常数的系数可精确复制。
            for (unsigned j = 0; j <= lines; ++j)
                for (unsigned i = 0; i <= target; ++i)
                    (along ? next.at(i, j) : next.at(j, i)) =
                        along ? current.at(0, j) : current.at(j, 0);
        }
        else
        {
            // 与表示常数 1 的次数差多项式相乘等价于直接升阶；按次数缓存向外舍入的二项式系数比。
            const auto &w = Weights(n, delta);
            for (unsigned j = 0; j <= lines; ++j)
                for (unsigned i = 0; i <= n; ++i)
                    for (unsigned k = 0; k <= delta; ++k)
                    {
                        I &out = along ? next.at(i + k, j) : next.at(j, i + k);
                        out = out + (along ? current.at(i, j) : current.at(j, i)) *
                                        w[i * (delta + 1) + k];
                    }
        }
        current = std::move(next);
    }
    return current;
}
BV Products::Subtract(BV a, const BV &b)
{
    for (auto pair :
         {std::make_pair(&a.x, &b.x), std::make_pair(&a.y, &b.y), std::make_pair(&a.z, &b.z)})
    {
        auto &x = *pair.first;
        const auto &y = *pair.second;
        unsigned u = std::max(x.u, y.u), v = std::max(x.v, y.v);
        if (x.u != u || x.v != v)
            x = Elevate(x, u, v);
        BPoly storage;
        const BPoly &aligned = y.u == u && y.v == v ? y : (storage = Elevate(y, u, v));
        for (size_t i = 0; i < x.c.size(); ++i)
            x.c[i] = x.c[i] - aligned.c[i];
    }
    return a;
}
namespace
{
I Length(IV x)
{
    return Sqrt(Sqr(x.x) + Sqr(x.y) + Sqr(x.z));
}
} // namespace
double PositionBoxUpper(const BV &error, const BPoly &weight)
{
    I w = weight.range();
    IV box = Range(error);
    if (!w.finite() || w.lo <= 0 || !box.x.finite() || !box.y.finite() || !box.z.finite())
        return INFINITY;
    // 即使次数不同，各多项式也分别位于自身系数凸包中。
    // 上界未通过只能说明无结论；只有严格通过才可提前返回。
    return (Length(box) / I(w.lo)).hi;
}
double PositionControlUpper(const BV &error, const BPoly &weight, Products &products)
{
    unsigned u = std::max({error.x.u, error.y.u, error.z.u, weight.u}),
             v = std::max({error.x.v, error.y.v, error.z.v, weight.v});
    std::array<BPoly, 4> storage;
    auto aligned = [&](const BPoly &p, unsigned slot) -> const BPoly & {
        return p.u == u && p.v == v ? p : (storage[slot] = products.Elevate(p, u, v));
    };
    const auto &x = aligned(error.x, 0), &y = aligned(error.y, 1), &z = aligned(error.z, 2),
               &w = aligned(weight, 3);
    double upper = 0;
    for (size_t i = 0; i < w.c.size(); ++i)
    {
        if (!w.c[i].finite() || w.c[i].lo <= 0)
            return INFINITY;
        I inverse = I(1) / w.c[i];
        upper = std::max(upper, Length({x.c[i] * inverse, y.c[i] * inverse, z.c[i] * inverse}).hi);
    }
    return upper;
}
} // namespace bspline2analytic::detail::paired
