/**
 * @file GeneralBoundaryEnvelope.cpp
 * @brief 一般兜底边界误差的双精度保守包络过滤。
 * 用系数中值和统一误差半径降低区间卷积开销；无法形成有效包络时由完整区间路径继续认证。
 */
#include "GeneralRationalBoundaryInternal.h"

namespace bspline2analytic::detail::rational_boundary
{
namespace
{
// 使用与 MotherSourcePreparation 相同的非负误差算术，但实现独立保留：一般验证器还支持有理源面，缓存合同不同。
double Add(double a, double b)
{
    return a == 0 ? b : b == 0 ? a : std::nextafter(a + b, INFINITY);
}
double ErrorMul(double a, double b)
{
    return a == 0 || b == 0 ? 0 : std::nextafter(a * b, INFINITY);
}
double Gamma(unsigned n)
{
    double nu = n * 0x1p-53;
    return std::nextafter(nu / std::nextafter(1 - nu, 0.), INFINITY);
}
double Round(double magnitude)
{
    return Add(ErrorMul(Gamma(1), magnitude), std::numeric_limits<double>::denorm_min());
}
I Choose(unsigned n, unsigned k)
{
    I out(1);
    k = std::min(k, n - k);
    for (unsigned i = 1; i <= k; ++i)
        out = out * I(n - k + i) / I(i);
    return out;
}
struct Unavailable
{
};
void Check(double x)
{
    if (!std::isfinite(x) || x < 0)
        throw Unavailable{};
}
} // namespace
const CurveEnvelopes::Plan &CurveEnvelopes::Weights(unsigned p, unsigned q)
{
    auto key = std::make_pair(p, q);
    auto found = plans.find(key);
    if (found != plans.end())
        return found->second;
    Plan plan;
    plan.weights.resize((p + 1) * (q + 1));
    std::vector<double> norm(p + q + 1), error(p + q + 1);
    for (unsigned i = 0; i <= p; ++i)
        for (unsigned j = 0; j <= q; ++j)
        {
            I f = (!p || !q) ? I(1) : Choose(p, i) * Choose(q, j) / Choose(p + q, i + j);
            if (!f.finite())
                throw Unavailable{};
            double mid = f.lo * .5 + f.hi * .5;
            plan.weights[i * (q + 1) + j] = mid;
            norm[i + j] = Add(norm[i + j], std::abs(mid));
            error[i + j] = Add(error[i + j], (f - I(mid)).absMax());
        }
    plan.norm = *std::max_element(norm.begin(), norm.end());
    plan.error = *std::max_element(error.begin(), error.end());
    // 每项有两次乘法和一次累加；3m 覆盖包括最终累加在内的舍入路径，不依赖融合乘加 FMA。
    plan.gamma = Gamma(3 * (std::min(p, q) + 1));
    return plans.emplace(key, std::move(plan)).first->second;
}
CurveEnvelopes::Polynomial CurveEnvelopes::Read(const BPoly &p)
{
    if (p.v || p.u > 64 || p.c.size() != p.u + 1)
        throw Unavailable{};
    Polynomial out;
    out.c.reserve(p.c.size());
    for (I x : p.c)
    {
        if (!x.finite())
            throw Unavailable{};
        double mid = x.lo * .5 + x.hi * .5;
        out.c.push_back(mid);
        out.magnitude = std::max(out.magnitude, std::abs(mid));
        out.radius = std::max(out.radius, (x - I(mid)).absMax());
    }
    Check(out.radius);
    return out;
}
CurveEnvelopes::Polynomial CurveEnvelopes::Mul(const Polynomial &a, const Polynomial &b)
{
    unsigned p = unsigned(a.c.size() - 1), q = unsigned(b.c.size() - 1);
    const auto &plan = Weights(p, q);
    Polynomial out;
    out.c.resize(p + q + 1);
    for (unsigned i = 0; i <= p; ++i)
        for (unsigned j = 0; j <= q; ++j)
            out.c[i + j] += (a.c[i] * b.c[j]) * plan.weights[i * (q + 1) + j];
    for (double x : out.c)
    {
        if (!std::isfinite(x))
            throw Unavailable{};
        out.magnitude = std::max(out.magnitude, std::abs(x));
    }
    // 精确输入系数分别被 a、b 的误差半径包围。
    // 每个输出 Bernstein 乘积系数的精确权重均为正且和为 1；这里仍使用缓存的行范数及其不确定性，独立覆盖权重构造的舍入误差。
    double amplitude = ErrorMul(a.magnitude, b.magnitude);
    double input = Add(Add(ErrorMul(a.magnitude, b.radius), ErrorMul(b.magnitude, a.radius)),
                       ErrorMul(a.radius, b.radius));
    // a*b 的下溢误差还会乘以权重 f；附加保护项覆盖这些项、与 f 的乘法以及所有加法，按渐进下溢处理。
    double tiny = ErrorMul(8. * (std::min(p, q) + 1), std::numeric_limits<double>::denorm_min());
    out.radius =
        Add(Add(ErrorMul(Add(plan.norm, plan.error), input), ErrorMul(plan.error, amplitude)),
            Add(ErrorMul(ErrorMul(plan.gamma, plan.norm), amplitude), ErrorMul(tiny, Add(1, plan.norm))));
    Check(out.radius);
    return out;
}
CurveEnvelopes::Polynomial CurveEnvelopes::Subtract(Polynomial a, Polynomial b)
{
    size_t n = std::max(a.c.size(), b.c.size());
    auto elevate = [&](Polynomial &x) {
        if (x.c.size() != n)
        {
            Polynomial one;
            one.c.assign(n - x.c.size() + 1, 1.);
            one.magnitude = 1;
            x = Mul(x, one);
        }
    };
    elevate(a);
    elevate(b);
    Polynomial out;
    out.c.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        out.c[i] = a.c[i] - b.c[i];
        if (!std::isfinite(out.c[i]))
            throw Unavailable{};
        out.magnitude = std::max(out.magnitude, std::abs(out.c[i]));
    }
    out.radius = Add(Add(a.radius, b.radius), Round(out.magnitude));
    Check(out.radius);
    return out;
}
BPoly CurveEnvelopes::Enclose(const Polynomial &p)
{
    BPoly out(unsigned(p.c.size() - 1), 0);
    for (size_t i = 0; i < p.c.size(); ++i)
    {
        out.c[i] = I(p.c[i]) + I(-p.radius, p.radius);
        if (!out.c[i].finite())
            throw Unavailable{};
    }
    return out;
}
std::optional<CurveError> CurveEnvelopes::Build(const BezierCell &s, const BezierCell &t)
{
    try
    {
        auto sw = Read(s.w), tw = Read(t.w);
        auto axis = [&](const BPoly &sp, const BPoly &tp) {
            return Enclose(Subtract(Mul(Read(tp), sw), Mul(Read(sp), tw)));
        };
        return CurveError{{axis(s.p.x, t.p.x), axis(s.p.y, t.p.y), axis(s.p.z, t.p.z)}, Enclose(Mul(sw, tw))};
    }
    catch (const Unavailable &)
    {
        return {};
    }
}
std::optional<BV> CurveEnvelopes::Tangent(const BezierCell &s)
{
    try
    {
        auto w = Read(s.w);
        auto dw = Read(Derivative(s.w, true));
        auto axis = [&](const BPoly &p) {
            auto first = Mul(Read(Derivative(p, true)), w);
            if (!s.w.u)
                return Enclose(first);
            return Enclose(Subtract(std::move(first), Mul(Read(p), dw)));
        };
        return BV{axis(s.p.x), axis(s.p.y), axis(s.p.z)};
    }
    catch (const Unavailable &)
    {
        return {};
    }
}
std::optional<BV> CurveEnvelopes::TangentResidual(const BV &s,const BV &t,const BPoly &scale)
{
    try
    {
        if(!scale.range().finite()||!(scale.range().lo>0))return {};
        auto lambda=Read(scale);
        auto axis=[&](const BPoly &a,const BPoly &b)
        {
            auto target=Read(b);
            Polynomial scaled;
            if(lambda.c.size()==1)
            {
                // 辅助比例为常数时，无需进行 Bernstein 卷积。
                scaled.c.reserve(target.c.size());
                for(double x:target.c)
                {
                    double y=x*lambda.c[0];if(!std::isfinite(y))throw Unavailable{};
                    scaled.c.push_back(y);scaled.magnitude=std::max(scaled.magnitude,std::abs(y));
                }
                scaled.radius=Add(Add(ErrorMul(target.magnitude,lambda.radius),
                                      ErrorMul(lambda.magnitude,target.radius)),
                                  Add(ErrorMul(target.radius,lambda.radius),Round(scaled.magnitude)));
                Check(scaled.radius);
            }
            else scaled=Mul(target,lambda);
            return Enclose(Subtract(Read(a),std::move(scaled)));
        };
        return BV{axis(s.x,t.x),axis(s.y,t.y),axis(s.z,t.z)};
    }
    catch (const Unavailable &)
    {
        return {};
    }
}
I CurveEnvelopes::Midpoint(const BPoly &input)
{
    auto p = Read(input);
    unsigned n = unsigned(p.c.size() - 1);
    for (unsigned k = n; k > 0; --k)
        for (unsigned i = 0; i < k; ++i)
            p.c[i] = p.c[i] * .5 + p.c[i + 1] * .5;
    double radius = Add(p.radius, Add(ErrorMul(Gamma(3 * n), p.magnitude),
                                      ErrorMul(8. * (n + 1), std::numeric_limits<double>::denorm_min())));
    return I(p.c[0]) + I(-radius, radius);
}
bool CurveEnvelopes::PositionExceeded(const CurveError &e, double tolerance)
{
    try
    {
        // 中点见证可以证明同一对应关系下的两个系数位置上界都不能通过；只跳过冗余区间乘积，仍保留原来的细分与映射修正。
        I w = e.denominator.range();
        if (!w.finite() || w.lo <= 0)
            return false;
        I x = Midpoint(e.numerator.x), y = Midpoint(e.numerator.y), z = Midpoint(e.numerator.z);
        return (Sqrt(Sqr(x) + Sqr(y) + Sqr(z)) / I(w.hi)).lo >= tolerance;
    }
    catch (const Unavailable &)
    {
        return false;
    }
}
bool CurveEnvelopes::TangentExceeded(const BV &s, const BV &t, double sineSquaredLimitUpper)
{
    try
    {
        IV a{Midpoint(s.x), Midpoint(s.y), Midpoint(s.z)}, b{Midpoint(t.x), Midpoint(t.y), Midpoint(t.z)};
        auto length2 = [](IV v) { return Sqr(v.x) + Sqr(v.y) + Sqr(v.z); };
        I a2=length2(a),b2=length2(b);
        if(!a2.finite()||!b2.finite()||!(a2.lo>0&&b2.lo>0)||
           !(sineSquaredLimitUpper>0&&sineSquaredLimitUpper<1))
            return false;
        I dot = Dot(a, b);
        if(!dot.finite())return false;
        if (dot.hi <= 0)
            return true;
        if(!(dot.lo>0))return false;
        I denominator=a2*b2;if(!denominator.finite()||!(denominator.lo>0))return false;
        // 容易产生相消的点值恒等式只用于证明失败；点值界无结论时，不能接受整个区间，也不能跳过细分。
        I q=I(1)-Sqr(dot)/denominator;
        return q.finite()&&q.lo>=sineSquaredLimitUpper;
    }
    catch (const Unavailable &)
    {
        return false;
    }
}
} // namespace bspline2analytic::detail::rational_boundary
