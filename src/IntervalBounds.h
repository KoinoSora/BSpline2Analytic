/**
 * @file IntervalBounds.h
 * @brief 区间标量、向量及齐次 Bézier 单元的数据结构。
 * 系数按 U 外层、V 内层存储；有理面为 p/w，所有接受判据须使用保守区间端点。
 */
#pragma once
#include "Numerics.h"
namespace bspline2analytic::detail
{
// 基础运算采用向外舍入；本模块须使用 /fp:strict 编译，保留浮点求界语义。
struct I
{
    double lo = 0, hi = 0;
    I() = default;
    I(double x) : lo(std::isnan(x) ? -INFINITY : x), hi(std::isnan(x) ? INFINITY : x) {}
    I(double a, double b) : lo(a <= b ? a : -INFINITY), hi(a <= b ? b : INFINITY) {}
    static I Whole()
    {
        return {-INFINITY, INFINITY};
    }
    bool finite() const
    {
        return std::isfinite(lo) && std::isfinite(hi) && lo <= hi;
    }
    double absMax() const
    {
        return std::max(std::abs(lo), std::abs(hi));
    }
};
I operator+(I, I);
I operator-(I, I);
I operator-(I);
I operator*(I, I);
I operator/(I, I);
I Sqr(I);
I Sqrt(I);
I Sin(I);
I Cos(I);
I Atan(I);
I Acos(I);
I IPi();
bool IntervalEnvironmentSupported();
struct IV
{
    I x, y, z;
};
IV operator+(IV, IV);
IV operator-(IV, IV);
IV operator*(IV, I);
IV Exact(V);
IV Normalize(IV);
I Dot(IV, IV);
IV Cross(IV, IV);
struct BPoly
{
    unsigned u = 0, v = 0;
    std::vector<I> c{I(0)};
    BPoly() = default;
    BPoly(unsigned a, unsigned b, I x = I(0)) : u(a), v(b), c((a + 1) * (b + 1), x) {}
    I &at(unsigned a, unsigned b)
    {
        return c[a * (v + 1) + b];
    }
    I at(unsigned a, unsigned b) const
    {
        return c[a * (v + 1) + b];
    }
    I range() const;
};
BPoly operator+(const BPoly &, const BPoly &);
BPoly operator-(const BPoly &, const BPoly &);
BPoly operator*(const BPoly &, const BPoly &);
BPoly operator*(const BPoly &, I);
BPoly Constant(const BPoly &, I);
BPoly Derivative(const BPoly &, bool);
BPoly Evaluate(const BPoly &, double, double);
std::pair<BPoly, BPoly> Split(const BPoly &, bool);
struct BV
{
    BPoly x, y, z;
};
BV operator+(const BV &, const BV &);
BV operator-(const BV &, const BV &);
BV operator*(const BV &, I);
BV operator*(const BV &, const BPoly &);
BV Constant(const BV &, IV);
BPoly Dot(const BV &, const BV &);
BPoly Dot(const BV &, IV);
BV Cross(const BV &, const BV &);
BV Derivative(const BV &, bool);
BV Evaluate(const BV &, double, double);
IV Range(const BV &);
struct BezierCell
{
    BV p;
    double u0 = 0, u1 = 1, v0 = 0, v1 = 1;
    // 有理输入中 p 为齐次分子；多项式输入的权重 w=1。
    BPoly w{0, 0, I(1)};
};
std::vector<BezierCell> ExtractCells(const sggk::BSplineSurface &);
std::pair<BezierCell, BezierCell> Split(const BezierCell &, bool);
BV Boundary(const BV &, BoundaryId);
// 与 Boundary 压成一维的返回值不同，此处保留源 U/V 轴及其参数域。
BezierCell BoundaryCell(const BezierCell &, BoundaryId);
I AngleRange(const BPoly &x, const BPoly &y, double reference);
} // namespace bspline2analytic::detail
