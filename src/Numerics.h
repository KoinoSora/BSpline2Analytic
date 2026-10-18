/**
 * @file Numerics.h
 * @brief 精确识别内部数值类型与辅助函数声明。
 * 提供向量、失败异常及双精度 Bernstein 表示；此表示本身不等同于区间证书。
 */
#pragma once
#include "BSpline2Analytic.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace bspline2analytic::detail
{
constexpr double Pi = 3.1415926535897932384626433832795;
struct V
{
    double x = 0, y = 0, z = 0;
    V operator+(V b) const
    {
        return {x + b.x, y + b.y, z + b.z};
    }
    V operator-(V b) const
    {
        return {x - b.x, y - b.y, z - b.z};
    }
    V operator*(double s) const
    {
        return {x * s, y * s, z * s};
    }
    V operator/(double s) const
    {
        return *this * (1 / s);
    }
};
inline V Get(const sggk::Point3D &p)
{
    return {p.X(), p.Y(), p.Z()};
}
inline V Get(const sggk::Vector3D &p)
{
    return {p.X(), p.Y(), p.Z()};
}
inline sggk::Point3D Point(V p)
{
    return {p.x, p.y, p.z};
}
inline double Dot(V a, V b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline V Cross(V a, V b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double Norm(V a)
{
    return std::sqrt(Dot(a, a));
}
inline V Unit(V a)
{
    double n = Norm(a);
    if (!(n > 0))
        throw std::runtime_error("Zero direction");
    return a / n;
}
inline sggk::Dir3D Direction(V a)
{
    a = Unit(a);
    return {a.x, a.y, a.z};
}
inline bool Finite(V a)
{
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}
inline V Perpendicular(V a)
{
    return Unit(Cross(a, std::abs(a.x) < 0.8 ? V{1, 0, 0} : V{0, 1, 0}));
}
inline bool Parallel(V a, V b, double tol)
{
    return Norm(Cross(a, b)) <= tol;
}
inline double Unwrap(double value, double reference)
{
    return value + 2 * Pi * std::round((reference - value) / (2 * Pi));
}

struct Failure : std::runtime_error
{
    Status status;
    Failure(Status s, const std::string &m) : std::runtime_error(m), status(s) {}
};

// 张量积 Bernstein 多项式；存储中 v 为变化最快的索引。
struct Poly
{
    unsigned u = 0, v = 0;
    std::vector<double> c;
    Poly(unsigned du = 0, unsigned dv = 0, double value = 0)
        : u(du), v(dv), c((du + 1) * (dv + 1), value)
    {
    }
    double &at(unsigned i, unsigned j)
    {
        return c[i * (v + 1) + j];
    }
    double at(unsigned i, unsigned j) const
    {
        return c[i * (v + 1) + j];
    }
    double maxAbs() const;
    double min() const
    {
        return *std::min_element(c.begin(), c.end());
    }
    double max() const
    {
        return *std::max_element(c.begin(), c.end());
    }
};
Poly operator+(const Poly &a, const Poly &b);
Poly operator-(const Poly &a, const Poly &b);
Poly operator*(const Poly &a, const Poly &b);
Poly operator*(const Poly &a, double s);
std::pair<Poly, Poly> Split(const Poly &p, bool alongU);
struct Homogeneous
{
    Poly x, y, z, w;
};
std::pair<Homogeneous, Homogeneous> Split(const Homogeneous &h, bool alongU);
std::vector<Homogeneous> CurvePieces(const sggk::BSplineCurve3D &curve, V origin, double scale);
std::vector<Homogeneous> SurfacePieces(const sggk::BSplineSurface &s, V origin, double scale);
Poly DotNumerator(const Homogeneous &h, V d);
Poly SquaredNorm(const Homogeneous &h);
bool LeastSquares3(std::vector<std::array<double, 3>> rows, std::vector<double> rhs,
                   double rankTolerance, std::array<double, 3> &result);
std::vector<double> SampleParameters(const sggk::BSplineCurve3D &c, unsigned perSpan);
std::vector<double> SampleParameters(const sggk::RealArray &knots, double low, double high,
                                     unsigned perSpan);
void ValidateOptions(const RecognitionOptions &o);
void ValidateSurface(const sggk::BSplineSurface &s);
double Roundoff(double scale) noexcept;
} // namespace bspline2analytic::detail
