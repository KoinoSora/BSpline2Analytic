/**
 * @file GeneralGeometry.cpp
 * @brief 一般有限解析面片的几何对象实现。
 * 关联母面参数、规范角点和 SDK PCurve；母面对象自身不含有限裁剪边界。
 */
#include "GeneralInternal.h"
namespace bspline2analytic::detail
{
V ProjectToMother(const AnalyticParameters &a, V world)
{
    auto f = FrameOf(a);
    V p = f.local(world), q;
    double rho = std::hypot(p.x, p.y);
    switch (a.type)
    {
    case sggk::SurfaceType::Plane:
        q = {p.x, p.y, 0}; break;
    case sggk::SurfaceType::Sphere:
        if (!(a.radius > 0) || !(Norm(p) > 0))
            throw Failure(Status::IllConditioned, "Sphere projection at center");
        q = Unit(p) * a.radius; break;
    case sggk::SurfaceType::Cylinder:
        if (!(rho > 0) || !(a.radius > 0))
            throw Failure(Status::IllConditioned, "Cylinder projection on axis");
        q = {p.x * a.radius / rho, p.y * a.radius / rho, p.z}; break;
    case sggk::SurfaceType::Cone: {
        if (!(a.halfAngle > 0 && a.halfAngle < Pi / 2))
            throw Failure(Status::InvalidInput, "Invalid cone angle");
        double c = std::cos(a.halfAngle), s = std::sin(a.halfAngle);
        double t = std::max(0.0, rho * s + p.z * c);
        if (t == 0) q = {};
        else
        {
            if (!(rho > 0))
                throw Failure(Status::IllConditioned, "Cone projection on positive axis");
            q = {p.x * t * s / rho, p.y * t * s / rho, t * c};
        }
        break;
    }
    case sggk::SurfaceType::Torus: {
        if (!(a.majorRadius > a.radius && a.radius > 0) || !(rho > 0))
            throw Failure(Status::IllConditioned, "Nonregular ring-torus projection");
        double dx = rho - a.majorRadius, length = std::hypot(dx, p.z);
        if (!(length > 0))
            throw Failure(Status::IllConditioned, "Torus projection on tube centerline");
        double r = a.majorRadius + a.radius * dx / length;
        q = {p.x * r / rho, p.y * r / rho, a.radius * p.z / length};
        break;
    }
    default: throw Failure(Status::Unsupported, "Unsupported general mother surface");
    }
    V out = f.global(q);
    if (!Finite(out)) throw Failure(Status::InvalidInput, "Nonfinite mother projection");
    return out;
}
}
namespace bspline2analytic
{
sggk::Point2D LinearPCurve::CalcPoint(double t) const
{
    if (!std::isfinite(t) || t < 0 || t > 1)
        throw std::out_of_range("PCurve fraction outside [0,1]");
    return {start.X() + (end.X() - start.X()) * t,
            start.Y() + (end.Y() - start.Y()) * t};
}
GeneralAnalyticPatch::GeneralAnalyticPatch(const AnalyticParameters &a,
    const std::array<double, 4> &sourceDomain, const std::array<sggk::Point2D, 4> &uvCorners,
    const std::array<LinearPCurve, 4> &edges)
    : parameters(a), domain(sourceDomain), corners(uvCorners), pcurves(edges)
{
    AnalyticPatch mother;
    mother.parameters = parameters;
    detail::MakeNative(mother);
    motherSurface = mother.surface;
    parameters = mother.parameters;
}
sggk::Point3D GeneralAnalyticPatch::CalcPoint(double u, double v) const
{
    detail::GeneralTrim trim; trim.parameters = parameters; trim.corners = corners; trim.edges = pcurves;
    if (!std::isfinite(u) || !std::isfinite(v) || !detail::InsideTrim(trim, {u, v}, 1e-12))
        throw std::out_of_range("Canonical UV outside finite analytic trim");
    return detail::Point(detail::FrameOf(parameters).global(detail::CanonicalPoint(parameters, u, v)));
}
sggk::Point3D GeneralAnalyticPatch::CalcBoundaryPoint(BoundaryId edge, double fraction) const
{
    if (!std::isfinite(fraction) || fraction < 0 || fraction > 1)
        throw std::out_of_range("Boundary fraction outside [0,1]");
    unsigned k = static_cast<unsigned>(edge);
    if (k > 3) throw std::out_of_range("Invalid boundary identifier");
    auto uv = pcurves[k].CalcPoint(fraction);
    return detail::Point(detail::FrameOf(parameters).global(
        detail::CanonicalPoint(parameters, uv.X(), uv.Y())));
}
}
