#include "AnalyticReconstruction.h"
#include "Numerics.h"
#include <Geometry/3D/Surface/ConeSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <Geometry/3D/Surface/TorusSurface.h>

namespace bspline2analytic::detail
{
static bool Near(V a, V b, double tol)
{
    return Norm(a - b) <= tol;
}
static bool Axial(const CircleData &c, V point, V axis, double tol, double angle)
{
    return Parallel(Get(c.normal), axis, angle) && Norm(Cross(Get(c.center) - point, axis)) <= tol;
}
static sggk::Ucs3D Frame(AnalyticParameters &p, V center, V axis, V x)
{
    axis = Unit(axis);
    x = Unit(x - axis * Dot(x, axis));
    p.origin = Point(center);
    p.axis = Direction(axis);
    p.xDirection = Direction(x);
    return {p.origin, p.xDirection, p.axis};
}
static void Add(std::vector<Candidate> &out, AnalyticParameters p, const sggk::Ucs3D &frame)
{
    sggk::SurfacePtr s;
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        s = std::make_shared<sggk::Plane>(frame);
        break;
    case sggk::SurfaceType::Cylinder:
        s = std::make_shared<sggk::CylinderSurface>(frame, p.radius);
        break;
    case sggk::SurfaceType::Cone:
        s = std::make_shared<sggk::ConeSurface>(
            sggk::ConeSurface::CreateFromHalfAngle(frame, p.halfAngle));
        break;
    case sggk::SurfaceType::Sphere:
        s = std::make_shared<sggk::SphereSurface>(frame, p.radius);
        break;
    case sggk::SurfaceType::Torus:
        s = std::make_shared<sggk::TorusSurface>(frame, p.majorRadius, p.radius);
        break;
    default:
        return;
    }
    out.push_back({p, s});
}
double Distance(const AnalyticParameters &p, const sggk::Point3D &point)
{
    V q = Get(point) - Get(p.origin), a = Get(p.axis);
    double z = Dot(q, a), rho = Norm(q - a * z);
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        return std::abs(z);
    case sggk::SurfaceType::Cylinder:
        return std::abs(rho - p.radius);
    case sggk::SurfaceType::Sphere:
        return std::abs(Norm(q) - p.radius);
    case sggk::SurfaceType::Cone: {
        double t = rho * std::sin(p.halfAngle) + z * std::cos(p.halfAngle);
        return t < 0 ? Norm(q) : std::abs(rho * std::cos(p.halfAngle) - z * std::sin(p.halfAngle));
    }
    case sggk::SurfaceType::Torus:
        return std::abs(std::hypot(rho - p.majorRadius, z) - p.radius);
    default:
        return std::numeric_limits<double>::infinity();
    }
}

std::vector<Candidate> Reconstruct(const std::array<BoundaryInfo, 4> &b,
                                   const RecognitionOptions &o)
{
    std::vector<Candidate> out;
    double tol = o.distanceTolerance, angle = o.angularTolerance;
    bool allLines = true, allCircles = true;
    for (const auto &e : b)
    {
        allLines &= e.kind == BoundaryKind::Line;
        allCircles &= e.kind == BoundaryKind::Circle;
    }
    if (allLines)
    {
        V x = Get(b[2].line.direction), y = Get(b[0].line.direction);
        if (Norm(Cross(x, y)) <= angle)
            return out;
        V n = Unit(Cross(x, y)), origin = Get(b[0].curve->CalcPoint(b[0].curve->MinParam()));
        // SGK plane UV is orthonormal, so a skew parallelogram is not this input class.
        if (std::abs(Dot(x, y)) > angle)
            return out;
        for (unsigned i = 0; i < 4; ++i)
            if (std::abs(Dot(Get(b[i].line.origin) - origin, n)) > tol ||
                !Parallel(Get(b[i].line.direction), i < 2 ? y : x, angle))
                return out;
        AnalyticParameters p;
        p.type = sggk::SurfaceType::Plane;
        auto frame = Frame(p, origin, n, x);
        Add(out, p, frame);
        return out;
    }
    for (unsigned g = 0; g < 2; ++g)
    {
        unsigned c = 2 * g, l = 2 * (1 - g);
        if (b[c].kind != BoundaryKind::Circle || b[c + 1].kind != BoundaryKind::Circle ||
            b[l].kind != BoundaryKind::Line || b[l + 1].kind != BoundaryKind::Line)
            continue;
        const auto &c0 = b[c].circle;
        const auto &c1 = b[c + 1].circle;
        V o0 = Get(c0.center), o1 = Get(c1.center), delta = o1 - o0;
        double h = Norm(delta);
        if (h <= Roundoff(std::max(Norm(o0), Norm(o1))))
            continue;
        V axis = delta / h;
        if (!Axial(c0, o0, axis, tol, angle) || !Axial(c1, o0, axis, tol, angle))
            continue;
        double dr = c1.radius - c0.radius;
        bool parallelLines = Parallel(Get(b[l].line.direction), axis, angle) &&
                             Parallel(Get(b[l + 1].line.direction), axis, angle);
        V first = Get(b[c].curve->CalcPoint(b[c].curve->MinParam())), x = first - o0;
        // A short arc has a less accurately determined radius after knot refinement.
        // Accept radius agreement at geometric tolerance only together with parallel generators.
        const bool parallelGenerators =
            Parallel(Get(b[l].line.direction), Get(b[l + 1].line.direction), angle);
        if (parallelLines && parallelGenerators && std::abs(dr) <= tol)
        {
            bool valid = true;
            for (unsigned i = l; i < l + 2; ++i)
                valid &= std::abs(Norm(Cross(Get(b[i].line.origin) - o0, axis)) - c0.radius) <= tol;
            if (valid)
            {
                AnalyticParameters p;
                p.type = sggk::SurfaceType::Cylinder;
                p.radius = (c0.radius + c1.radius) * 0.5;
                auto f = Frame(p, o0, axis, x);
                Add(out, p, f);
            }
        }
        if (!parallelGenerators && std::abs(dr) > Roundoff(std::max(c0.radius, c1.radius)))
        {
            double slope = dr / h;
            V apex = o0 - axis * (c0.radius / slope), coneAxis = slope > 0 ? axis : axis * (-1);
            bool valid = true;
            for (unsigned i = l; i < l + 2; ++i)
                valid &= Norm(Cross(Get(b[i].line.origin) - apex, Get(b[i].line.direction))) <=
                         tol + Roundoff(Norm(apex));
            if (valid)
            {
                AnalyticParameters p;
                p.type = sggk::SurfaceType::Cone;
                p.halfAngle = std::atan(std::abs(slope));
                auto f = Frame(p, apex, coneAxis, x);
                Add(out, p, f);
            }
        }
    }
    if (!allCircles)
        return out;
    // Try both assignments: meridians vs parallels. Circle-plane normals have arbitrary signs.
    for (unsigned g = 0; g < 2; ++g)
    {
        unsigned m = 2 * g, v = 2 * (1 - g);
        const auto &m0 = b[m].circle;
        const auto &m1 = b[m + 1].circle;
        const auto &v0 = b[v].circle;
        const auto &v1 = b[v + 1].circle;
        V a = Get(v0.normal), p0 = Get(v0.center), M0 = Get(m0.center), M1 = Get(m1.center);
        if (!Axial(v1, p0, a, tol, angle))
            continue;
        if (std::abs(Dot(Get(m0.normal), a)) > angle || std::abs(Dot(Get(m1.normal), a)) > angle)
            continue;
        if (Near(M0, M1, tol) && std::abs(m0.radius - m1.radius) <= tol)
        {
            V center = (M0 + M1) * 0.5;
            double radius = (m0.radius + m1.radius) * 0.5;
            bool valid = true;
            for (const auto &e : b)
            {
                V d = Get(e.circle.center) - center;
                valid &= Norm(Cross(d, Get(e.circle.normal))) <= tol;
                valid &= std::abs(std::hypot(Norm(d), e.circle.radius) - radius) <= tol;
            }
            if (valid)
            {
                V x = Get(b[m].curve->CalcPoint(b[m].curve->MinParam())) - center;
                x = x - a * Dot(x, a);
                if (Norm(x) > tol)
                {
                    AnalyticParameters p;
                    p.type = sggk::SurfaceType::Sphere;
                    p.radius = radius;
                    auto f = Frame(p, center, a, x);
                    Add(out, p, f);
                }
            }
        }
        V center = p0 + a * Dot(M0 - p0, a), other = p0 + a * Dot(M1 - p0, a);
        V radial0 = M0 - center, radial1 = M1 - center;
        if (!Near(center, other, tol) || std::abs(m0.radius - m1.radius) > tol ||
            std::abs(Norm(radial0) - Norm(radial1)) > tol ||
            Norm(radial0) <= Roundoff(Norm(center)))
            continue;
        if (std::abs(Dot(radial0, Get(m0.normal))) > tol ||
            std::abs(Dot(radial1, Get(m1.normal))) > tol)
            continue;
        V point = Get(b[m].curve->CalcPoint(b[m].curve->MinParam())), x = point - center;
        x = x - a * Dot(x, a);
        if (Norm(x) <= tol)
            continue;
        x = Unit(x);
        double major = Dot(radial0, x), tube = (m0.radius + m1.radius) * 0.5;
        bool valid = true;
        for (unsigned i = v; i < v + 2; ++i)
        {
            double z = Dot(Get(b[i].circle.center) - center, a);
            valid &= std::abs(std::hypot(b[i].circle.radius - major, z) - tube) <= tol;
        }
        if (valid)
        {
            AnalyticParameters p;
            p.type = sggk::SurfaceType::Torus;
            p.radius = tube;
            p.majorRadius = major;
            auto f = Frame(p, center, a, x);
            Add(out, p, f);
        }
    }
    return out;
}
} // namespace bspline2analytic::detail
