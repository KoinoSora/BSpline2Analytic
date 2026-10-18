/**
 * @file ReferenceCases.cpp
 * @brief 精确识别测试的五类解析参考曲面构造。
 * 只由测试侧持有真值，转换器从实际 NURBS 输入恢复几何。
 */
#include "ReferenceCases.h"
#include <Geometry/3D/Surface/BezierSurface.h>
#include <Geometry/3D/Surface/ConeSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <Geometry/3D/Surface/TorusSurface.h>

namespace bspline2analytic::tests
{
std::vector<ReferenceCase> MakeCases()
{
    std::vector<ReferenceCase> cases;
    sggk::Ucs3D f(sggk::Point3D(13, -21, 8), sggk::Dir3D(1, 0, 1), sggk::Dir3D(-1, 2, 1));
    auto add = [&](std::string n, sggk::SurfacePtr a, double u0, double u1, double v0, double v1) {
        auto s = a->ToBSpline(sggk::UVRange(sggk::Interval(u0, u1), sggk::Interval(v0, v1)));
        cases.push_back({n, a, s});
    };
    add("plane", std::make_shared<sggk::Plane>(f), -3, 4, -2, 5);
    add("cylinder", std::make_shared<sggk::CylinderSurface>(f, 3), 0.2, 1.4, -2, 3);
    add("cone", std::make_shared<sggk::ConeSurface>(sggk::ConeSurface::CreateFromHalfAngle(f, 0.4)),
        0.2, 1.3, 2, 5);
    add("sphere", std::make_shared<sggk::SphereSurface>(f, 4), 0.3, 1.5, 0.5, 2.2);
    add("torus", std::make_shared<sggk::TorusSurface>(f, 5, 1), 0.2, 1.4, 0.2, 2.8);
    add("torus-concentric", std::make_shared<sggk::TorusSurface>(f, 3, 1), 0.2, 1.1, 0,
        3.141592653589793);
    add("sphere-antipodal-meridians", std::make_shared<sggk::SphereSurface>(f, 4), 0.2,
        3.341592653589793, 0.4, 2.3);
    add("shallow-cone",
        std::make_shared<sggk::ConeSurface>(sggk::ConeSurface::CreateFromHalfAngle(f, 0.005)), 0.2,
        1.1, 200, 400);
    add("short-cylinder-arc", std::make_shared<sggk::CylinderSurface>(f, 3), 0.3, 0.32, 1, 5);
    add("torus-apple", std::make_shared<sggk::TorusSurface>(f, 1, 3), 0.2, 1.1, 0.2, 0.9);
    add("torus-lemon", std::make_shared<sggk::TorusSurface>(f, -1, 3), 0.2, 1.1, 0.2, 0.9);
    add("cylinder-full-turn", std::make_shared<sggk::CylinderSurface>(f, 3), 0, 6.283185307179586,
        1, 5);
    add("torus-cross-seams", std::make_shared<sggk::TorusSurface>(f, 5, 1), 5.8, 6.8, 5.5, 6.9);
    add("torus-full-turns", std::make_shared<sggk::TorusSurface>(f, 5, 1), 0, 6.283185307179586, 0,
        6.283185307179586);
    sggk::Ucs3D big(sggk::Point3D(1e6, -2e6, 3e6), sggk::Dir3D(1, 0, 1), sggk::Dir3D(-1, 2, 1));
    add("large-origin-sphere", std::make_shared<sggk::SphereSurface>(big, 4), 0.3, 1.5, 0.5, 2.2);
    return cases;
}

sggk::BSplineSurfacePtr ElevateBezierPatch(const sggk::BSplineSurface &s)
{
    auto b = s.ToBezierSurfaces().front().front();
    const auto &p = b->ControlPoints();
    const auto &w = b->Weights();
    using H = std::array<double, 4>;
    std::vector<std::vector<H>> h(p.size(), std::vector<H>(p.front().size()));
    for (std::size_t i = 0; i < p.size(); ++i)
        for (std::size_t j = 0; j < p[i].size(); ++j)
        {
            double weight = w.empty() ? 1 : w[i][j];
            h[i][j] = {p[i][j].X() * weight, p[i][j].Y() * weight, p[i][j].Z() * weight, weight};
        }
    auto blend = [](const H &a, const H &b, double t) {
        H r;
        for (unsigned k = 0; k < 4; ++k)
            r[k] = t * a[k] + (1 - t) * b[k];
        return r;
    };
    unsigned nu = static_cast<unsigned>(h.size()), nv = static_cast<unsigned>(h[0].size());
    std::vector<std::vector<H>> hu(nu + 1, std::vector<H>(nv));
    hu.front() = h.front();
    hu.back() = h.back();
    for (unsigned i = 1; i < nu; ++i)
        for (unsigned j = 0; j < nv; ++j)
            hu[i][j] = blend(h[i - 1][j], h[i][j], double(i) / nu);
    std::vector<std::vector<H>> hv(nu + 1, std::vector<H>(nv + 1));
    for (unsigned i = 0; i <= nu; ++i)
    {
        hv[i].front() = hu[i].front();
        hv[i].back() = hu[i].back();
        for (unsigned j = 1; j < nv; ++j)
            hv[i][j] = blend(hu[i][j - 1], hu[i][j], double(j) / nv);
    }
    sggk::Point3DMatrix points(nu + 1, std::vector<sggk::Point3D>(nv + 1));
    sggk::RealMatrix weights(nu + 1, std::vector<double>(nv + 1));
    for (unsigned i = 0; i <= nu; ++i)
        for (unsigned j = 0; j <= nv; ++j)
        {
            const auto &q = hv[i][j];
            points[i][j] = sggk::Point3D(q[0] / q[3], q[1] / q[3], q[2] / q[3]);
            weights[i][j] = q[3];
        }
    return std::make_shared<sggk::BSplineSurface>(
        nu, nv, points, weights, sggk::RealArray{0, 1}, sggk::RealArray{0, 1},
        sggk::UIntArray{nu + 1, nu + 1}, sggk::UIntArray{nv + 1, nv + 1});
}
} // namespace bspline2analytic::tests
