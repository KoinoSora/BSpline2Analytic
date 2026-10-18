/**
 * @file ConversionValidation.cpp
 * @brief 已废弃的等参投影验证器实现。
 * 仅保留历史算法与接口，不运行其测试、预热或性能对照；当前配对及母面验证器分别独立实现。
 */
#include "ApproximateInternal.h"
#include "ValidationTiming.h"
#include <functional>
namespace bspline2analytic::detail
{
static I Abs(I a)
{
    return {a.lo <= 0 && a.hi >= 0 ? 0 : std::min(std::abs(a.lo), std::abs(a.hi)), a.absMax()};
}
static BV Local(const BV &p, const AnalyticParameters &a)
{
    IV z = Normalize(Exact(Get(a.axis))), x = Exact(Get(a.xDirection));
    x = Normalize(x - z * Dot(x, z));
    IV y = Cross(z, x);
    BV d = p - Constant(p, Exact(Get(a.origin)));
    return {Dot(d, x), Dot(d, y), Dot(d, z)};
}
struct CertifiedEdge
{
    EdgeGeometry nominal;
    IV o, d, ex, ey, n0, nx, ny;
    I radius, g0, g1;
};
static CertifiedEdge CertifyEdge(const AnalyticPatch &patch, const EdgeGeometry &e)
{
    CertifiedEdge r;
    r.nominal = e;
    const auto &a = patch.parameters;
    I c = Cos(I(e.value)), s = Sin(I(e.value)), R(a.radius);
    IV x{I(1), I(0), I(0)}, y{I(0), I(1), I(0)}, z{I(0), I(0), I(1)};
    if (a.type == sggk::SurfaceType::Plane)
    {
        r.o = e.fixed == 0 ? x * I(e.value) : y * I(e.value);
        r.d = e.fixed == 0 ? y : x;
        r.n0 = z;
    }
    else if (!e.circle)
    {
        IV radial = x * c + y * s;
        if (a.type == sggk::SurfaceType::Cylinder)
        {
            r.o = radial * R;
            r.d = z;
            r.n0 = radial;
        }
        else
        {
            I ca = Cos(I(a.halfAngle)), sa = Sin(I(a.halfAngle));
            r.d = radial * (sa / ca) + z;
            r.n0 = radial * ca - z * sa;
        }
    }
    else if (e.fixed == 1)
    {
        r.ex = x;
        r.ey = y;
        switch (a.type)
        {
        case sggk::SurfaceType::Cylinder:
            r.radius = R;
            r.o = z * I(e.value);
            r.nx = x;
            r.ny = y;
            break;
        case sggk::SurfaceType::Cone: {
            I ca = Cos(I(a.halfAngle)), sa = Sin(I(a.halfAngle));
            r.radius = I(e.value) * sa / ca;
            r.o = z * I(e.value);
            r.n0 = z * (-sa);
            r.nx = x * ca;
            r.ny = y * ca;
            break;
        }
        case sggk::SurfaceType::Sphere:
            r.radius = R * s;
            r.o = z * (-R * c);
            r.n0 = z * (-c);
            r.nx = x * s;
            r.ny = y * s;
            break;
        case sggk::SurfaceType::Torus:
            r.radius = I(a.majorRadius) + R * c;
            r.o = z * (R * s);
            r.n0 = z * s;
            r.nx = x * c;
            r.ny = y * c;
            break;
        default:
            throw Failure(Status::Unsupported, "Unsupported boundary chart");
        }
    }
    else if (a.type == sggk::SurfaceType::Sphere)
    {
        r.radius = R;
        r.ex = z * I(-1);
        r.ey = x * c + y * s;
        r.nx = r.ex;
        r.ny = r.ey;
    }
    else
    {
        r.radius = R;
        r.ex = x * c + y * s;
        r.ey = z;
        r.o = r.ex * I(a.majorRadius);
        r.nx = r.ex;
        r.ny = r.ey;
    }
    return r;
}
static I Coordinate(const BV &p, const CertifiedEdge &e, double reference)
{
    BV d = p - Constant(p, e.o);
    if (e.nominal.circle)
        return AngleRange(Dot(d, e.ex), Dot(d, e.ey), reference);
    return Dot(d, e.d).range() / Dot(e.d, e.d);
}
static double Reference(double mid, double lo, double hi, bool periodic)
{
    if (periodic)
        mid = Unwrap(mid, (lo + hi) * 0.5);
    // 两个端点都接近时选更近者，使稳定的横向坐标公式在其零点附近更紧。
    // 固定选择下端点会使狭窄管截面角区间的上边界过度细分。
    if (std::abs(mid - lo) < Pi / 4 && std::abs(mid - lo) <= std::abs(mid - hi))
        return lo;
    if (std::abs(mid - hi) < Pi / 4)
        return hi;
    return mid;
}
static double Excess(I a, double lo, double hi)
{
    if (!a.finite())
        return INFINITY;
    return std::max({0.0, (I(lo) - a).hi, (a - I(hi)).hi});
}
static I TubeAngle(const BPoly &radial, const BPoly &z, I rho, I major, double reference)
{
    I c = Cos(I(reference)), s = Sin(I(reference));
    I forward = (rho - major) * c + z.range() * s;
    if (forward.lo <= 0)
        return I::Whole();
    // 旋转后的横向坐标为 b=z*cos(ref)+R*sin(ref)-rho*sin(ref)。
    // 仅对分子平方可保留管截面角边界附近 z 与 rho² 的相关性；直接旋转独立的 (rho,z) 区间盒会丢失相关性，导致过度细分。
    BPoly line = z * c + Constant(z, major * s);
    I denominator = line.range() + rho * s;
    I cross;
    if (denominator.lo > 0 || denominator.hi < 0)
        cross = (line * line - radial * (s * s)).range() / denominator;
    else
        cross = line.range() - rho * s;
    return I(reference) + Atan(cross / forward);
}
static I SurfaceDistance(const BV &p, const AnalyticPatch &patch)
{
    const auto &a = patch.parameters;
    BPoly radial = p.x * p.x + p.y * p.y, norm = radial + p.z * p.z;
    I rho = Sqrt(radial.range()), n = Sqrt(norm.range()), z = p.z.range(), r(a.radius);
    I distance;
    switch (a.type)
    {
    case sggk::SurfaceType::Plane:
        distance = Abs(z);
        break;
    case sggk::SurfaceType::Cylinder:
        distance = Abs((radial - Constant(radial, r * r)).range()) / (rho + r);
        break;
    case sggk::SurfaceType::Sphere:
        distance = Abs((norm - Constant(norm, r * r)).range()) / (n + r);
        break;
    case sggk::SurfaceType::Cone: {
        I ca = Cos(I(a.halfAngle)), sa = Sin(I(a.halfAngle));
        I divisor = rho * ca + z * sa;
        if (divisor.lo <= 0 || (rho * sa + z * ca).lo <= 0)
            return {0, INFINITY};
        distance = Abs((radial * (ca * ca) - (p.z * p.z) * (sa * sa)).range()) / divisor;
        break;
    }
    case sggk::SurfaceType::Torus: {
        I R(a.majorRadius);
        BPoly sum = norm + Constant(norm, R * R - r * r);
        BPoly f = sum * sum - radial * (I(4) * R * R);
        I other = Abs(Sqr(rho + R) + Sqr(z) - r * r);
        if (other.lo <= 0)
            return {0, INFINITY};
        I q = Sqrt(Sqr(rho - R) + Sqr(z));
        distance = Abs(f.range()) / (other * (q + r));
        break;
    }
    default:
        return {0, INFINITY};
    }
    double du = 0, dv = 0, speedU = 1, speedV = 1;
    auto point = Range(Evaluate(p, 0.5, 0.5));
    V mid{(point.x.lo + point.x.hi) * 0.5, (point.y.lo + point.y.hi) * 0.5,
          (point.z.lo + point.z.hi) * 0.5};
    auto uv = MotherParameters(a, mid);
    if (a.type == sggk::SurfaceType::Plane)
    {
        du = Excess(p.x.range(), patch.domain[0], patch.domain[1]);
        dv = Excess(p.y.range(), patch.domain[2], patch.domain[3]);
    }
    else
    {
        if ((I(patch.domain[1]) - I(patch.domain[0])).lo < (IPi() * I(2)).hi)
        {
            I u = AngleRange(p.x, p.y, Reference(uv[0], patch.domain[0], patch.domain[1], true));
            du = Excess(u, patch.domain[0], patch.domain[1]);
        }
        I v;
        switch (a.type)
        {
        case sggk::SurfaceType::Cylinder:
            v = z;
            speedU = a.radius;
            break;
        case sggk::SurfaceType::Cone: {
            I t = Sin(I(a.halfAngle)) / Cos(I(a.halfAngle));
            v = (rho * t + z) / (I(1) + t * t);
            speedU = ((I(std::max(std::abs(v.lo), std::abs(v.hi))) +
                       I(std::max(std::abs(patch.domain[2]), std::abs(patch.domain[3])))) *
                      t)
                         .hi;
            speedV = Sqrt(I(1) + t * t).hi;
            break;
        }
        case sggk::SurfaceType::Sphere:
            v = Acos(-z / n);
            speedU = speedV = a.radius;
            break;
        case sggk::SurfaceType::Torus: {
            double ref = Reference(uv[1], patch.domain[2], patch.domain[3], true);
            v = TubeAngle(radial, p.z, rho, I(a.majorRadius), ref);
            speedU = (I(std::abs(a.majorRadius)) + r).hi;
            speedV = a.radius;
            break;
        }
        default:
            break;
        }
        if (!(a.type == sggk::SurfaceType::Torus &&
              (I(patch.domain[3]) - I(patch.domain[2])).lo >= (IPi() * I(2)).hi))
            dv = Excess(v, patch.domain[2], patch.domain[3]);
    }
    if (!std::isfinite(du) || !std::isfinite(dv))
        return {std::max(0.0, distance.lo), INFINITY};
    I upper = I(distance.hi) + I(du) * I(speedU) + I(dv) * I(speedV);
    return {std::max(0.0, distance.lo), upper.hi};
}
static I DirectionAngle(const BV &a, const BV &b)
{
    I d = Dot(a, b).range();
    auto c = Cross(a, b);
    I cross = Sqrt(Dot(c, c).range());
    if (d.lo <= 0)
        return {d.hi <= 0 ? IPi().lo / 2 : 0, IPi().hi};
    return Atan(cross / d);
}
struct EdgeBounds
{
    I distance, tangent, normal;
    bool monotone = false;
};
static EdgeBounds BoundEdge(const BV &p, const BV &t, const BV &normal, const CertifiedEdge &e,
                            double fraction, int sign, double halfFraction = 0)
{
    EdgeBounds result{{0, INFINITY}, {0, IPi().hi}, {0, IPi().hi}, false};
    auto d = p - Constant(p, e.o);
    I g = Coordinate(p, e, e.nominal.g0 + (e.nominal.g1 - e.nominal.g0) * fraction);
    I extent = e.g1 - e.g0;
    if (!g.finite() || extent.lo <= 0 && extent.hi >= 0)
        return result;
    if (e.nominal.circle)
    {
        // 相对连续参考角选择唯一的周期提升，避免相邻已认证边界区间之间隐藏 2*pi 跳变。
        I reference = I(e.nominal.g0) + I(e.nominal.g1 - e.nominal.g0) *
                                            I(fraction - halfFraction, fraction + halfFraction);
        if ((g - reference).absMax() >= IPi().lo)
            return result;
    }
    I ratio = (I(e.nominal.end) - I(e.nominal.start)) / extent;
    // 与 phi(g)-g 代数等价，但保留小的端点修正量，避免两个大相位区间相减。
    I delta = (I(e.nominal.start) - e.g0) + (ratio - I(1)) * (g - e.g0);
    BPoly derivative;
    BV targetT, targetN;
    I distance, correction(0);
    double tangentCorrection = 0;
    if (e.nominal.circle)
    {
        BPoly x = Dot(d, e.ex), y = Dot(d, e.ey);
        BPoly rho2 = x * x + y * y;
        I rho = Sqrt(rho2.range());
        if (rho.lo <= 0)
            return result;
        BPoly h = Dot(d, Cross(e.ex, e.ey));
        I dr = Abs((rho2 - Constant(rho2, e.radius * e.radius)).range()) / (rho + e.radius);
        I support2 = Sqr(h.range()) + Sqr(dr);
        I upper2 = I(support2.hi) + rho * e.radius * Sqr(delta);
        distance = {Sqrt(support2).lo, Sqrt(upper2).hi};
        targetT = Constant(p, e.ey) * x - Constant(p, e.ex) * y;
        targetN = Constant(p, e.n0 * e.radius) + Constant(p, e.nx) * x + Constant(p, e.ny) * y;
        derivative = x * Dot(t, e.ey) - y * Dot(t, e.ex);
        I targetLength = Sqrt(Dot(targetN, targetN).range());
        I error = Sqrt(Dot(e.n0, e.n0)) * I(dr.hi);
        if (targetLength.lo <= error.hi)
            return result;
        correction = Atan(I(error.hi) / Sqrt(I(targetLength.lo) * I(targetLength.lo) -
                                             I(error.hi) * I(error.hi)));
        tangentCorrection = delta.absMax();
    }
    else
    {
        I speed = Sqrt(Dot(e.d, e.d));
        auto cross = Cross(d, Constant(d, e.d));
        I support2 = Dot(cross, cross).range() / Dot(e.d, e.d);
        distance = {Sqrt(support2).lo, Sqrt(I(std::max(0.0, support2.hi)) + Sqr(delta * speed)).hi};
        targetT = Constant(p, e.d);
        targetN = Constant(p, e.n0);
        derivative = Dot(t, e.d);
    }
    I monotonic = (derivative * I(extent.lo > 0 ? 1 : -1)).range();
    result.monotone = monotonic.lo > 0;
    targetT = targetT * I(e.nominal.end > e.nominal.start ? 1 : -1);
    targetN = targetN * I(sign);
    I ta = DirectionAngle(t, targetT), na = DirectionAngle(normal, targetN);
    I tc(tangentCorrection), nc = correction + tc;
    result.distance = distance;
    result.tangent = {std::max(0.0, (ta - tc).lo), std::min(IPi().hi, (ta + tc).hi)};
    result.normal = {std::max(0.0, (na - nc).lo), std::min(IPi().hi, (na + nc).hi)};
    return result;
}
static void Record(ConstraintReport &r, I value, double tolerance,
                   const std::array<double, 4> &region, BoundaryId edge, bool witness)
{
    r.upperBound = std::max(r.upperBound, value.hi);
    if (witness && value.lo >= tolerance)
    {
        if (r.state != CheckState::Exceeded)
        {
            r.parameterRegion = region;
            r.boundary = edge;
            r.hasWitness = true;
        }
        r.state = CheckState::Exceeded;
    }
    else if (!(value.hi < tolerance) && r.state != CheckState::Exceeded)
    {
        r.state = CheckState::Unresolved;
        r.parameterRegion = region;
        r.boundary = edge;
    }
}
static void CheckPatch(const AnalyticPatch &p)
{
    for (double v : p.domain)
        if (!std::isfinite(v))
            throw Failure(Status::InvalidInput, "Nonfinite analytic domain");
    if (!(p.domain[0] < p.domain[1] && p.domain[2] < p.domain[3]) ||
        !(p.normalSign == 1 || p.normalSign == -1))
        throw Failure(Status::InvalidInput, "Invalid analytic domain or orientation");
    auto f = FrameOf(p.parameters);
    if (!Finite(f.o) || !Finite(f.x) || !Finite(f.z))
        throw Failure(Status::InvalidInput, "Nonfinite analytic frame");
    auto &a = p.parameters;
    if (a.type != sggk::SurfaceType::Plane && a.type != sggk::SurfaceType::Cone &&
        !(a.radius > 0 && std::isfinite(a.radius)))
        throw Failure(Status::InvalidInput, "Invalid analytic radius");
    if (a.type == sggk::SurfaceType::Cone &&
        !(a.halfAngle > 0 && a.halfAngle < Pi / 2 && p.domain[2] > 0))
        throw Failure(Status::Unsupported, "Cone apex or invalid cone chart");
    if (a.type == sggk::SurfaceType::Sphere && !(p.domain[2] > 0 && p.domain[3] < Pi))
        throw Failure(Status::Unsupported, "Sphere pole chart");
    if (a.type != sggk::SurfaceType::Plane && p.domain[1] - p.domain[0] > 2 * Pi + 1e-10)
        throw Failure(Status::Unsupported, "Multiple-turn U chart");
    if (a.type == sggk::SurfaceType::Torus)
    {
        if (!std::isfinite(a.majorRadius) || p.domain[3] - p.domain[2] > 2 * Pi + 1e-10)
            throw Failure(Status::Unsupported, "Invalid or multiple-turn torus chart");
        if ((I(a.majorRadius) + I(a.radius) * Cos(I(p.domain[2], p.domain[3]))).lo <= 0)
        {
            // 长圆弧上的三角函数区间可能过宽；对常规环形圆环面，R-r>0 已足以证明相应正则性。
            if (!(a.majorRadius > a.radius))
                throw Failure(Status::Unsupported, "Torus chart not certified regular");
        }
    }
}
} // namespace bspline2analytic::detail
namespace bspline2analytic
{
ConversionValidationReport ValidateConversion(const sggk::BSplineSurface &s,
                                              const AnalyticPatch &patch,
                                              const ValidationOptions &o)
{
    using namespace detail;
    ConversionValidationReport report;
    ValidationTimer timer(report.timing);
    auto record = [&](ConstraintReport &r, I value, double tolerance,
                      const std::array<double, 4> &region, BoundaryId edge, bool witness) {
        ValidationTimer::Scope stage(timer, report.timing.reporting);
        Record(r, value, tolerance, region, edge, witness);
    };
    try
    {
        if (!IntervalEnvironmentSupported())
            throw Failure(Status::VerificationFailed,
                          "Interval validation requires round-to-nearest and gradual underflow");
        CheckApproximateInput(s, o);
        CheckPatch(patch);
        const size_t nu = s.KnotsU().size() - 1, nv = s.KnotsV().size() - 1;
        if (nu > o.maxValidationPieces || nv > o.maxValidationPieces / nu)
            throw Failure(Status::VerificationFailed,
                          "Initial knot cells exceed the validation budget");
        auto cells = ExtractCells(s);
        if (cells.empty())
            throw Failure(Status::InvalidInput, "No nonzero knot rectangles");
        auto edges = BoundaryMaps(s, patch);
        std::array<CertifiedEdge, 4> certified;
        std::array<BV, 4> corners;
        std::array<bool, 4> foundCorner{};
        for (const auto &cell : cells)
        {
            auto local = Local(cell.p, patch.parameters);
            if (cell.u0 == s.MinParamU() && cell.v0 == s.MinParamV())
            {
                corners[0] = Evaluate(local, 0, 0);
                foundCorner[0] = true;
            }
            if (cell.u1 == s.MaxParamU() && cell.v0 == s.MinParamV())
            {
                corners[1] = Evaluate(local, 1, 0);
                foundCorner[1] = true;
            }
            if (cell.u0 == s.MinParamU() && cell.v1 == s.MaxParamV())
            {
                corners[2] = Evaluate(local, 0, 1);
                foundCorner[2] = true;
            }
            if (cell.u1 == s.MaxParamU() && cell.v1 == s.MaxParamV())
            {
                corners[3] = Evaluate(local, 1, 1);
                foundCorner[3] = true;
            }
        }
        if (!std::all_of(foundCorner.begin(), foundCorner.end(), [](bool b) { return b; }))
            throw Failure(Status::VerificationFailed,
                          "SDK domain extends outside the polynomial knot domain");
        const unsigned ends[4][2] = {{0, 2}, {1, 3}, {0, 1}, {2, 3}};
        for (unsigned k = 0; k < 4; ++k)
        {
            certified[k] = CertifyEdge(patch, edges[k]);
            certified[k].g0 = Coordinate(corners[ends[k][0]], certified[k], edges[k].g0);
            certified[k].g1 = Coordinate(corners[ends[k][1]], certified[k], edges[k].g1);
            if (!certified[k].g0.finite() || !certified[k].g1.finite())
                throw Failure(Status::VerificationFailed,
                              "Unresolved boundary endpoint phase at edge " + std::to_string(k) +
                                  " (references " + std::to_string(edges[k].g0) + ", " +
                                  std::to_string(edges[k].g1) + ")");
        }
        {
            ValidationTimer::Scope stage(timer, report.timing.reporting);
            for (auto *r : {&report.geometry, &report.boundaryTangent, &report.boundaryNormal})
            {
                r->state = CheckState::Passed;
                r->upperBound = 0;
            }
            report.correspondence = CheckState::Passed;
        }
        const I degree = I(180) / IPi();
        timer.Switch(report.timing.subdivision);
        auto region = [](const BezierCell &c) {
            return std::array<double, 4>{c.u0, c.u1, c.v0, c.v1};
        };
        std::function<void(const BezierCell &, unsigned)> surface;
        surface = [&](const BezierCell &c, unsigned depth) {
            if (report.piecesVisited >= o.maxValidationPieces)
            {
                record(report.geometry, I(0, INFINITY), o.distanceTolerance, region(c),
                       BoundaryId::Umin, false);
                return;
            }
            ++report.piecesVisited;
            I distance, witness;
            {
                ValidationTimer::Scope stage(timer, report.timing.surfaceBounds);
                BV p = Local(c.p, patch.parameters);
                distance = SurfaceDistance(p, patch);
                witness = SurfaceDistance(Evaluate(p, 0.5, 0.5), patch);
            }
            {
                ValidationTimer::Scope stage(timer, report.timing.reporting);
                report.geometry.maxSampledUpperBound =
                    std::max(report.geometry.maxSampledUpperBound, witness.hi);
            }
            if (distance.hi < o.distanceTolerance)
                record(report.geometry, distance, o.distanceTolerance, region(c), BoundaryId::Umin,
                       false);
            else if (witness.lo >= o.distanceTolerance)
                record(report.geometry, I(witness.lo, distance.hi), o.distanceTolerance,
                       {(c.u0 + c.u1) / 2, (c.u0 + c.u1) / 2, (c.v0 + c.v1) / 2, (c.v0 + c.v1) / 2},
                       BoundaryId::Umin, true);
            else if (depth >= o.maxSubdivisionDepth)
                record(report.geometry, distance, o.distanceTolerance, region(c), BoundaryId::Umin,
                       false);
            else if (report.geometry.state == CheckState::Exceeded)
                record(report.geometry, distance, o.distanceTolerance, region(c), BoundaryId::Umin,
                       false);
            else
            {
                auto split = Split(c, depth % 2 == 0);
                surface(split.first, depth + 1);
                surface(split.second, depth + 1);
            }
        };
        for (const auto &cell : cells)
            surface(cell, 0);
        timer.Switch(report.timing.boundaryChecks);
        std::function<void(const BezierCell &, unsigned, unsigned)> boundary;
        boundary = [&](const BezierCell &c, unsigned edge, unsigned depth) {
            auto area = region(c);
            if (edge < 2)
                area[0] = area[1] = edge == 0 ? s.MinParamU() : s.MaxParamU();
            else
                area[2] = area[3] = edge == 2 ? s.MinParamV() : s.MaxParamV();
            if (report.piecesVisited >= o.maxValidationPieces)
            {
                for (auto *r : {&report.geometry, &report.boundaryTangent, &report.boundaryNormal})
                    record(*r, I(0, INFINITY), r == &report.geometry ? o.distanceTolerance : 0,
                           area, static_cast<BoundaryId>(edge), false);
                report.correspondence = CheckState::Unresolved;
                return;
            }
            ++report.piecesVisited;
            BV p = Local(c.p, patch.parameters), du = Derivative(p, true),
               dv = Derivative(p, false);
            // 先求导再限制到边界，保留横跨边界方向的偏导信息。
            BV curve = Boundary(p, static_cast<BoundaryId>(edge));
            BV tangent = Boundary(edge < 2 ? dv : du, static_cast<BoundaryId>(edge));
            BV normal = Cross(Boundary(du, static_cast<BoundaryId>(edge)),
                              Boundary(dv, static_cast<BoundaryId>(edge)));
            double t = edge < 2 ? (c.v0 + c.v1) * 0.5 : (c.u0 + c.u1) * 0.5;
            double fraction = edge < 2 ? (t - s.MinParamV()) / (s.MaxParamV() - s.MinParamV())
                                       : (t - s.MinParamU()) / (s.MaxParamU() - s.MinParamU());
            double halfFraction = edge < 2 ? (c.v1 - c.v0) / (2 * (s.MaxParamV() - s.MinParamV()))
                                           : (c.u1 - c.u0) / (2 * (s.MaxParamU() - s.MinParamU()));
            auto value = BoundEdge(curve, tangent, normal, certified[edge], fraction,
                                   patch.normalSign, halfFraction);
            auto point =
                BoundEdge(Evaluate(curve, 0.5, 0), Evaluate(tangent, 0.5, 0),
                          Evaluate(normal, 0.5, 0), certified[edge], fraction, patch.normalSign);
            value.tangent = value.tangent * degree;
            value.normal = value.normal * degree;
            point.tangent = point.tangent * degree;
            point.normal = point.normal * degree;
            {
                ValidationTimer::Scope stage(timer, report.timing.reporting);
                report.geometry.maxSampledUpperBound =
                    std::max(report.geometry.maxSampledUpperBound, point.distance.hi);
                report.boundaryTangent.maxSampledUpperBound =
                    std::max(report.boundaryTangent.maxSampledUpperBound, point.tangent.hi);
                report.boundaryNormal.maxSampledUpperBound =
                    std::max(report.boundaryNormal.maxSampledUpperBound, point.normal.hi);
            }
            bool pass = (report.geometry.state == CheckState::Exceeded ||
                         value.distance.hi < o.distanceTolerance) &&
                        (report.boundaryTangent.state == CheckState::Exceeded ||
                         value.tangent.hi < o.boundaryTangentToleranceDeg) &&
                        (report.boundaryNormal.state == CheckState::Exceeded ||
                         value.normal.hi < o.boundaryNormalToleranceDeg) &&
                        value.monotone;
            bool fail = point.distance.lo >= o.distanceTolerance ||
                        point.tangent.lo >= o.boundaryTangentToleranceDeg ||
                        point.normal.lo >= o.boundaryNormalToleranceDeg;
            if (!pass && !fail && depth < o.maxSubdivisionDepth)
            {
                auto split = Split(c, edge >= 2);
                boundary(split.first, edge, depth + 1);
                boundary(split.second, edge, depth + 1);
                return;
            }
            if (!value.monotone)
                report.correspondence = CheckState::Unresolved;
            auto witness = area;
            if (edge < 2)
                witness[2] = witness[3] = t;
            else
                witness[0] = witness[1] = t;
            auto save = [&](ConstraintReport &r, I bound, I pt, double tol) {
                if (pt.lo >= tol)
                    record(r, I(pt.lo, bound.hi), tol, witness, static_cast<BoundaryId>(edge),
                           true);
                else
                    record(r, bound, tol, area, static_cast<BoundaryId>(edge), false);
            };
            save(report.geometry, value.distance, point.distance, o.distanceTolerance);
            save(report.boundaryTangent, value.tangent, point.tangent,
                 o.boundaryTangentToleranceDeg);
            save(report.boundaryNormal, value.normal, point.normal, o.boundaryNormalToleranceDeg);
        };
        for (const auto &c : cells)
        {
            if (c.u0 == s.MinParamU())
                boundary(c, 0, 0);
            if (c.u1 == s.MaxParamU())
                boundary(c, 1, 0);
            if (c.v0 == s.MinParamV())
                boundary(c, 2, 0);
            if (c.v1 == s.MaxParamV())
                boundary(c, 3, 0);
        }
        timer.Switch(report.timing.reporting);
        report.message = report.Passed() ? "Outward interval/Bernstein bounds passed for "
                                           "finite-patch distance and both boundary angles"
                                         : "A constraint exceeded its limit or could not be "
                                           "bounded within the validation budget";
    }
    catch (const std::exception &e)
    {
        timer.Switch(report.timing.reporting);
        for (auto *r : {&report.geometry, &report.boundaryTangent, &report.boundaryNormal})
            if (r->state != CheckState::Exceeded)
            {
                r->state = CheckState::Unresolved;
                r->upperBound = INFINITY;
            }
        report.correspondence = CheckState::Unresolved;
        report.message = e.what();
    }
    timer.Finish();
    return report;
}
} // namespace bspline2analytic
