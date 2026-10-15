#include "ApproximateInternal.h"
#include <Geometry/3D/Surface/ConeSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <Geometry/3D/Surface/TorusSurface.h>
namespace bspline2analytic::detail
{
Frame3 FrameOf(const AnalyticParameters &p)
{
    V z = Unit(Get(p.axis)), x = Get(p.xDirection);
    x = Unit(x - z * Dot(x, z));
    return {Get(p.origin), x, Cross(z, x), z};
}
void MakeNative(AnalyticPatch &p)
{
    auto f = FrameOf(p.parameters);
    auto &a = p.parameters;
    a.axis = Direction(f.z);
    a.xDirection = Direction(f.x);
    sggk::Ucs3D frame(a.origin, a.xDirection, a.axis);
    switch (a.type)
    {
    case sggk::SurfaceType::Plane:
        p.surface = std::make_shared<sggk::Plane>(frame);
        break;
    case sggk::SurfaceType::Cylinder:
        p.surface = std::make_shared<sggk::CylinderSurface>(frame, a.radius);
        break;
    case sggk::SurfaceType::Cone:
        p.surface = std::make_shared<sggk::ConeSurface>(
            sggk::ConeSurface::CreateFromHalfAngle(frame, a.halfAngle));
        break;
    case sggk::SurfaceType::Sphere:
        p.surface = std::make_shared<sggk::SphereSurface>(frame, a.radius);
        break;
    case sggk::SurfaceType::Torus:
        p.surface = std::make_shared<sggk::TorusSurface>(frame, a.majorRadius, a.radius);
        break;
    default:
        throw Failure(Status::Unsupported, "Unsupported analytic type");
    }
    p.sdkDomain = sggk::UVRange(sggk::Interval(p.domain[0], p.domain[1]),
                                sggk::Interval(p.domain[2], p.domain[3]));
}
V CanonicalPoint(const AnalyticParameters &p, double u, double v)
{
    double c = std::cos(u), s = std::sin(u);
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        return {u, v, 0};
    case sggk::SurfaceType::Cylinder:
        return {p.radius * c, p.radius * s, v};
    case sggk::SurfaceType::Cone:
        return {v * std::tan(p.halfAngle) * c, v * std::tan(p.halfAngle) * s, v};
    case sggk::SurfaceType::Sphere:
        return {p.radius * std::sin(v) * c, p.radius * std::sin(v) * s, -p.radius * std::cos(v)};
    case sggk::SurfaceType::Torus:
        return {(p.majorRadius + p.radius * std::cos(v)) * c,
                (p.majorRadius + p.radius * std::cos(v)) * s, p.radius * std::sin(v)};
    default:
        throw Failure(Status::Unsupported, "Unknown analytic type");
    }
}
V CanonicalNormal(const AnalyticParameters &p, double u, double v)
{
    double c = std::cos(u), s = std::sin(u);
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        return {0, 0, 1};
    case sggk::SurfaceType::Cylinder:
        return {c, s, 0};
    case sggk::SurfaceType::Cone:
        return {std::cos(p.halfAngle) * c, std::cos(p.halfAngle) * s, -std::sin(p.halfAngle)};
    case sggk::SurfaceType::Sphere:
        return {std::sin(v) * c, std::sin(v) * s, -std::cos(v)};
    case sggk::SurfaceType::Torus:
        return {std::cos(v) * c, std::cos(v) * s, std::sin(v)};
    default:
        throw Failure(Status::Unsupported, "Unknown analytic normal");
    }
}
std::array<double, 2> MotherParameters(const AnalyticParameters &p, V a)
{
    double u = std::atan2(a.y, a.x), rho = std::hypot(a.x, a.y);
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        return {a.x, a.y};
    case sggk::SurfaceType::Cylinder:
        return {u, a.z};
    case sggk::SurfaceType::Cone: {
        double t = std::tan(p.halfAngle);
        return {u, std::max(0.0, (rho * t + a.z) / (1 + t * t))};
    }
    case sggk::SurfaceType::Sphere:
        return {u, std::atan2(rho, -a.z)};
    case sggk::SurfaceType::Torus:
        return {u, std::atan2(a.z, rho - p.majorRadius)};
    default:
        throw Failure(Status::Unsupported, "Unknown analytic projection");
    }
}
EdgeGeometry MakeEdge(const AnalyticPatch &p, unsigned edge)
{
    EdgeGeometry e;
    unsigned dimensionU = p.mapping.swappedUV ? 1 : 0;
    e.fixed = edge < 2 ? dimensionU : 1 - dimensionU;
    bool sourceHigh = (edge % 2) == 1;
    bool reverseFixed = edge < 2 ? p.mapping.reversedU : p.mapping.reversedV;
    e.value = p.domain[2 * e.fixed + (sourceHigh != reverseFixed ? 1 : 0)];
    unsigned varying = 1 - e.fixed;
    bool reverseVarying = edge < 2 ? p.mapping.reversedV : p.mapping.reversedU;
    e.start = p.domain[2 * varying + (reverseVarying ? 1 : 0)];
    e.end = p.domain[2 * varying + (reverseVarying ? 0 : 1)];
    const auto &a = p.parameters;
    double c = std::cos(e.value), s = std::sin(e.value);
    if (a.type == sggk::SurfaceType::Plane)
    {
        e.o = e.fixed == 0 ? V{e.value, 0, 0} : V{0, e.value, 0};
        e.d = e.fixed == 0 ? V{0, 1, 0} : V{1, 0, 0};
        e.n0 = {0, 0, 1};
    }
    else if (e.fixed == 0 &&
             (a.type == sggk::SurfaceType::Cylinder || a.type == sggk::SurfaceType::Cone))
    {
        if (a.type == sggk::SurfaceType::Cylinder)
        {
            e.o = {a.radius * c, a.radius * s, 0};
            e.d = {0, 0, 1};
        }
        else
            e.d = {std::tan(a.halfAngle) * c, std::tan(a.halfAngle) * s, 1};
        e.n0 = CanonicalNormal(a, e.value, 0);
    }
    else
    {
        e.circle = true;
        if (e.fixed == 1)
        {
            e.ex = {1, 0, 0};
            e.ey = {0, 1, 0};
            switch (a.type)
            {
            case sggk::SurfaceType::Cylinder:
                e.radius = a.radius;
                e.o = {0, 0, e.value};
                e.nx = e.ex;
                e.ny = e.ey;
                break;
            case sggk::SurfaceType::Cone:
                e.radius = e.value * std::tan(a.halfAngle);
                e.o = {0, 0, e.value};
                e.n0 = {0, 0, -std::sin(a.halfAngle)};
                e.nx = e.ex * std::cos(a.halfAngle);
                e.ny = e.ey * std::cos(a.halfAngle);
                break;
            case sggk::SurfaceType::Sphere:
                e.radius = a.radius * s;
                e.o = {0, 0, -a.radius * c};
                e.n0 = {0, 0, -c};
                e.nx = e.ex * s;
                e.ny = e.ey * s;
                break;
            case sggk::SurfaceType::Torus:
                e.radius = a.majorRadius + a.radius * c;
                e.o = {0, 0, a.radius * s};
                e.n0 = {0, 0, s};
                e.nx = e.ex * c;
                e.ny = e.ey * c;
                break;
            default:
                break;
            }
        }
        else if (a.type == sggk::SurfaceType::Sphere)
        {
            e.radius = a.radius;
            e.ex = {0, 0, -1};
            e.ey = {c, s, 0};
            e.nx = e.ex;
            e.ny = e.ey;
        }
        else
        {
            e.radius = a.radius;
            e.o = {a.majorRadius * c, a.majorRadius * s, 0};
            e.ex = {c, s, 0};
            e.ey = {0, 0, 1};
            e.nx = e.ex;
            e.ny = e.ey;
        }
        if (!(e.radius > 0))
            throw Failure(
                Status::Degenerate,
                "Boundary circle is collapsed or outside the regular positive-radial chart");
    }
    return e;
}
double EdgeCoordinateRaw(const EdgeGeometry &e, V p)
{
    p = p - e.o;
    return e.circle ? std::atan2(Dot(p, e.ey), Dot(p, e.ex)) : Dot(p, e.d) / Dot(e.d, e.d);
}
double EdgeCoordinate(const EdgeGeometry &e, V p, double reference)
{
    double raw = EdgeCoordinateRaw(e, p);
    return e.circle ? Unwrap(raw, reference) : raw;
}
std::array<EdgeGeometry, 4> BoundaryMaps(const sggk::BSplineSurface &s, const AnalyticPatch &p)
{
    auto f = FrameOf(p.parameters);
    std::array<EdgeGeometry, 4> edges;
    for (unsigned k = 0; k < 4; ++k)
    {
        auto e = MakeEdge(p, k);
        auto point = [&](double t) {
            return f.local(Get(s.CalcPoint(k < 2 ? (k == 0 ? s.MinParamU() : s.MaxParamU()) : t,
                                           k < 2 ? t : (k == 2 ? s.MinParamV() : s.MaxParamV()))));
        };
        double a = k < 2 ? s.MinParamV() : s.MinParamU(), b = k < 2 ? s.MaxParamV() : s.MaxParamU();
        e.g0 = EdgeCoordinate(e, point(a), e.start);
        // Track the input order, rather than selecting a shortest arc at the endpoint.
        double last = e.g0;
        auto params = SampleParameters(k < 2 ? s.KnotsV() : s.KnotsU(), a, b, 8);
        for (double t : params)
            last = EdgeCoordinate(e, point(t), last);
        e.g1 = last;
        if (std::abs(e.g1 - e.g0) < 1e-12)
            throw Failure(Status::Degenerate, "Boundary correspondence has zero extent");
        edges[k] = e;
    }
    return edges;
}
EdgePoint MapBoundaryPoint(const AnalyticPatch &p, const EdgeGeometry &e, V world, double fraction)
{
    auto f = FrameOf(p.parameters);
    double g = EdgeCoordinate(e, f.local(world), e.g0 + (e.g1 - e.g0) * fraction);
    double s = e.start + (e.end - e.start) * (g - e.g0) / (e.g1 - e.g0);
    V point, tangent, normal;
    if (e.circle)
    {
        point = e.o + (e.ex * std::cos(s) + e.ey * std::sin(s)) * e.radius;
        tangent = e.ex * (-std::sin(s)) + e.ey * std::cos(s);
        normal = e.n0 + e.nx * std::cos(s) + e.ny * std::sin(s);
    }
    else
    {
        point = e.o + e.d * s;
        tangent = e.d;
        normal = e.n0;
    }
    return {f.global(point), Unit(f.vector(tangent)) * (e.end > e.start ? 1 : -1),
            Unit(f.vector(normal)) * double(p.normalSign), s};
}
V ProjectInteriorPoint(const AnalyticPatch &p, V world)
{
    auto f = FrameOf(p.parameters);
    V a = f.local(world);
    auto uv = MotherParameters(p.parameters, a);
    if (p.parameters.type != sggk::SurfaceType::Plane)
        uv[0] = Unwrap(uv[0], (p.domain[0] + p.domain[1]) * 0.5);
    if (p.parameters.type == sggk::SurfaceType::Torus)
        uv[1] = Unwrap(uv[1], (p.domain[2] + p.domain[3]) * 0.5);
    V best{};
    double error = INFINITY;
    auto accept = [&](V q) {
        double d = Dot(q - a, q - a);
        if (d < error)
        {
            error = d;
            best = q;
        }
    };
    if (uv[0] >= p.domain[0] && uv[0] <= p.domain[1] && uv[1] >= p.domain[2] &&
        uv[1] <= p.domain[3])
        return f.global(CanonicalPoint(p.parameters, uv[0], uv[1]));
    for (unsigned k = 0; k < 4; ++k)
    {
        auto e = MakeEdge(p, k);
        double t = EdgeCoordinate(e, a, (e.start + e.end) * 0.5);
        t = std::clamp(t, std::min(e.start, e.end), std::max(e.start, e.end));
        auto q = [&](double s) {
            return e.circle ? e.o + (e.ex * std::cos(s) + e.ey * std::sin(s)) * e.radius
                            : e.o + e.d * s;
        };
        accept(q(t));
        accept(q(e.start));
        accept(q(e.end));
    }
    return f.global(best);
}
V ProjectInteriorPointLocal(const AnalyticPatch &p, V a, const std::array<EdgeGeometry, 4> &edges,
                            double coneTangent)
{
    auto uv = p.parameters.type == sggk::SurfaceType::Cone
                  ? std::array<double, 2>{std::atan2(a.y, a.x),
                                          std::max(0.0, (std::hypot(a.x, a.y) * coneTangent + a.z) /
                                                            (1 + coneTangent * coneTangent))}
                  : MotherParameters(p.parameters, a);
    if (p.parameters.type != sggk::SurfaceType::Plane)
        uv[0] = Unwrap(uv[0], (p.domain[0] + p.domain[1]) * 0.5);
    if (p.parameters.type == sggk::SurfaceType::Torus)
        uv[1] = Unwrap(uv[1], (p.domain[2] + p.domain[3]) * 0.5);
    V best{};
    double error = INFINITY;
    auto accept = [&](V q) {
        double d = Dot(q - a, q - a);
        if (d < error)
        {
            error = d;
            best = q;
        }
    };
    if (uv[0] >= p.domain[0] && uv[0] <= p.domain[1] && uv[1] >= p.domain[2] &&
        uv[1] <= p.domain[3])
        return p.parameters.type == sggk::SurfaceType::Cone
                   ? V{uv[1] * coneTangent * std::cos(uv[0]), uv[1] * coneTangent * std::sin(uv[0]),
                       uv[1]}
                   : CanonicalPoint(p.parameters, uv[0], uv[1]);
    for (unsigned k = 0; k < 4; ++k)
    {
        const auto &e = edges[k];
        double t = EdgeCoordinate(e, a, (e.start + e.end) * 0.5);
        t = std::clamp(t, std::min(e.start, e.end), std::max(e.start, e.end));
        auto q = [&](double s) {
            return e.circle ? e.o + (e.ex * std::cos(s) + e.ey * std::sin(s)) * e.radius
                            : e.o + e.d * s;
        };
        accept(q(t));
        accept(q(e.start));
        accept(q(e.end));
    }
    return best;
}
EdgePoint MapBoundaryPointAndTangent(const Frame3 &f, const EdgeGeometry &e, double g)
{
    double s = e.start + (e.end - e.start) * (g - e.g0) / (e.g1 - e.g0);
    V point, tangent;
    if (e.circle)
    {
        double c = std::cos(s), sn = std::sin(s);
        point = e.o + (e.ex * c + e.ey * sn) * e.radius;
        tangent = e.ex * (-sn) + e.ey * c;
    }
    else
    {
        point = e.o + e.d * s;
        tangent = e.d;
    }
    return {f.global(point), Unit(f.vector(tangent)) * (e.end > e.start ? 1 : -1), {}, s};
}
void MaterializePatch(AnalyticPatch &p)
{
    MakeNative(p);
    for (unsigned i = 0; i < 4; ++i)
    {
        unsigned dim = p.mapping.swappedUV ? 1 : 0;
        std::array<double, 2> uv{};
        uv[dim] = p.domain[2 * dim + (bool(i % 2) != p.mapping.reversedU ? 1 : 0)];
        uv[1 - dim] = p.domain[2 * (1 - dim) + (bool(i / 2) != p.mapping.reversedV ? 1 : 0)];
        p.mapping.analyticCorners[i] = sggk::Point2D(uv[0], uv[1]);
    }
}

double Angle(V a, V b)
{
    return std::atan2(Norm(Cross(a, b)), Dot(a, b));
}
static double Median(std::vector<double> a)
{
    std::sort(a.begin(), a.end());
    return a[a.size() / 2];
}
ClassificationSamples BuildClassificationSamples(const sggk::BSplineSurface &s,
                                                 ClassificationDiagnostics *diagnostics,
                                                 SurfaceSampleContext *shared)
{
    SurfaceSampleContext local(s);
    auto &samples = shared ? *shared : local;
    ClassificationSamples cache;
    cache.sourceDomain = {s.MinParamU(), s.MaxParamU(), s.MinParamV(), s.MaxParamV()};
    unsigned pointCalls = 0;
    auto point = [&](double u, double v) {
        ++pointCalls;
        return samples.At(u, v,sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
    };
    for (unsigned k = 0; k < 4; ++k)
    {
        auto ts =
            SampleParameters(k < 2 ? s.KnotsV() : s.KnotsU(), k < 2 ? s.MinParamV() : s.MinParamU(),
                             k < 2 ? s.MaxParamV() : s.MaxParamU(), 8);
        cache.boundaryParameters[k] = ts;
        for (double t : ts)
            cache.boundaries[k].push_back(
                point(k < 2 ? (k == 0 ? s.MinParamU() : s.MaxParamU()) : t,
                      k < 2 ? t : (k == 2 ? s.MinParamV() : s.MaxParamV())));
    }
    const auto &greville = samples.Greville(false);
    for (size_t i : samples.Scores())
        cache.scorePoints.push_back(greville[i].point);
    double u = (s.MinParamU() + s.MaxParamU()) * .5, v = (s.MinParamV() + s.MaxParamV()) * .5;
    const auto &d = samples.At(u, v);
    cache.centerPoint = point(u, v);
    cache.centerDu = d.du;
    cache.centerDv = d.dv;
    if (diagnostics)
    {
        ++diagnostics->sampleCacheBuilds;
        diagnostics->cachedPointEvaluations += pointCalls;
        ++diagnostics->cachedDerivativeEvaluations;
    }
    return cache;
}
AnalyticPatch RecoverApproxPatch(const sggk::BSplineSurface &s, const AnalyticParameters &a)
{
    return RecoverApproxPatch(BuildClassificationSamples(s), a);
}
PatchRecoveryTrace TraceApproximateBoundaries(const ClassificationSamples &cache,
                                             const AnalyticParameters &a)
{
    PatchRecoveryTrace result;
    auto &paths = result.paths;
    auto f = FrameOf(a);
    auto corner = MotherParameters(a, f.local(cache.boundaries[0].front()));
    auto trace = [&](unsigned k, std::array<double, 2> last) {
        for (V point : cache.boundaries[k])
        {
            auto q = MotherParameters(a, f.local(point));
            if (a.type != sggk::SurfaceType::Plane)
                q[0] = Unwrap(q[0], last[0]);
            if (a.type == sggk::SurfaceType::Torus)
                q[1] = Unwrap(q[1], last[1]);
            paths[k].push_back(q);
            last = q;
        }
    };
    trace(0, corner);
    trace(2, corner);
    trace(1, paths[2].back());
    trace(3, paths[0].back());
    return result;
}
bool SameMotherParameters(const AnalyticParameters &a, const AnalyticParameters &b)
{
    auto same = [](V x, V y) { return x.x == y.x && x.y == y.y && x.z == y.z; };
    return a.type == b.type && same(Get(a.origin), Get(b.origin)) &&
           same(Get(a.axis), Get(b.axis)) && same(Get(a.xDirection), Get(b.xDirection)) &&
           a.radius == b.radius && a.majorRadius == b.majorRadius && a.halfAngle == b.halfAngle;
}
AnalyticPatch RecoverApproxPatch(const ClassificationSamples &cache, const AnalyticParameters &a,
                                 PatchRecoveryTrace *savedTrace)
{
    AnalyticPatch p;
    p.parameters = a;
    auto f = FrameOf(a);
    auto trace = TraceApproximateBoundaries(cache, a);
    const auto &paths = trace.paths;
    std::array<unsigned, 4> fixed{};
    std::array<double, 4> value{};
    for (unsigned k = 0; k < 4; ++k)
    {
        std::array<double, 2> range{};
        for (unsigned j = 0; j < 2; ++j)
        {
            double lo = INFINITY, hi = -INFINITY;
            for (auto q : paths[k])
            {
                lo = std::min(lo, q[j]);
                hi = std::max(hi, q[j]);
            }
            range[j] = hi - lo;
        }
        fixed[k] = range[0] < range[1] ? 0 : 1;
        std::vector<double> q;
        for (auto x : paths[k])
            q.push_back(x[fixed[k]]);
        value[k] = Median(q);
    }
    if (fixed[0] != fixed[1] || fixed[2] != fixed[3] || fixed[0] == fixed[2])
        throw Failure(Status::InconsistentBoundary,
                      "Approximate boundaries do not form a rectangular analytic chart");
    p.mapping.swappedUV = fixed[0] == 1;
    p.mapping.reversedU = value[1] < value[0];
    p.mapping.reversedV = value[3] < value[2];
    for (unsigned group = 0; group < 2; ++group)
    {
        unsigned j = fixed[2 * group];
        p.domain[2 * j] = std::min(value[2 * group], value[2 * group + 1]);
        p.domain[2 * j + 1] = std::max(value[2 * group], value[2 * group + 1]);
    }
    if (p.domain[1] - p.domain[0] < 1e-10 || p.domain[3] - p.domain[2] < 1e-10)
        throw Failure(Status::Degenerate, "Zero analytic parameter extent");
    if (a.type != sggk::SurfaceType::Plane && p.domain[1] - p.domain[0] > 2 * Pi + 1e-5)
        throw Failure(Status::Unsupported, "Multiple-turn patch");
    for (unsigned i = 0; i < 4; ++i)
    {
        unsigned uHigh = i % 2, vHigh = i / 2;
        std::array<double, 2> uv{};
        unsigned dim = p.mapping.swappedUV ? 1 : 0;
        uv[dim] = p.domain[2 * dim + (bool(uHigh) != p.mapping.reversedU ? 1 : 0)];
        uv[1 - dim] = p.domain[2 * (1 - dim) + (bool(vHigh) != p.mapping.reversedV ? 1 : 0)];
        p.mapping.analyticCorners[i] = sggk::Point2D(uv[0], uv[1]);
    }
    auto uv = MotherParameters(a, f.local(cache.centerPoint));
    p.normalSign =
        Dot(Cross(cache.centerDu, cache.centerDv), f.vector(CanonicalNormal(a, uv[0], uv[1]))) >= 0
            ? 1
            : -1;
    MakeNative(p);
    if (savedTrace)
        *savedTrace = std::move(trace);
    return p;
}
void CheckApproximateInput(const sggk::BSplineSurface &s, const ValidationOptions &o)
{
    if (!(o.distanceTolerance > 0 && std::isfinite(o.distanceTolerance)) ||
        !(o.boundaryTangentToleranceDeg > 0 && o.boundaryTangentToleranceDeg < 90) ||
        !(o.boundaryNormalToleranceDeg > 0 && o.boundaryNormalToleranceDeg < 90) ||
        o.maxSubdivisionDepth > 30 || !o.maxValidationPieces)
        throw Failure(Status::InvalidInput, "Invalid conversion tolerances or budgets");
    if (s.IsRational())
        throw Failure(Status::Unsupported, "ConvertApproximate requires a non-rational B-spline");
    ValidateSurface(s);
    if (!s.DegreeU() || !s.DegreeV())
        throw Failure(Status::Degenerate, "Zero degree surface direction");
    auto knots = [](const sggk::RealArray &k, const sggk::UIntArray &m, unsigned degree,
                    size_t points) {
        if (k.size() < 2 || k.size() != m.size() || degree >= points)
            throw Failure(Status::InvalidInput, "Invalid knot/control dimensions");
        size_t count = 0;
        for (size_t i = 0; i < k.size(); ++i)
        {
            if (!std::isfinite(k[i]) || !m[i] || m[i] > degree + 1 || (i && k[i] <= k[i - 1]))
                throw Failure(Status::InvalidInput, "Invalid knot sequence");
            count += m[i];
        }
        if (count != points + degree + 1)
            throw Failure(Status::InvalidInput, "Invalid full knot count");
    };
    knots(s.KnotsU(), s.MultsU(), s.DegreeU(), s.ControlPoints().size());
    knots(s.KnotsV(), s.MultsV(), s.DegreeV(), s.ControlPoints().front().size());
    // Bernstein extraction and finite differences need a continuous input surface.
    for (size_t i = 1; i + 1 < s.MultsU().size(); ++i)
        if (s.MultsU()[i] > s.DegreeU())
            throw Failure(Status::Unsupported, "Discontinuous U knot");
    for (size_t i = 1; i + 1 < s.MultsV().size(); ++i)
        if (s.MultsV()[i] > s.DegreeV())
            throw Failure(Status::Unsupported, "Discontinuous V knot");
}
} // namespace bspline2analytic::detail
