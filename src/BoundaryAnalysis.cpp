/**
 * @file BoundaryAnalysis.cpp
 * @brief 四条源边的直线、圆弧及退化分类。
 * 使用正权重控制凸包和误差界检查几何，有序遍历处理周期接缝，不仅依据端点判断圆弧。
 */
#include "BoundaryAnalysis.h"
#include "Numerics.h"

namespace bspline2analytic::detail
{
static double BoundPiece(const Homogeneous &h, const BoundaryInfo &info, double scale,
                         const RecognitionOptions &o, unsigned depth, unsigned &visited)
{
    if (++visited > o.maxValidationPieces)
        throw Failure(Status::VerificationFailed, "Boundary validation budget exceeded");
    double bound = std::numeric_limits<double>::infinity(), wm = h.w.min();
    if (wm > 0)
    {
        if (info.kind == BoundaryKind::Line)
        {
            // 有理 Bernstein 权重全为正时，控制凸包可给出到直线的距离上界。
            bound = 0;
            V d = Get(info.line.direction);
            for (std::size_t i = 0; i < h.w.c.size(); ++i)
            {
                V p{h.x.c[i], h.y.c[i], h.z.c[i]};
                bound = std::max(bound, Norm(Cross(p / h.w.c[i], d)));
            }
            bound += Roundoff(1);
        }
        else
        {
            Poly z = DotNumerator(h, Get(info.circle.normal));
            Poly norm2 = SquaredNorm(h), ww = h.w * h.w;
            // 先投影至圆平面，再估计径向距离误差。
            Poly radial = norm2 - z * z - ww;
            double normalError = (z.maxAbs() + Roundoff(z.maxAbs())) / wm;
            double radialError =
                (radial.maxAbs() + Roundoff(norm2.maxAbs() + ww.maxAbs())) / (wm * wm);
            bound = std::hypot(normalError, radialError);
        }
        bound *= scale;
        if (bound <= o.distanceTolerance)
            return bound;
    }
    if (depth >= o.maxSubdivisionDepth)
        throw Failure(Status::VerificationFailed, "Boundary does not satisfy distance bound");
    auto children = Split(h, true);
    return std::max(BoundPiece(children.first, info, scale, o, depth + 1, visited),
                    BoundPiece(children.second, info, scale, o, depth + 1, visited));
}

BoundaryInfo AnalyzeCurve(const sggk::Curve3DPtr &original, BoundaryId id,
                          const RecognitionOptions &o)
{
    if (!original)
        throw Failure(Status::InvalidInput, "Null isoparametric boundary");
    auto curve = std::dynamic_pointer_cast<sggk::BSplineCurve3D>(original);
    if (!curve)
        curve = original->ToBSpline();
    if (!curve)
        throw Failure(Status::Unsupported, "Boundary cannot be represented as an exact B-spline");
    const auto &pts = curve->ControlPoints();
    if (pts.size() < 2)
        throw Failure(Status::Degenerate, "Collapsed boundary");
    V origin = Get(pts.front()), farthest = origin;
    double scale = 0, world = Norm(origin);
    for (const auto &p : pts)
    {
        V q = Get(p);
        if (!Finite(q))
            throw Failure(Status::InvalidInput, "Nonfinite boundary control point");
        double r = Norm(q - origin);
        world = std::max(world, Norm(q));
        if (r > scale)
        {
            scale = r;
            farthest = q;
        }
    }
    double numeric = Roundoff(world + scale);
    if (scale <= numeric)
        throw Failure(Status::Degenerate, "Boundary extent is below numerical resolution");
    V x = Unit(farthest - origin), perp{};
    double height = 0;
    for (const auto &p : pts)
    {
        V delta = Get(p) - origin;
        V v = delta - x * Dot(delta, x);
        double h = Norm(v);
        if (h > height)
        {
            height = h;
            perp = v;
        }
    }
    BoundaryInfo out;
    out.id = id;
    out.curve = curve;
    if (height <= numeric)
    {
        out.kind = BoundaryKind::Line;
        out.line.origin = Point(origin);
        out.line.direction = Direction(x);
    }
    else
    {
        V normal = Unit(Cross(x, perp)), y = Cross(normal, x);
        double planeError = 0;
        for (const auto &p : pts)
            planeError = std::max(planeError, std::abs(Dot(Get(p) - origin, normal)));
        if (planeError > numeric)
            throw Failure(Status::InconsistentBoundary, "Boundary control points are not planar");
        std::vector<std::array<double, 3>> rows;
        std::vector<double> rhs;
        for (const auto &h : CurvePieces(*curve, origin, scale))
        {
            Poly px = DotNumerator(h, x), py = DotNumerator(h, y);
            Poly xx = px * px + py * py, xw = px * h.w, yw = py * h.w, ww = h.w * h.w;
            for (std::size_t i = 0; i < xx.c.size(); ++i)
            {
                rows.push_back({xw.c[i], yw.c[i], ww.c[i]});
                rhs.push_back(-xx.c[i]);
            }
        }
        std::array<double, 3> solution{};
        if (!LeastSquares3(rows, rhs, o.rankTolerance, solution))
            throw Failure(Status::IllConditioned, "Circle recovery is rank deficient");
        double cx = -solution[0] * 0.5, cy = -solution[1] * 0.5,
               r2 = cx * cx + cy * cy - solution[2];
        if (!(r2 > 0) || !std::isfinite(r2))
            throw Failure(Status::InconsistentBoundary, "No finite positive circle radius");
        out.kind = BoundaryKind::Circle;
        out.circle.center = Point(origin + (x * cx + y * cy) * scale);
        out.circle.normal = Direction(normal);
        out.circle.xDirection = Direction(x);
        out.circle.radius = std::sqrt(r2) * scale;
    }
    V center = out.kind == BoundaryKind::Line ? Get(out.line.origin) : Get(out.circle.center);
    double boundScale = out.kind == BoundaryKind::Line ? scale : out.circle.radius;
    // 递归求界前，先快速拒绝不满足圆几何的圆锥曲线。
    for (double t : SampleParameters(*curve, o.samplesPerSpan))
    {
        V delta = Get(curve->CalcPoint(t)) - center;
        double error;
        if (out.kind == BoundaryKind::Line)
            error = Norm(Cross(delta, Get(out.line.direction)));
        else
        {
            V n = Get(out.circle.normal);
            double z = Dot(delta, n);
            error = std::hypot(z, Norm(delta - n * z) - out.circle.radius);
        }
        if (error > o.distanceTolerance)
            throw Failure(Status::InconsistentBoundary,
                          "Boundary samples disagree with the recovered line or circle");
    }
    unsigned visited = 0;
    out.errorBound = 0;
    for (const auto &h : CurvePieces(*curve, center, boundScale))
        out.errorBound = std::max(out.errorBound, BoundPiece(h, out, boundScale, o, 0, visited));

    // 沿源曲线依次遍历结点与接缝，不仅凭端点推断圆弧范围。
    const auto parameters = SampleParameters(*curve, std::max(8u, o.samplesPerSpan));
    double first = 0, last = 0, sign = 0;
    for (std::size_t i = 0; i < parameters.size(); ++i)
    {
        V p = Get(curve->CalcPoint(parameters[i])) - center;
        double value;
        if (out.kind == BoundaryKind::Line)
            value = Dot(p, Get(out.line.direction));
        else
        {
            V ex = Get(out.circle.xDirection), ey = Cross(Get(out.circle.normal), ex);
            value = std::atan2(Dot(p, ey), Dot(p, ex));
            if (i)
                value = Unwrap(value, last);
        }
        if (i == 0)
            first = value;
        else
        {
            double delta = value - last;
            const double traversalTolerance =
                out.kind == BoundaryKind::Line ? o.distanceTolerance : o.angularTolerance;
            if (std::abs(delta) > traversalTolerance)
            {
                if (sign && sign * delta < 0)
                    throw Failure(Status::InconsistentBoundary,
                                  "Boundary retraces its analytic support");
                sign = std::copysign(1.0, delta);
            }
        }
        last = value;
    }
    if (out.kind == BoundaryKind::Line)
    {
        out.line.start = first;
        out.line.end = last;
    }
    else
    {
        out.circle.startAngle = first;
        out.circle.sweepAngle = last - first;
    }
    return out;
}
} // namespace bspline2analytic::detail

namespace bspline2analytic
{
BoundaryAnalysisResult AnalyzeBoundaries(const sggk::BSplineSurface &s, const RecognitionOptions &o)
{
    BoundaryAnalysisResult r;
    try
    {
        detail::ValidateOptions(o);
        detail::ValidateSurface(s);
        const std::array<sggk::Curve3DPtr, 4> curves{
            s.CalcVCurve(s.MinParamU()), s.CalcVCurve(s.MaxParamU()), s.CalcUCurve(s.MinParamV()),
            s.CalcUCurve(s.MaxParamV())};
        for (unsigned i = 0; i < 4; ++i)
        {
            try {r.boundaries[i] = detail::AnalyzeCurve(curves[i], static_cast<BoundaryId>(i), o);}
            catch (const detail::Failure &e)
            {
                r.boundaries[i].id=static_cast<BoundaryId>(i);r.boundaries[i].curve=curves[i];
                r.boundaries[i].kind=e.status==Status::Degenerate ? BoundaryKind::Degenerate : BoundaryKind::Unknown;
                r.boundaries[i].message=e.what();throw;
            }
        }
        r.status = Status::Success;
    }
    catch (const detail::Failure &e)
    {
        r.status = e.status;
        r.message = e.what();
        r.failureReason = e.status == Status::Degenerate
                              ? ClassificationFailureReason::DegenerateBoundary
                              : (e.status == Status::InvalidInput
                                     ? ClassificationFailureReason::None
                                     : ClassificationFailureReason::NonLineOrCircleBoundary);
    }
    catch (const std::exception &e)
    {
        r.status = Status::InvalidInput;
        r.message = e.what();
    }
    return r;
}
} // namespace bspline2analytic
