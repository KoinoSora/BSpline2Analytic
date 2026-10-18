/**
 * @file RecognitionValidation.cpp
 * @brief 精确识别候选的母面与边界误差检查。
 * 使用采样快速检查及带浮点保护的 Bernstein 残差界，报告实际达到的验证等级。
 */
#include "RecognitionValidation.h"
#include "Numerics.h"

namespace bspline2analytic::detail
{
static double Guard(double magnitude)
{
    return 512 * std::numeric_limits<double>::epsilon() * magnitude;
}
static std::pair<double, double> RatioRange(const Poly &p, const Poly &w)
{
    double a = p.min() / w.min(), b = p.min() / w.max(), c = p.max() / w.min(),
           d = p.max() / w.max();
    return {std::min({a, b, c, d}), std::max({a, b, c, d})};
}
static std::pair<double, double> SquareRange(double lo, double hi)
{
    return {lo <= 0 && hi >= 0 ? 0 : std::min(lo * lo, hi * hi), std::max(lo * lo, hi * hi)};
}
static double Bound(const Homogeneous &h, const AnalyticParameters &p, double scale)
{
    double wm = h.w.min();
    if (!(wm > 0))
        return std::numeric_limits<double>::infinity();
    V a = Get(p.axis);
    Poly z = DotNumerator(h, a);
    if (p.type == sggk::SurfaceType::Plane)
        return scale * (z.maxAbs() + Guard(h.x.maxAbs() + h.y.maxAbs() + h.z.maxAbs())) / wm;
    Poly ww = h.w * h.w, norm = SquaredNorm(h);
    double r = p.radius / scale;
    if (p.type == sggk::SurfaceType::Sphere)
    {
        Poly rr = ww * (r * r), f = norm - rr;
        return scale * (f.maxAbs() + Guard(norm.maxAbs() + rr.maxAbs())) / (wm * wm * r);
    }
    // 先投影后平方，避免极浅圆锥中的灾难性相消。
    Poly x = h.x - z * a.x, y = h.y - z * a.y, t = h.z - z * a.z;
    Poly radial = x * x + y * y + t * t;
    if (p.type == sggk::SurfaceType::Cylinder)
    {
        Poly rr = ww * (r * r), f = radial - rr;
        return scale * (f.maxAbs() + Guard(radial.maxAbs() + rr.maxAbs())) / (wm * wm * r);
    }
    auto zr = RatioRange(z, h.w);
    double rhoMin = std::sqrt(std::max(0.0, radial.min())) / h.w.max();
    double rhoMax = std::sqrt(std::max(0.0, radial.max())) / wm;
    if (p.type == sggk::SurfaceType::Cone)
    {
        if (zr.first < 0)
            return std::numeric_limits<double>::infinity();
        double ca = std::cos(p.halfAngle), sa = std::sin(p.halfAngle);
        double divisor = rhoMin * ca + zr.first * sa;
        if (!(divisor > 0))
            return std::numeric_limits<double>::infinity();
        Poly rterm = radial * (ca * ca), zterm = (z * z) * (sa * sa), f = rterm - zterm;
        return scale * (f.maxAbs() + Guard(rterm.maxAbs() + zterm.maxAbs())) / (wm * wm * divisor);
    }
    if (p.type == sggk::SurfaceType::Torus)
    {
        double R = p.majorRadius / scale;
        Poly sum = norm + ww * (R * R - r * r), left = sum * sum,
             right = (ww * radial) * (4 * R * R), f = left - right;
        // 利用 F=(q_minus²-r²)(q_plus²-r²)，对另一个因子建立严格远离零的下界。
        auto sr = SquareRange(rhoMin + R, rhoMax + R), sz = SquareRange(zr.first, zr.second);
        double lo = sr.first + sz.first - r * r, hi = sr.second + sz.second - r * r;
        double factor = lo > 0 ? lo : (hi < 0 ? -hi : 0);
        if (!(factor > 0))
            return std::numeric_limits<double>::infinity();
        return scale * (f.maxAbs() + Guard(left.maxAbs() + right.maxAbs())) /
               (wm * wm * wm * wm * factor * r);
    }
    return std::numeric_limits<double>::infinity();
}
static double ValidatePiece(const Homogeneous &h, const AnalyticParameters &p, double scale,
                            const RecognitionOptions &o, unsigned depth, unsigned &count, const ExecutionDeadline *deadline)
{
    if (deadline) deadline->Check("exact_surface_bound");
    if (++count > o.maxValidationPieces)
        throw Failure(Status::VerificationFailed, "Surface bound subdivision budget exceeded");
    double bound = Bound(h, p, scale);
    if (bound <= o.distanceTolerance)
        return bound;
    if (depth >= o.maxSubdivisionDepth)
        throw Failure(Status::VerificationFailed,
                      "Unable to establish the requested surface distance bound");
    auto children = Split(h, (depth % 2) == 0);
    return std::max(ValidatePiece(children.first, p, scale, o, depth + 1, count, deadline),
                    ValidatePiece(children.second, p, scale, o, depth + 1, count, deadline));
}
ValidationReport ValidateCandidate(const sggk::BSplineSurface &s, const Candidate &c,
                                   const std::array<BoundaryInfo, 4> &b,
                                   const RecognitionOptions &o, SurfaceSampleContext *shared)
{
    SurfaceSampleContext local(s);
    auto &samples = shared ? *shared : local;
    ValidationReport report;
    for (const auto &edge : b)
    {
        report.maxBoundaryBound = std::max(report.maxBoundaryBound, edge.errorBound);
        auto bc = std::dynamic_pointer_cast<sggk::BSplineCurve3D>(edge.curve);
        for (double t : SampleParameters(*bc, o.samplesPerSpan))
        {
            unsigned id=static_cast<unsigned>(edge.id);
            V p=samples.At(id<2?(id?s.MaxParamU():s.MinParamU()):t,
                id<2?t:(id==2?s.MinParamV():s.MaxParamV()),sggk::DerivType::NONE,
                sggk::DerivType::NONE,false).point;
            if (Distance(c.parameters, Point(p)) > o.distanceTolerance)
                throw Failure(Status::VerificationFailed,
                              "Recovered analytic surface does not contain a boundary");
        }
    }
    for (const auto &sample : samples.Greville(false))
        {
            double d = Distance(c.parameters, Point(sample.point));
            if (!std::isfinite(d))
                throw Failure(Status::VerificationFailed, "Nonfinite surface distance");
            report.maxSampledDistance = std::max(report.maxSampledDistance, d);
        }
    report.level = VerificationLevel::Sampled;
    if (report.maxSampledDistance > o.distanceTolerance)
        throw Failure(Status::VerificationFailed,
                      "Interior samples disagree with the recovered surface");
    if (!o.requireBound)
    {
        report.message = "Sampled verification only; no global distance guarantee.";
        return report;
    }
    V origin = Get(c.parameters.origin);
    double scale = std::max(c.parameters.radius, std::abs(c.parameters.majorRadius));
    for (const auto &row : s.ControlPoints())
        for (const auto &p : row)
            scale = std::max(scale, Norm(Get(p) - origin));
    if (!(scale > 0))
        throw Failure(Status::Degenerate, "Zero validation scale");
    report.surfaceDistanceBound = 0;
    for (const auto &h : SurfacePieces(s, origin, scale))
        report.surfaceDistanceBound =
            std::max(report.surfaceDistanceBound,
                     ValidatePiece(h, c.parameters, scale, o, 0, report.piecesVisited, samples.deadline));
    report.level = VerificationLevel::BernsteinBound;
    report.message = "Bernstein bound to analytic mother surface with floating-point guard; not a "
                     "directed-rounding interval certificate or trimmed-patch Hausdorff bound.";
    return report;
}
} // namespace bspline2analytic::detail
