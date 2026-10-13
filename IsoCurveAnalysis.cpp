#include "ApproximateInternal.h"
#include <numeric>
namespace bspline2analytic::detail
{
static std::array<V, 3> EigenFrame(double a[3][3], std::array<double, 3> &eigen)
{
    double q[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (unsigned iter = 0; iter < 40; ++iter)
    {
        unsigned p = 0, r = 1;
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = i + 1; j < 3; ++j)
                if (std::abs(a[i][j]) > std::abs(a[p][r]))
                {
                    p = i;
                    r = j;
                }
        if (std::abs(a[p][r]) < 1e-16 * (std::abs(a[0][0]) + std::abs(a[1][1]) + std::abs(a[2][2])))
            break;
        double phi = 0.5 * std::atan2(2 * a[p][r], a[r][r] - a[p][p]), c = std::cos(phi),
               s = std::sin(phi);
        for (unsigned i = 0; i < 3; ++i)
        {
            double x = a[i][p], y = a[i][r];
            a[i][p] = c * x - s * y;
            a[i][r] = s * x + c * y;
        }
        for (unsigned j = 0; j < 3; ++j)
        {
            double x = a[p][j], y = a[r][j];
            a[p][j] = c * x - s * y;
            a[r][j] = s * x + c * y;
        }
        for (unsigned i = 0; i < 3; ++i)
        {
            double x = q[i][p], y = q[i][r];
            q[i][p] = c * x - s * y;
            q[i][r] = s * x + c * y;
        }
    }
    std::array<unsigned, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(),
              [&](unsigned i, unsigned j) { return a[i][i] < a[j][j]; });
    std::array<V, 3> out;
    for (unsigned i = 0; i < 3; ++i)
    {
        unsigned k = order[i];
        out[i] = Unit({q[0][k], q[1][k], q[2][k]});
        eigen[i] = a[k][k];
    }
    return out;
}
struct CurveFits
{
    BoundaryInfo line, circle;
    bool hasLine = false, hasCircle = false;
};
static CurveFits FitCurve(const sggk::Curve3DPtr &curve, BoundaryId id, double tolerance,
                          unsigned maxSamples)
{
    auto bc = std::dynamic_pointer_cast<sggk::BSplineCurve3D>(curve);
    if (!bc)
        throw Failure(Status::Unsupported, "Expected spline isocurve");
    auto ts = SampleParameters(*bc, 8);
    if (ts.size() > maxSamples)
        throw Failure(Status::VerificationFailed, "Isocurve sample budget exceeded");
    std::vector<V> points;
    for (double t : ts)
        points.push_back(Get(curve->CalcPoint(t)));
    std::vector<double> weights(points.size(), 0);
    for (size_t i = 1; i < points.size(); ++i)
    {
        double length = Norm(points[i] - points[i - 1]);
        // Midpoint chord quadrature is representation-independent to sampling accuracy.
        V mid = Get(curve->CalcPoint((ts[i] + ts[i - 1]) * 0.5));
        length = std::max(length, Norm(mid - points[i - 1]) + Norm(points[i] - mid));
        weights[i - 1] += length * 0.5;
        weights[i] += length * 0.5;
    }
    double length = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (!(length > 1e-12))
        throw Failure(Status::Degenerate, "Collapsed isocurve");
    for (auto &w : weights)
        w /= length;
    V center{};
    for (size_t i = 0; i < points.size(); ++i)
        center = center + points[i] * weights[i];
    double scale = 0;
    for (V x : points)
        scale = std::max(scale, Norm(x - center));
    std::vector<V> a;
    for (V x : points)
        a.push_back((x - center) / scale);
    double cov[3][3]{};
    for (size_t i = 0; i < a.size(); ++i)
    {
        double z[3] = {a[i].x, a[i].y, a[i].z};
        for (unsigned j = 0; j < 3; ++j)
            for (unsigned k = 0; k < 3; ++k)
                cov[j][k] += weights[i] * z[j] * z[k];
    }
    std::array<double, 3> eigen;
    auto frame = EigenFrame(cov, eigen);
    CurveFits out;
    out.line.id = id;
    out.line.curve = curve;
    out.line.kind = BoundaryKind::Line;
    V direction = frame[2];
    if (Dot(direction, points.back() - points.front()) < 0)
        direction = direction * (-1);
    double lineError = 0;
    for (V x : points)
        lineError = std::max(lineError, Norm(Cross(x - center, direction)));
    out.line.line = {Point(center), Direction(direction), Dot(points.front() - center, direction),
                     Dot(points.back() - center, direction)};
    out.line.errorBound = lineError;
    out.line.message = "Sampled approximate fit; not a certified curve bound";
    out.hasLine = lineError <= tolerance;
    if (eigen[1] < 1e-10 * std::max(eigen[2], 1e-30))
        return out;
    V normal = frame[0], ex = frame[2], ey = Cross(normal, ex);
    std::vector<std::array<double, 3>> rows;
    std::vector<double> rhs;
    for (size_t i = 0; i < a.size(); ++i)
    {
        double x = Dot(a[i], ex), y = Dot(a[i], ey), w = std::sqrt(weights[i]);
        rows.push_back({w * x, w * y, w});
        rhs.push_back(-w * (x * x + y * y));
    }
    std::array<double, 3> sol;
    if (!LeastSquares3(rows, rhs, 1e-11, sol))
        return out;
    V o = ex * (-sol[0] * 0.5) + ey * (-sol[1] * 0.5);
    double r2 = (sol[0] * sol[0] + sol[1] * sol[1]) * 0.25 - sol[2];
    if (!(r2 > 0))
        return out;
    double radius = std::sqrt(r2);
    auto residual = [&](V c, V n, double r) {
        std::vector<double> f;
        for (size_t i = 0; i < a.size(); ++i)
        {
            V d = a[i] - c;
            double h = Dot(d, n), rho = Norm(d - n * h), w = std::sqrt(weights[i]);
            f.push_back(w * h);
            f.push_back(w * (rho - r));
        }
        return f;
    };
    for (unsigned iteration = 0; iteration < 8; ++iteration)
    {
        V nx = Perpendicular(normal), ny = Cross(normal, nx);
        auto f = residual(o, normal, radius);
        std::vector<std::vector<double>> J(f.size() + 6, std::vector<double>(6));
        const double h = 1e-5;
        auto moved = [&](unsigned j, double step) {
            V c = o, n = normal;
            double r = radius;
            if (j < 3)
            {
                V d{};
                if (j == 0)
                    d.x = step;
                if (j == 1)
                    d.y = step;
                if (j == 2)
                    d.z = step;
                c = c + d;
            }
            else if (j < 5)
                n = Unit(n + (j == 3 ? nx : ny) * step);
            else
                r += step;
            return residual(c, n, r);
        };
        for (unsigned j = 0; j < 6; ++j)
        {
            auto fp = moved(j, h), fm = moved(j, -h);
            for (size_t i = 0; i < f.size(); ++i)
                J[i][j] = (fp[i] - fm[i]) / (2 * h);
            J[f.size() + j][j] = 1e-5;
        }
        std::vector<double> b(f.size() + 6, 0), step;
        for (size_t i = 0; i < f.size(); ++i)
            b[i] = -f[i];
        if (!SolveLeastSquares(J, b, step))
            break;
        V trial = o + V{step[0], step[1], step[2]}, nn = Unit(normal + nx * step[3] + ny * step[4]);
        double rr = radius + step[5];
        if (!(rr > 0))
            break;
        auto next = residual(trial, nn, rr);
        double old = 0, now = 0;
        for (double x : f)
            old += x * x;
        for (double x : next)
            now += x * x;
        if (now >= old)
            break;
        o = trial;
        normal = nn;
        radius = rr;
        double norm = 0;
        for (double x : step)
            norm += x * x;
        if (norm < 1e-24)
            break;
    }
    ex = Unit(ex - normal * Dot(ex, normal));
    ey = Cross(normal, ex);
    V worldCenter = center + o * scale;
    double maxError = 0, first = 0, last = 0;
    for (size_t i = 0; i < points.size(); ++i)
    {
        V d = points[i] - worldCenter;
        double z = Dot(d, normal), rho = Norm(d - normal * z);
        maxError = std::max(maxError, std::hypot(z, rho - radius * scale));
        double angle = std::atan2(Dot(d, ey), Dot(d, ex));
        if (i == 0)
            first = last = angle;
        else
            last = Unwrap(angle, last);
    }
    out.circle.id = id;
    out.circle.curve = curve;
    out.circle.kind = BoundaryKind::Circle;
    out.circle.circle = {
        Point(worldCenter), Direction(normal), Direction(ex), radius * scale, first, last - first};
    out.circle.errorBound = maxError;
    out.circle.message = "Sampled geometric circle fit; not a certified curve bound";
    out.hasCircle = maxError <= tolerance && std::abs(last - first) > 1e-3;
    return out;
}
std::vector<AnalyticPatch> ClassifyApproximate(const sggk::BSplineSurface &s,
                                               const ConversionOptions &o)
{
    std::array<std::vector<CurveFits>, 2> groups;
    unsigned n = o.sampling.interiorIsoCurves;
    if (n > 5)
        throw Failure(Status::InvalidInput, "At most five internal isocurves per direction");
    double tolerance = o.validation.distanceTolerance * o.refinement.candidateToleranceFactor;
    for (unsigned g = 0; g < 2; ++g)
        for (unsigned j = 0; j < n + 2; ++j)
        {
            double a = g == 0 ? s.MinParamU() : s.MinParamV(),
                   b = g == 0 ? s.MaxParamU() : s.MaxParamV();
            double t = a + (b - a) * double(j) / (n + 1);
            auto curve = g == 0 ? s.CalcVCurve(t) : s.CalcUCurve(t);
            groups[g].push_back(FitCurve(curve,
                                         static_cast<BoundaryId>(2 * g + (j == n + 1 ? 1 : 0)),
                                         tolerance, o.sampling.maxSamples));
        }
    struct Scored
    {
        AnalyticPatch patch;
        double score;
        bool fromBoundary = false;
    };
    std::vector<Scored> best, boundarySeeds;
    std::vector<std::pair<unsigned, unsigned>> pairs{{0, n + 1}};
    if (n >= 2)
        pairs.push_back({1, n});
    if (n >= 1)
    {
        pairs.push_back({0, (n + 1) / 2});
        pairs.push_back({(n + 1) / 2, n + 1});
    }
    RecognitionOptions ro;
    ro.distanceTolerance = tolerance;
    ro.angularTolerance = 1e-3;
    for (unsigned kind0 = 0; kind0 < 2; ++kind0)
        for (unsigned kind1 = 0; kind1 < 2; ++kind1)
            for (auto pair0 : pairs)
                for (auto pair1 : pairs)
                {
                    std::array<BoundaryInfo, 4> edges;
                    bool valid = true;
                    for (unsigned g = 0; g < 2; ++g)
                    {
                        unsigned kind = g ? kind1 : kind0;
                        auto pair = g ? pair1 : pair0;
                        for (unsigned j = 0; j < 2; ++j)
                        {
                            const auto &fits = groups[g][j ? pair.second : pair.first];
                            valid &= kind ? fits.hasCircle : fits.hasLine;
                            edges[2 * g + j] = kind ? fits.circle : fits.line;
                        }
                    }
                    if (!valid)
                        continue;
                    for (const auto &c : Reconstruct(edges, ro))
                    {
                        try
                        {
                            auto patch = RecoverApproxPatch(s, c.parameters);
                            double score = 0;
                            for (unsigned i = 0; i < 5; ++i)
                                for (unsigned j = 0; j < 5; ++j)
                                    score = std::max(
                                        score,
                                        Distance(c.parameters,
                                                 s.CalcPoint(
                                                     s.MinParamU() +
                                                         (s.MaxParamU() - s.MinParamU()) * i / 4,
                                                     s.MinParamV() +
                                                         (s.MaxParamV() - s.MinParamV()) * j / 4)));
                            const bool fromBoundary =
                                pair0 == pairs.front() && pair1 == pairs.front();
                            if (fromBoundary)
                            {
                                auto seed = std::find_if(boundarySeeds.begin(), boundarySeeds.end(),
                                                         [&](const Scored &x) {
                                                             return x.patch.parameters.type ==
                                                                    c.parameters.type;
                                                         });
                                if (seed == boundarySeeds.end())
                                    boundarySeeds.push_back({patch, score, true});
                                else if (score < seed->score)
                                    *seed = {patch, score, true};
                            }
                            auto found =
                                std::find_if(best.begin(), best.end(), [&](const Scored &x) {
                                    return x.patch.parameters.type == c.parameters.type;
                                });
                            if (found == best.end())
                                best.push_back({patch, score, fromBoundary});
                            else if (score < found->score)
                                *found = {patch, score, fromBoundary};
                        }
                        catch (const Failure &)
                        {
                        }
                    }
                }
    // Refitted spherical meridians can have slightly different fitted centers,
    // also creating a torus with an almost zero major radius. Its marginally
    // smaller sampled score should not put this singular limit ahead of an
    // already competitive sphere. Retain it as a fallback; only full validation
    // can accept either candidate, and no provenance/type label is consulted.
    const bool hasSphere = std::any_of(best.begin(), best.end(), [&](const Scored &c) {
        return c.patch.parameters.type == sggk::SurfaceType::Sphere &&
               c.score < o.validation.distanceTolerance;
    });
    auto deferTorus = [&](const Scored &c) {
        return hasSphere && c.patch.parameters.type == sggk::SurfaceType::Torus &&
               std::abs(c.patch.parameters.majorRadius) <= o.validation.distanceTolerance;
    };
    std::sort(best.begin(), best.end(), [&](const Scored &a, const Scored &b) {
        if (deferTorus(a) != deferTorus(b))
            return !deferTorus(a);
        return a.score < b.score;
    });
    std::vector<AnalyticPatch> result;
    for (auto &b : best)
    {
        // A low mother-surface score on internal isocurves can hide a poor
        // finite boundary chart. Retain the four-boundary initialization rather
        // than permanently replacing it with that scored alternative.
        auto seed = std::find_if(boundarySeeds.begin(), boundarySeeds.end(), [&](const Scored &x) {
            return x.patch.parameters.type == b.patch.parameters.type;
        });
        if (!b.fromBoundary && seed != boundarySeeds.end())
            result.push_back(std::move(seed->patch));
        result.push_back(std::move(b.patch));
    }
    return result;
}
} // namespace bspline2analytic::detail
