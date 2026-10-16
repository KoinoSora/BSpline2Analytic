#include "ApproximateInternal.h"
#include <numeric>
#include <sstream>
#include <chrono>
namespace bspline2analytic
{
const char *CircleFitStateName(CircleFitState state)
{
    switch (state)
    {
    case CircleFitState::NotAttempted: return "NotAttempted";
    case CircleFitState::Accepted: return "Accepted";
    case CircleFitState::NumericallyLinear: return "NumericallyLinear";
    case CircleFitState::RankDeficient: return "RankDeficient";
    case CircleFitState::NonFinite: return "NonFinite";
    case CircleFitState::InvalidRadius: return "InvalidRadius";
    case CircleFitState::UnresolvedExtent: return "UnresolvedExtent";
    case CircleFitState::ResidualExceeded: return "ResidualExceeded";
    }
    return "Unknown";
}
}
namespace bspline2analytic::detail
{
// Small, contiguous column-pivoted Householder QR. Reuses caller storage and
// keeps the rank check; normal equations would square the condition number.
static bool SolveCircleQR(std::vector<std::array<double, 6>> &a, std::vector<double> &b,
                          std::array<double, 6> &x, std::vector<double> &v)
{
    const size_t m = a.size();
    std::array<unsigned, 6> perm{0,1,2,3,4,5};
    double first = 0;
    v.resize(m);
    for (size_t k = 0; k < 6; ++k)
    {
        size_t pivot = k;
        double best = -1;
        for (size_t j = k; j < 6; ++j)
        {
            double norm = 0;
            for (size_t i = k; i < m; ++i) norm = std::hypot(norm, a[i][j]);
            if (norm > best) { best = norm; pivot = j; }
        }
        if (!k) first = best;
        if (!(best > std::max(1e-30, first * 1e-12))) return false;
        if (pivot != k)
        {
            for (auto &row : a) std::swap(row[k], row[pivot]);
            std::swap(perm[k], perm[pivot]);
        }
        double alpha = -std::copysign(best, a[k][k]);
        for (size_t i = k; i < m; ++i) v[i] = a[i][k];
        v[k] -= alpha;
        double norm = 0;
        for (size_t i = k; i < m; ++i) norm = std::hypot(norm, v[i]);
        for (size_t i = k; i < m; ++i) v[i] /= norm;
        for (size_t j = k; j < 6; ++j)
        {
            double dot = 0;
            for (size_t i = k; i < m; ++i) dot += v[i] * a[i][j];
            for (size_t i = k; i < m; ++i) a[i][j] -= 2 * v[i] * dot;
        }
        double dot = 0;
        for (size_t i = k; i < m; ++i) dot += v[i] * b[i];
        for (size_t i = k; i < m; ++i) b[i] -= 2 * v[i] * dot;
        a[k][k] = alpha;
    }
    std::array<double, 6> y{};
    for (size_t j = 6; j-- > 0;)
    {
        double value = b[j];
        for (size_t k = j + 1; k < 6; ++k) value -= a[j][k] * y[k];
        y[j] = value / a[j][j];
    }
    for (size_t j = 0; j < 6; ++j) x[perm[j]] = y[j];
    return std::all_of(x.begin(), x.end(), [](double z) { return std::isfinite(z); });
}
struct IsoMeasure
{
    using Clock = std::chrono::steady_clock;
    double *elapsed;
    Clock::time_point start = Clock::now();
    explicit IsoMeasure(double &value) : elapsed(&value) {}
    void Switch(double &value)
    {
        auto end = Clock::now();
        *elapsed += std::chrono::duration<double, std::milli>(end - start).count();
        elapsed = &value; start = end;
    }
    ~IsoMeasure() { *elapsed += std::chrono::duration<double, std::milli>(Clock::now() - start).count(); }
};
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
CurveFits FitApproximateIsoCurve(const sggk::Curve3DPtr &curve, BoundaryId id, double tolerance,
                                unsigned maxSamples, IsoCurveFitDiagnostic *diagnostic,
                                const std::function<V(double)> &sample, const IsoCurveFitOptions &options)
{
    IsoCurveFitDiagnostic local;
    auto &diag = diagnostic ? *diagnostic : local;
    auto direction = diag.fixedDirection;
    auto fixedParameter = diag.fixedParameter;
    auto boundary = diag.boundary;
    diag = {};
    diag.fixedDirection = direction;
    diag.fixedParameter = fixedParameter;
    diag.boundary = boundary;
    IsoMeasure phase(diag.samplingMs);
    diag.subdivisionsPerSpan = options.subdivisionsPerSpan;
    diag.edge = id;
    diag.candidateTolerance = tolerance;
    if (!(tolerance > 0 && std::isfinite(tolerance)) || !maxSamples)
        throw Failure(Status::InvalidInput, "Invalid isocurve fitting options");
    auto bc = std::dynamic_pointer_cast<sggk::BSplineCurve3D>(curve);
    if (!bc)
        throw Failure(Status::Unsupported, "Expected spline isocurve");
    if (!options.subdivisionsPerSpan || options.subdivisionsPerSpan > 8)
        throw Failure(Status::InvalidInput, "Isocurve subdivisions must be between 1 and 8");
    auto ts = SampleParameters(*bc, options.subdivisionsPerSpan);
    diag.sampleCount = static_cast<unsigned>(ts.size());
    if (ts.size() > maxSamples)
        throw Failure(Status::VerificationFailed, "Isocurve sample budget exceeded");
    std::vector<V> points;
    for (double t : ts)
        points.push_back(sample ? sample(t) : Get(curve->CalcPoint(t)));
    double coordinateUlp = 0;
    for (V p : points)
    {
        if (!Finite(p))
        {
            diag.circleState = CircleFitState::NonFinite;
            throw Failure(Status::IllConditioned, "Nonfinite evaluated isocurve point");
        }
        for (double x : {p.x, p.y, p.z})
            coordinateUlp = std::max(coordinateUlp, std::max(
                std::nextafter(x, INFINITY)-x, x-std::nextafter(x, -INFINITY)));
    }
    std::vector<double> weights(points.size(), 0);
    for (size_t i = 1; i < points.size(); ++i)
    {
        double length = Norm(points[i] - points[i - 1]);
        // Reuse the fitting samples: chord weights require no extra evaluations.
        weights[i - 1] += length * 0.5;
        weights[i] += length * 0.5;
    }
    phase.Switch(diag.algebraicMs);
    double length = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (!std::isfinite(length))
        throw Failure(Status::IllConditioned, "Nonfinite isocurve extent");
    if (!(length > 0))
        throw Failure(Status::Degenerate, "Collapsed isocurve");
    for (auto &w : weights)
        w /= length;
    // Accumulate local offsets, avoiding a sum of large world coordinates.
    V offset{}, anchor = points.front();
    for (size_t i = 0; i < points.size(); ++i)
        offset = offset + (points[i] - anchor) * weights[i];
    V center = anchor + offset;
    double scale = 0;
    for (V x : points)
        scale = std::max(scale, Norm(x - center));
    if (!(scale > 0 && std::isfinite(scale)) || !Finite(center))
        throw Failure(Status::IllConditioned, "Unresolved isocurve normalization");
    // Heuristic numerical-resolution estimate, not a geometric certificate.
    diag.numericFloor = 8 * coordinateUlp + 64 * std::numeric_limits<double>::epsilon() * scale;
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
    V lineDirection = frame[2];
    if (Dot(lineDirection, points.back() - points.front()) < 0)
        lineDirection = lineDirection * (-1);
    double lineError = 0;
    for (V x : points)
        lineError = std::max(lineError, Norm(Cross(x - center, lineDirection)));
    out.line.line = {Point(center), Direction(lineDirection), Dot(points.front() - center, lineDirection),
                     Dot(points.back() - center, lineDirection)};
    out.line.errorBound = lineError;
    out.line.message = "Sampled approximate fit; not a certified curve bound";
    diag.lineError = lineError;
    out.hasLine = diag.hasLine = std::isfinite(lineError) && lineError <= tolerance;
    diag.transverseRatio = std::max(0.0, eigen[1]) / std::max(eigen[2], 1e-30);
    // A small covariance ratio alone does not establish geometric straightness.
    // Skip circles only when an admissible line is numerically indistinguishable.
    if (out.hasLine && lineError <= diag.numericFloor)
    {
        diag.circleState = CircleFitState::NumericallyLinear;
        return out;
    }
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
    {
        diag.circleState = CircleFitState::RankDeficient;
        return out;
    }
    V o = ex * (-sol[0] * 0.5) + ey * (-sol[1] * 0.5);
    double r2 = (sol[0] * sol[0] + sol[1] * sol[1]) * 0.25 - sol[2];
    if (!(r2 > 0 && std::isfinite(r2)) || !Finite(o))
    {
        diag.circleState = std::isfinite(r2) && Finite(o) ?
            CircleFitState::InvalidRadius : CircleFitState::NonFinite;
        return out;
    }
    double radius = std::sqrt(r2);
    auto maxCircleError = [&](V c, V n, double r) {
        double error = 0;
        for (V point : a)
        {
            V d = point - c;
            double h = Dot(d, n);
            error = std::max(error, std::hypot(h, Norm(d - n*h) - r));
        }
        return error * scale;
    };
    diag.skippedGeometricFit = !options.forceGeometric &&
        maxCircleError(o, normal, radius) <= options.fastTolerance;
    phase.Switch(diag.geometricMs);
    if (!diag.skippedGeometricFit)
    {
        const size_t rowsCount = 2*a.size();
        std::vector<std::array<double, 6>> J(rowsCount + 6);
        std::vector<double> f(rowsCount), next(rowsCount), b(rowsCount + 6), workspace;
        std::vector<double> rootWeight;
        for (double w : weights) rootWeight.push_back(std::sqrt(w));
        for (unsigned iteration = 0; iteration < 8; ++iteration)
        {
            ++diag.geometricIterations;
            V nx = Perpendicular(normal), ny = Cross(normal, nx);
            bool regular = true;
            for (size_t i = 0; i < a.size(); ++i)
            {
                V d = a[i] - o;
                double h = Dot(d, normal), rho = Norm(d - normal*h), w = rootWeight[i];
                if (!(rho > 1e-30)) { regular = false; break; }
                V q = (d - normal*h) / rho;
                double dx = Dot(d, nx), dy = Dot(d, ny);
                f[2*i] = w*h; f[2*i+1] = w*(rho-radius);
                J[2*i] = {-w*normal.x, -w*normal.y, -w*normal.z, w*dx, w*dy, 0};
                J[2*i+1] = {-w*q.x, -w*q.y, -w*q.z, -w*h*dx/rho, -w*h*dy/rho, -w};
            }
            if (!regular) break;
            for (size_t i = 0; i < rowsCount; ++i) b[i] = -f[i];
            for (size_t j = 0; j < 6; ++j)
            {
                J[rowsCount+j].fill(0); J[rowsCount+j][j] = 1e-5;
                b[rowsCount+j] = 0;
            }
            std::array<double, 6> step{};
            if (!SolveCircleQR(J, b, step, workspace)) break;
            V trial = o + V{step[0], step[1], step[2]};
            V nn = Unit(normal + nx*step[3] + ny*step[4]);
            double rr = radius + step[5];
            if (!(rr > 0 && std::isfinite(rr)) || !Finite(trial) || !Finite(nn)) break;
            double old = 0, now = 0, maxError = 0;
            for (size_t i = 0; i < a.size(); ++i)
            {
                V d = a[i] - trial;
                double h = Dot(d, nn), dr = Norm(d-nn*h)-rr;
                next[2*i] = rootWeight[i]*h; next[2*i+1] = rootWeight[i]*dr;
                maxError = std::max(maxError, std::hypot(h, dr)*scale);
            }
            for (double x : f) old += x*x;
            for (double x : next) now += x*x;
            if (!std::isfinite(now) || now >= old) break;
            o = trial; normal = nn; radius = rr;
            double norm = 0;
            for (double x : step) norm += x*x;
            if (norm < 1e-24 || old-now <= old*1e-10 ||
                (!options.forceGeometric && maxError <= options.fastTolerance)) break;
        }
    }
    phase.Switch(diag.algebraicMs);
    ex = Unit(ex - normal * Dot(ex, normal));
    ey = Cross(normal, ex);
    V worldCenter = center + o * scale;
    if (!Finite(worldCenter) || !std::isfinite(radius * scale))
    {
        diag.circleState = CircleFitState::NonFinite;
        return out;
    }
    double maxError = 0, first = 0, last = 0;
    for (size_t i = 0; i < points.size(); ++i)
    {
        V d = points[i] - worldCenter;
        double z = Dot(d, normal), rho = Norm(d - normal * z);
        double error = std::hypot(z, rho - radius * scale);
        double angle = std::atan2(Dot(d, ey), Dot(d, ex));
        if (!std::isfinite(error) || !std::isfinite(angle))
        {
            diag.circleState = CircleFitState::NonFinite;
            return out;
        }
        maxError = std::max(maxError, error);
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
    diag.circleError = maxError;
    diag.radius = radius * scale;
    diag.sweepAngle = last - first;
    double sine = std::sin(std::abs(diag.sweepAngle) * .25);
    diag.sagitta = diag.radius * (2 * sine * sine);
    double extent = diag.radius * std::abs(diag.sweepAngle);
    if (!(std::isfinite(extent) && extent > diag.numericFloor))
        diag.circleState = CircleFitState::UnresolvedExtent;
    else if (maxError > tolerance)
        diag.circleState = CircleFitState::ResidualExceeded;
    else
    {
        diag.circleState = CircleFitState::Accepted;
        out.hasCircle = diag.hasCircle = true;
    }
    return out;
}
static std::optional<AnalyticParameters> DirectionalCone(const std::array<BoundaryInfo, 4> &edges,
                                                         const AnalyticParameters &legacy,
                                                         const RecognitionOptions &o)
{
    if (legacy.type != sggk::SurfaceType::Cone)
        return {};
    std::vector<V> normals, generators;
    V sum{}, reference = Get(legacy.axis);
    for (const auto &edge : edges)
    {
        if (edge.kind == BoundaryKind::Circle)
        {
            V n = Unit(Get(edge.circle.normal));
            if (Dot(n, reference) < 0)
                n = n * (-1);
            normals.push_back(n);
            sum = sum + n;
        }
        else if (edge.kind == BoundaryKind::Line)
            generators.push_back(Unit(Get(edge.line.direction)));
    }
    if (normals.size() != 2 || generators.size() != 2 || Norm(sum) < 1e-12)
        return {};
    V axis = Unit(sum);
    for (V n : normals)
        if (Angle(n, axis) > o.angularTolerance)
            return {};
    double halfAngle = 0;
    std::array<double, 2> angles;
    for (size_t i = 0; i < generators.size(); ++i)
    {
        angles[i] =
            std::atan2(Norm(Cross(generators[i], axis)), std::abs(Dot(generators[i], axis)));
        halfAngle += .5 * angles[i];
    }
    if (!(halfAngle > 0 && halfAngle < Pi / 2) ||
        std::abs(angles[0] - angles[1]) > o.angularTolerance)
        return {};
    V apex = Get(legacy.origin);
    for (const auto &edge : edges)
    {
        if (edge.kind == BoundaryKind::Line &&
            Norm(Cross(Get(edge.line.origin) - apex, Get(edge.line.direction))) >
                o.distanceTolerance + Roundoff(Norm(apex)))
            return {};
        if (edge.kind == BoundaryKind::Circle && Norm(Cross(Get(edge.circle.center) - apex, axis)) >
                                                     o.distanceTolerance + Roundoff(Norm(apex)))
            return {};
    }
    auto candidate = legacy;
    candidate.axis = Direction(axis);
    candidate.halfAngle = halfAngle;
    return candidate;
}
std::vector<AnalyticPatch> RecoverRankedCandidates(const ClassificationSamples &cache,
                                                   const std::vector<ApproximateSeed> &seeds,
                                                   double distanceTolerance,
                                                   ClassificationDiagnostics *diagnostics,
                                                   std::vector<RecoveredCandidate> *recovered)
{
    if (recovered)
        recovered->clear();
    const size_t missing = seeds.size();
    struct Selection
    {
        size_t best, boundary;
    };
    struct Family
    {
        sggk::SurfaceType type;
        Selection legacy, directional;
        double score;
    };
    std::vector<std::optional<AnalyticPatch>> patches(seeds.size());
    std::vector<PatchRecoveryTrace> traces(recovered ? seeds.size() : 0);
    std::vector<bool> attempted(seeds.size(), false);
    auto recover = [&](size_t i) {
        if (!attempted[i])
        {
            attempted[i] = true;
            if (diagnostics)
                ++diagnostics->patchRecoveryAttempts;
            try
            {
                patches[i] = RecoverApproxPatch(cache, seeds[i].parameters,
                                                recovered ? &traces[i] : nullptr);
            }
            catch (const std::runtime_error &)
            {
            }
        }
        return patches[i].has_value();
    };
    auto select = [&](sggk::SurfaceType type, bool directional) {
        std::vector<size_t> indices;
        for (size_t i = 0; i < seeds.size(); ++i)
            if (seeds[i].parameters.type == type && seeds[i].directionalCone == directional &&
                std::isfinite(seeds[i].score))
                indices.push_back(i);
        // Stable ties preserve the original enumeration order.
        std::stable_sort(indices.begin(), indices.end(),
                         [&](size_t a, size_t b) { return seeds[a].score < seeds[b].score; });
        Selection chosen{missing, missing};
        for (size_t i : indices)
            if (recover(i))
            {
                chosen.best = i;
                break;
            }
        if (chosen.best != missing && seeds[chosen.best].fromBoundary)
            chosen.boundary = chosen.best;
        else
            for (size_t i : indices)
                if (seeds[i].fromBoundary && recover(i))
                {
                    chosen.boundary = i;
                    break;
                }
        return chosen;
    };
    std::vector<Family> families;
    for (const auto &seed : seeds)
    {
        auto type = seed.parameters.type;
        if (std::any_of(families.begin(), families.end(),
                        [&](const Family &f) { return f.type == type; }))
            continue;
        auto legacy = select(type, false), directional = select(type, true);
        double score = INFINITY;
        for (auto selection : {legacy, directional})
            if (selection.best != missing)
                score = std::min(score, seeds[selection.best].score);
        families.push_back({type, legacy, directional, score});
    }
    const bool hasSphere = std::any_of(families.begin(), families.end(), [&](const Family &f) {
        return f.type == sggk::SurfaceType::Sphere && f.score < distanceTolerance;
    });
    // This margin changes candidate priority only; acceptance still uses the
    // requested distance tolerance. Preserve ordinary ring-torus priority.
    constexpr double nearSphereRadiusFactor = 5.0;
    const double nearSphereRadiusTolerance = nearSphereRadiusFactor * distanceTolerance;
    auto deferTorus = [&](const Family &f) {
        if (!hasSphere || f.type != sggk::SurfaceType::Torus || f.legacy.best == missing)
            return false;
        const auto &parameters = seeds[f.legacy.best].parameters;
        const double majorRadius = std::abs(parameters.majorRadius);
        return majorRadius <= nearSphereRadiusTolerance && majorRadius < parameters.radius;
    };
    std::stable_sort(families.begin(), families.end(), [&](const Family &a, const Family &b) {
        if (deferTorus(a) != deferTorus(b))
            return !deferTorus(a);
        return a.score < b.score;
    });
    std::vector<AnalyticPatch> result;
    auto append = [&](size_t i) {
        result.push_back(*patches[i]);
        if (recovered)
            recovered->push_back({*patches[i], traces[i], seeds[i].parameters, true, cache.stamp});
    };
    for (const auto &family : families)
        for (auto selection : {family.directional, family.legacy})
        {
            // New directional cone seeds are tried first; both original fallback
            // roles survive. A failed chart never discards the next ranked seed.
            if (selection.boundary != missing && selection.boundary != selection.best)
                append(selection.boundary);
            if (selection.best != missing)
                append(selection.best);
        }
    return result;
}
std::vector<AnalyticPatch> ClassifyApproximate(const sggk::BSplineSurface &s,
                                               const ConversionOptions &o,
                                               ClassificationDiagnostics *diagnostics,
                                               bool directionalCones, SurfaceSampleContext *shared,
                                               std::vector<RecoveredCandidate> *recovered,
                                               ClassificationSamples *savedSamples)
{
    std::array<std::vector<CurveFits>, 2> groups;
    unsigned n = o.sampling.interiorIsoCurves;
    if (n > 5)
        throw Failure(Status::InvalidInput, "At most five internal isocurves per direction");
    double tolerance = o.validation.distanceTolerance * o.refinement.candidateToleranceFactor;
    for (auto &group:groups) group.resize(n+2);
    auto fit=[&](unsigned g,unsigned j,bool boundary) {
        {
            double a = g == 0 ? s.MinParamU() : s.MinParamV(),
                   b = g == 0 ? s.MaxParamU() : s.MaxParamV();
            double t = a + (b - a) * double(j) / (n + 1);
            auto curve = g == 0 ? s.CalcVCurve(t) : s.CalcUCurve(t);
            std::function<V(double)> point;
            if (shared) point=[&,g,t](double varying) {
                return shared->At(g==0?t:varying,g==0?varying:t,
                    sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
            };
            IsoCurveFitDiagnostic local;
            if (diagnostics) diagnostics->isoCurveFits.emplace_back();
            auto &fitDiagnostic = diagnostics ? diagnostics->isoCurveFits.back() : local;
            fitDiagnostic.fixedDirection = g;
            fitDiagnostic.fixedParameter = t;
            fitDiagnostic.boundary = boundary;
            try { groups[g][j]=FitApproximateIsoCurve(curve,
                                         static_cast<BoundaryId>(2 * g + (j == n + 1 ? 1 : 0)),
                                         tolerance, o.sampling.maxSamples,&fitDiagnostic,point,
                                         IsoCurveFitOptions{8, 0, true}); }
            catch (const Failure &e)
            {
                if (boundary || e.status==Status::InvalidInput || e.status==Status::VerificationFailed) throw;
                return;
            }
            // Interior families must never bypass an unsupported outer edge.
            if (boundary && !groups[g][j].hasLine && !groups[g][j].hasCircle)
            {
                bool numerical = fitDiagnostic.circleState != CircleFitState::ResidualExceeded;
                std::ostringstream message;
                message << "Boundary fit failed: edge=" << unsigned(fitDiagnostic.edge)
                        << ", circle=" << CircleFitStateName(fitDiagnostic.circleState)
                        << ", line_error=" << fitDiagnostic.lineError
                        << ", circle_error=" << fitDiagnostic.circleError;
                throw Failure(numerical ? Status::IllConditioned : Status::InconsistentBoundary,
                              message.str());
            }
        }
    };
    // Gate all four outside edges before sampling any internal families.
    for (unsigned g=0;g<2;++g) {fit(g,0,true);fit(g,n+1,true);}
    for (unsigned g=0;g<2;++g) for (unsigned j=1;j<=n;++j) fit(g,j,false);
    std::vector<ApproximateSeed> seeds;
    std::optional<ClassificationSamples> cache;
    auto add = [&](const AnalyticParameters &parameters, bool boundary, bool directional) {
        if (!cache)
            cache = BuildClassificationSamples(s, diagnostics, shared);
        double score = 0;
        for (V point : cache->scorePoints)
            score = std::max(score, Distance(parameters, Point(point)));
        seeds.push_back({parameters, score, boundary, directional});
        if (diagnostics)
        {
            ++diagnostics->rawCandidateCount;
            diagnostics->directionalConeCandidates += directional ? 1 : 0;
        }
    };
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
                    const bool boundary = pair0 == pairs.front() && pair1 == pairs.front();
                    for (const auto &c : Reconstruct(edges, ro))
                    {
                        add(c.parameters, boundary, false);
                        if (directionalCones)
                        {
                            try
                            {
                                if (auto cone = DirectionalCone(edges, c.parameters, ro))
                                    add(*cone, boundary, true);
                            }
                            catch (const std::runtime_error &)
                            {
                                // The legacy seed remains available when the new frame fails.
                            }
                        }
                    }
                }
    if (recovered)
        recovered->clear();
    if (!cache)
        return {};
    auto result = RecoverRankedCandidates(*cache, seeds, o.validation.distanceTolerance,
                                          diagnostics, recovered);
    if (savedSamples)
        *savedSamples = std::move(*cache);
    return result;
}

struct IsoClassificationSearch::Impl
{
    SurfaceSampleContext &samples;
    ConversionOptions options;
    ClassificationDiagnostics &diag;
    struct Curve
    {
        sggk::Curve3DPtr curve;
        CurveFits fits;
        unsigned density = 0;
        bool refined = false;
    };
    std::array<std::vector<Curve>, 2> curves;
    struct Score
    {
        AnalyticParameters parameters;
        size_t count = 0;
        double maximum = 0;
    };
    struct Seed
    {
        AnalyticParameters parameters;
        size_t score;
        bool consumed = false;
    };
    std::vector<Score> scores;
    std::vector<Seed> seeds;
    std::vector<std::pair<AnalyticParameters, unsigned>> tried;
    std::vector<std::array<double, 2>> scoreUV;
    std::vector<V> scorePoints;
    ClassificationSamples cache;
    unsigned cacheDensity = 0, stage = 0, density = 2, n;
    bool stageReady = false;
    double tolerance;

    Impl(SurfaceSampleContext &context, const ConversionOptions &o, ClassificationDiagnostics &d)
        : samples(context), options(o), diag(d), n(o.sampling.interiorIsoCurves),
          tolerance(o.validation.distanceTolerance*o.refinement.candidateToleranceFactor)
    {
        if (n > 5 || !o.sampling.maxSamples ||
            !(o.validation.distanceTolerance > 0 && std::isfinite(o.validation.distanceTolerance)) ||
            !(tolerance > 0 && std::isfinite(tolerance)))
            throw Failure(Status::InvalidInput, "Invalid incremental classification options");
        for (auto &group : curves) group.resize(n+2);
    }
    double Fixed(unsigned g, unsigned j) const
    {
        const auto &s = samples.source;
        double a = g ? s.MinParamV() : s.MinParamU();
        double b = g ? s.MaxParamV() : s.MaxParamU();
        return a+(b-a)*double(j)/(n+1);
    }
    V CurvePoint(unsigned g, unsigned j, double t)
    {
        bool inserted = false;
        double fixed = Fixed(g, j);
        V point = samples.AtCurvePoint(g ? t : fixed, g ? fixed : t,
            [&] { return Get(curves[g][j].curve->CalcPoint(t)); }, inserted).point;
        if (inserted) ++diag.curvePointEvaluations;
        else ++diag.pointCacheHits;
        return point;
    }
    void Fit(unsigned g, unsigned j, bool force)
    {
        auto &c = curves[g][j];
        if (c.density == density && (c.refined || !force)) return;
        if (!c.curve)
        {
            IsoMeasure time(diag.extractionSamplingMs);
            c.curve = g ? samples.source.CalcUCurve(Fixed(g,j)) : samples.source.CalcVCurve(Fixed(g,j));
            ++diag.curveExtractions;
        }
        if (c.density && density > c.density) ++diag.samplingRefinements;
        IsoCurveFitDiagnostic d;
        d.fixedDirection = g; d.fixedParameter = Fixed(g,j);
        d.boundary = j == 0 || j == n+1;
        IsoCurveFitOptions fitOptions;
        fitOptions.subdivisionsPerSpan = density;
        fitOptions.fastTolerance = options.validation.distanceTolerance*.25;
        fitOptions.forceGeometric = force;
        c.fits = {};
        try
        {
            c.fits = FitApproximateIsoCurve(c.curve,
                static_cast<BoundaryId>(2*g+(j==n+1)), tolerance, options.sampling.maxSamples,
                &d, [&](double t) { return CurvePoint(g,j,t); }, fitOptions);
        }
        catch (const Failure &e)
        {
            if (e.status == Status::InvalidInput || e.status == Status::VerificationFailed) throw;
            // Sparse rank/extent failures request further samples, not acceptance
            // or a permanent rejection of a possibly short circular boundary.
        }
        c.density = density; c.refined = force;
        ++diag.fitPasses;
        diag.geometricIterations += d.geometricIterations;
        diag.fastCircleFits += d.skippedGeometricFit ? 1 : 0;
        diag.extractionSamplingMs += d.samplingMs;
        diag.algebraicFitMs += d.algebraicMs;
        diag.geometricFitMs += d.geometricMs;
        diag.isoCurveFits.push_back(d);
    }
    static bool SameScoreMother(const AnalyticParameters &a, const AnalyticParameters &b)
    {
        auto equal = [](V x, V y) { return x.x==y.x && x.y==y.y && x.z==y.z; };
        if (a.type != b.type || !equal(Get(a.origin),Get(b.origin)) ||
            a.radius != b.radius || a.majorRadius != b.majorRadius || a.halfAngle != b.halfAngle)
            return false;
        if (a.type == sggk::SurfaceType::Sphere) return true;
        // No approximate merges: distinct small errors may decide certification.
        return equal(Get(a.axis),Get(b.axis));
    }
    void Add(const AnalyticParameters &p)
    {
        ++diag.rawCandidateCount;
        if (std::any_of(tried.begin(), tried.end(), [&](const auto &x) {
                return x.second == density && SameMotherParameters(x.first,p); }) ||
            std::any_of(seeds.begin(), seeds.end(), [&](const Seed &x) {
                return SameMotherParameters(x.parameters,p); }))
        { ++diag.duplicateCandidates; return; }
        size_t i = 0;
        for (; i < scores.size(); ++i) if (SameScoreMother(scores[i].parameters,p)) break;
        if (i == scores.size()) scores.push_back({p,0,0});
        else ++diag.sharedMotherScores;
        seeds.push_back({p,i,false});
    }
    void BuildScores()
    {
        if (!scoreUV.empty()) return;
        IsoMeasure time(diag.extractionSamplingMs);
        auto uv = samples.ScoreParameters();
        if (uv.empty()) throw Failure(Status::Degenerate,"No independent classification score points");
        std::vector<double> nearest(uv.size(), INFINITY);
        std::vector<bool> selected(uv.size(), false);
        double du = samples.source.MaxParamU()-samples.source.MinParamU();
        double dv = samples.source.MaxParamV()-samples.source.MinParamV();
        size_t next = 0;
        for (size_t k = 0; k < std::min(size_t(9),uv.size()); ++k)
        {
            selected[next] = true; scoreUV.push_back(uv[next]);
            for (size_t i = 0; i < uv.size(); ++i)
            {
                double u = (uv[i][0]-uv[next][0])/du, v = (uv[i][1]-uv[next][1])/dv;
                nearest[i] = selected[i] ? -1 : std::min(nearest[i],u*u+v*v);
            }
            next = size_t(std::max_element(nearest.begin(),nearest.end())-nearest.begin());
        }
        for (size_t i=0;i<uv.size();++i) if (!selected[i]) scoreUV.push_back(uv[i]);
    }
    void CompleteScore(size_t index, size_t count)
    {
        auto &score = scores[index];
        if (score.count >= count) return;
        {
            IsoMeasure time(diag.extractionSamplingMs);
            while (scorePoints.size() < count)
            {
                auto uv = scoreUV[scorePoints.size()];
                unsigned before = samples.Evaluations();
                scorePoints.push_back(samples.At(uv[0],uv[1],sggk::DerivType::NONE,
                    sggk::DerivType::NONE,false).point);
                if (samples.Evaluations() > before) ++diag.sourcePointEvaluations;
                else ++diag.pointCacheHits;
            }
        }
        IsoMeasure time(diag.candidateScoringMs);
        while (score.count < count)
        {
            double d = Distance(score.parameters,Point(scorePoints[score.count++]));
            ++diag.distanceEvaluations;
            score.maximum = std::isfinite(d) ? std::max(score.maximum,d) : INFINITY;
        }
    }
    bool BuildStage()
    {
        density = stage < 4 ? 2 : stage == 4 ? 4 : 8;
        ++diag.classificationStages;
        for (unsigned g=0;g<2;++g) { Fit(g,0,stage>0); Fit(g,n+1,stage>0); }
        for (const auto &group : curves)
            for (unsigned j : {0u,n+1})
                if (!group[j].fits.hasLine && !group[j].fits.hasCircle) return false;
        std::vector<unsigned> active{0,n+1};
        if (stage >= 2 && n) active.push_back((n+1)/2);
        if (stage >= 3) for (unsigned j=1;j<=n;++j)
            if (std::find(active.begin(),active.end(),j)==active.end()) active.push_back(j);
        for (unsigned g=0;g<2;++g)
            for (unsigned j:active) Fit(g,j,stage>0);
        std::vector<std::pair<unsigned,unsigned>> pairs{{0,n+1}};
        if (stage >= 2 && n)
        {
            pairs.push_back({0,(n+1)/2}); pairs.push_back({(n+1)/2,n+1});
        }
        if (stage >= 3 && n>=2) pairs.push_back({1,n});
        seeds.clear();
        {
            IsoMeasure time(diag.candidateScoringMs);
            RecognitionOptions ro; ro.distanceTolerance=tolerance; ro.angularTolerance=1e-3;
            for (unsigned k0=0;k0<2;++k0) for (unsigned k1=0;k1<2;++k1)
                for (auto pair0:pairs) for (auto pair1:pairs)
                {
                    std::array<BoundaryInfo,4> edges;
                    bool valid = true;
                    for (unsigned g=0;g<2;++g)
                    {
                        auto pair=g?pair1:pair0; unsigned kind=g?k1:k0;
                        for (unsigned j=0;j<2;++j)
                        {
                            const auto &f=curves[g][j?pair.second:pair.first].fits;
                            valid &= kind?f.hasCircle:f.hasLine;
                            edges[2*g+j]=kind?f.circle:f.line;
                        }
                    }
                    if (!valid) continue;
                    for (const auto &candidate:Reconstruct(edges,ro,false))
                    {
                        Add(candidate.parameters);
                        try
                        {
                            if (auto cone=DirectionalCone(edges,candidate.parameters,ro))
                            { Add(*cone); ++diag.directionalConeCandidates; }
                        }
                        catch (const std::runtime_error &) { /* Preserve the original cone seed. */ }
                    }
                }
        }
        if (seeds.empty()) return false;
        BuildScores();
        for (const auto &seed:seeds) CompleteScore(seed.score,std::min(size_t(9),scoreUV.size()));
        return true;
    }
    void BuildBoundaryCache()
    {
        if (cacheDensity == density) return;
        IsoMeasure time(diag.extractionSamplingMs);
        ClassificationSamples next; // New immutable sampling identity, even for the same source.
        const auto &s=samples.source;
        next.sourceDomain={s.MinParamU(),s.MaxParamU(),s.MinParamV(),s.MaxParamV()};
        for (unsigned g=0;g<2;++g) for (unsigned last=0;last<2;++last)
        {
            unsigned e=2*g+last, j=last?n+1:0;
            auto bc=std::dynamic_pointer_cast<sggk::BSplineCurve3D>(curves[g][j].curve);
            next.boundaryParameters[e]=SampleParameters(*bc,density);
            for (double t:next.boundaryParameters[e]) next.boundaries[e].push_back(CurvePoint(g,j,t));
        }
        double u=(s.MinParamU()+s.MaxParamU())*.5, v=(s.MinParamV()+s.MaxParamV())*.5;
        const auto &center=samples.At(u,v);
        next.centerPoint=center.point; next.centerDu=center.du; next.centerDv=center.dv;
        ++diag.cachedDerivativeEvaluations; ++diag.sampleCacheBuilds;
        cache=std::move(next); cacheDensity=density;
    }
    std::optional<RecoveredCandidate> Next(ClassificationSamples &saved)
    {
        while (stage < 6)
        {
            samples.CheckDeadline("iso_search");
            if (!stageReady)
            {
                if (!BuildStage()) { ++stage; continue; }
                stageReady=true;
            }
            // Completing the current smallest lower bound avoids a full score
            // of candidates that cannot win. Equal maxima keep stable order.
            size_t best=seeds.size();
            bool sphere=false;
            for (const auto &s:seeds)
                if (!s.consumed && s.parameters.type==sggk::SurfaceType::Sphere &&
                    scores[s.score].maximum < options.validation.distanceTolerance) sphere=true;
            auto deferred=[&](const Seed &s) {
                return sphere && s.parameters.type==sggk::SurfaceType::Torus &&
                    std::abs(s.parameters.majorRadius)<=5*options.validation.distanceTolerance &&
                    std::abs(s.parameters.majorRadius)<s.parameters.radius;
            };
            auto better=[&](size_t a,size_t b) {
                if (deferred(seeds[a]) != deferred(seeds[b])) return !deferred(seeds[a]);
                double x=scores[seeds[a].score].maximum,y=scores[seeds[b].score].maximum;
                return x<y || (x==y && a<b);
            };
            for (;;)
            {
                best=seeds.size();
                for (size_t i=0;i<seeds.size();++i)
                    if (!seeds[i].consumed && (best==seeds.size() || better(i,best))) best=i;
                if (best==seeds.size()) break;
                if (scores[seeds[best].score].count == scoreUV.size()) break;
                CompleteScore(seeds[best].score,scoreUV.size());
                sphere=false;
                for (const auto &s:seeds)
                    if (!s.consumed && s.parameters.type==sggk::SurfaceType::Sphere &&
                        scores[s.score].maximum<options.validation.distanceTolerance) sphere=true;
            }
            if (best==seeds.size()) { ++stage; stageReady=false; continue; }
            auto &seed=seeds[best]; seed.consumed=true;
            tried.push_back({seed.parameters,density});
            // Scores select seeds only. Every returned candidate must still be
            // continuously certified; a poor seed can be replaced by the next stage.
            if (scores[seed.score].maximum > tolerance) continue;
            BuildBoundaryCache();
            IsoMeasure time(diag.patchRecoveryMs);
            ++diag.patchRecoveryAttempts;
            try
            {
                RecoveredCandidate result;
                result.patch=RecoverApproxPatch(cache,seed.parameters,&result.trace);
                result.traceParameters=seed.parameters; result.hasTrace=true;
                result.sourceStamp=cache.stamp; saved=cache;
                return result;
            }
            catch (const Failure &e)
            {
                if (e.status==Status::InvalidInput) throw;
            }
            catch (const std::runtime_error &) { /* Try the next chart. */ }
        }
        return {};
    }
};
IsoClassificationSearch::IsoClassificationSearch(SurfaceSampleContext &samples,
    const ConversionOptions &options, ClassificationDiagnostics &diagnostics)
    : impl(std::make_unique<Impl>(samples,options,diagnostics)) {}
IsoClassificationSearch::~IsoClassificationSearch() = default;
std::optional<RecoveredCandidate> IsoClassificationSearch::Next(ClassificationSamples &samples)
{
    return impl->Next(samples);
}
} // namespace bspline2analytic::detail
