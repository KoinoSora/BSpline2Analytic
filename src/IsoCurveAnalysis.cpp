#include "ApproximateInternal.h"
#include <numeric>
#include <sstream>
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
                                const std::function<V(double)> &sample)
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
    diag.edge = id;
    diag.candidateTolerance = tolerance;
    if (!(tolerance > 0 && std::isfinite(tolerance)) || !maxSamples)
        throw Failure(Status::InvalidInput, "Invalid isocurve fitting options");
    auto bc = std::dynamic_pointer_cast<sggk::BSplineCurve3D>(curve);
    if (!bc)
        throw Failure(Status::Unsupported, "Expected spline isocurve");
    auto ts = SampleParameters(*bc, 8);
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
        // Midpoint chord quadrature is representation-independent to sampling accuracy.
        double parameter=(ts[i]+ts[i-1])*.5;
        V mid=sample ? sample(parameter) : Get(curve->CalcPoint(parameter));
        if (!Finite(mid))
        {
            diag.circleState = CircleFitState::NonFinite;
            throw Failure(Status::IllConditioned, "Nonfinite evaluated isocurve midpoint");
        }
        length = std::max(length, Norm(mid - points[i - 1]) + Norm(points[i] - mid));
        weights[i - 1] += length * 0.5;
        weights[i] += length * 0.5;
    }
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
        if (!(rr > 0 && std::isfinite(rr)) || !Finite(trial) || !Finite(nn))
            break;
        auto next = residual(trial, nn, rr);
        double old = 0, now = 0;
        for (double x : f)
            old += x * x;
        for (double x : next)
            now += x * x;
        if (!std::isfinite(now) || now >= old)
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
                                         tolerance, o.sampling.maxSamples,&fitDiagnostic,point); }
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
} // namespace bspline2analytic::detail
