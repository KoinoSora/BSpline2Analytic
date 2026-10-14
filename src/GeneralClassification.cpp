#include "GeneralInternal.h"
#include <random>
namespace bspline2analytic::detail
{
using Type = sggk::SurfaceType;
static unsigned Dof(Type t)
{
    switch (t)
    {
    case Type::Plane: return 3;
    case Type::Sphere: return 4;
    case Type::Cylinder: return 5;
    case Type::Cone: return 6;
    case Type::Torus: return 7;
    default: return 0;
    }
}
static AnalyticParameters Frame(Type type, V center, V axis)
{
    AnalyticParameters a;
    a.type = type; a.origin = Point(center); a.axis = Direction(axis);
    a.xDirection = Direction(Perpendicular(Unit(axis)));
    return a;
}
static AnalyticParameters PlaneSeed(const std::array<SurfaceSample, 4> &p, double scale)
{
    V a = (p[1].point - p[0].point) / scale, b = (p[2].point - p[0].point) / scale;
    if (Norm(Cross(a, b)) < 1e-8 * Norm(a) * Norm(b))
        throw Failure(Status::IllConditioned, "Collinear plane seed");
    return Frame(Type::Plane, p[0].point, Cross(a, b));
}
static AnalyticParameters SphereSeed(const std::array<SurfaceSample, 4> &p, double scale)
{
    V n = p[0].normal, m = p[1].normal;
    double b = Dot(n, m), divisor = 1 - b * b;
    if (divisor < 1e-8) throw Failure(Status::IllConditioned, "Parallel sphere normals");
    V delta = (p[1].point - p[0].point) / scale;
    double t = (Dot(n, delta) - b * Dot(m, delta)) / divisor;
    double u = (b * Dot(n, delta) - Dot(m, delta)) / divisor;
    V c = (p[0].point + n * (t * scale) + p[1].point + m * (u * scale)) * .5;
    auto a = Frame(Type::Sphere, c, {0, 0, 1});
    a.radius = (Norm(p[0].point - c) + Norm(p[1].point - c)) * .5;
    return a;
}
static AnalyticParameters CylinderSeed(const std::array<SurfaceSample, 4> &p, double scale)
{
    V cross = Cross(p[0].normal, p[1].normal);
    if (Norm(cross) < 1e-4) throw Failure(Status::IllConditioned, "Parallel cylinder normals");
    V axis = Unit(cross), delta = (p[1].point - p[0].point) / scale;
    delta = delta - axis * Dot(delta, axis);
    double b = Dot(p[0].normal, p[1].normal), divisor = 1 - b * b;
    double t = (Dot(p[0].normal, delta) - b * Dot(p[1].normal, delta)) / divisor;
    V c = p[0].point + p[0].normal * (t * scale);
    auto a = Frame(Type::Cylinder, c, axis);
    a.radius = std::abs(t * scale);
    return a;
}
static AnalyticParameters ConeSeed(const std::array<SurfaceSample, 4> &p, double scale)
{
    std::vector<std::array<double, 3>> rows;
    std::vector<double> rhs;
    V origin = p[0].point;
    for (unsigned i = 0; i < 3; ++i)
    {
        V n = p[i].normal;
        rows.push_back({n.x, n.y, n.z});
        rhs.push_back(Dot(n, (p[i].point - origin) / scale));
    }
    std::array<double, 3> x;
    if (!LeastSquares3(rows, rhs, 1e-10, x))
        throw Failure(Status::IllConditioned, "Singular cone tangent planes");
    V apex = origin + V{x[0], x[1], x[2]} * scale;
    V d0 = Unit(p[0].point - apex), d1 = Unit(p[1].point - apex), d2 = Unit(p[2].point - apex);
    V cross = Cross(d1 - d0, d2 - d0);
    if (Norm(cross) < 1e-8) throw Failure(Status::IllConditioned, "Singular cone directions");
    V axis = Unit(cross);
    if (Dot(axis, d0 + d1 + d2) < 0) axis = axis * -1;
    auto a = Frame(Type::Cone, apex, axis);
    a.halfAngle = (std::acos(std::clamp(Dot(axis, d0), -1.0, 1.0)) +
                   std::acos(std::clamp(Dot(axis, d1), -1.0, 1.0)) +
                   std::acos(std::clamp(Dot(axis, d2), -1.0, 1.0))) / 3;
    if (!(a.halfAngle > 1e-8 && a.halfAngle < Pi / 2 - 1e-8))
        throw Failure(Status::IllConditioned, "Invalid cone seed angle");
    return a;
}
// Solve the four-normal-line transversal equations in normalized coordinates.
// All axes from the starts are retained for scoring; this is a bounded numerical
// alternative to the paper's algebraic transversal construction, not an exact
// enumeration guarantee. For each axis, fit a meridian circle to the four points.
static std::vector<AnalyticParameters> TorusSeeds(const std::array<SurfaceSample, 4> &p, double scale)
{
    V mean{};
    for (const auto &s : p) mean = mean + s.point * .25;
    std::array<V, 4> points;
    for (unsigned i = 0; i < 4; ++i) points[i] = (p[i].point - mean) / scale;
    std::vector<V> starts{{1,0,0}, {0,1,0}, {0,0,1}, {1,1,1}, {1,-1,1}, {-1,1,1}};
    for (unsigned i = 1; i < 4; ++i)
    {
        V x = Cross(p[0].normal, p[i].normal);
        if (Norm(x) > 1e-4) starts.push_back(x);
    }
    std::vector<AnalyticParameters> out;
    for (V start : starts)
    {
        V axis = Unit(start), center{};
        auto residual = [&](V c, V a) {
            std::vector<double> r;
            for (unsigned i = 0; i < 4; ++i)
                r.push_back(Dot(points[i] - c, Cross(p[i].normal, a)));
            return r;
        };
        auto cost = [](const std::vector<double> &r) { double v = 0; for (double x : r) v += x*x; return v; };
        double lambda = 1e-6;
        for (unsigned iteration = 0; iteration < 40; ++iteration)
        {
            auto r = residual(center, axis);
            if (cost(r) < 1e-20) break;
            V ex = Perpendicular(axis), ey = Cross(axis, ex);
            auto move = [&](const std::vector<double> &d) {
                V a = Unit(axis + ex * d[2] + ey * d[3]);
                V c = center + ex * d[0] + ey * d[1];
                return std::make_pair(c - a * Dot(c, a), a);
            };
            std::vector<std::vector<double>> jac(4, std::vector<double>(4));
            for (unsigned j = 0; j < 4; ++j)
            {
                std::vector<double> d(4, 0); d[j] = 1e-6;
                auto plus = move(d); d[j] = -1e-6; auto minus = move(d);
                auto rp = residual(plus.first, plus.second), rm = residual(minus.first, minus.second);
                for (unsigned i = 0; i < 4; ++i) jac[i][j] = (rp[i] - rm[i]) / 2e-6;
            }
            std::vector<double> rhs;
            for (double x : r) rhs.push_back(-x);
            for (unsigned j = 0; j < 4; ++j)
            {
                std::vector<double> row(4, 0); row[j] = std::sqrt(lambda);
                jac.push_back(row); rhs.push_back(0);
            }
            std::vector<double> d;
            if (!SolveLeastSquares(jac, rhs, d)) break;
            double length = 0; for (double x : d) length = std::max(length, std::abs(x));
            if (length > .5) for (double &x : d) x *= .5 / length;
            bool accepted = false;
            for (unsigned line = 0; line < 8; ++line)
            {
                auto trial = move(d);
                if (cost(residual(trial.first, trial.second)) < cost(r))
                { center = trial.first; axis = trial.second; accepted = true; break; }
                for (double &x : d) x *= .5;
            }
            lambda = accepted ? std::max(1e-12, lambda * .3) : lambda * 10;
            if (lambda > 1e8) break;
        }
        if (cost(residual(center, axis)) > 1e-8) continue;
        std::vector<std::array<double, 3>> rows;
        std::vector<double> rhs;
        for (V point : points)
        {
            V d = point - center;
            double z = Dot(d, axis), rho = Norm(d - axis * z);
            rows.push_back({2 * rho, 2 * z, 1}); rhs.push_back(rho*rho + z*z);
        }
        std::array<double, 3> circle;
        if (!LeastSquares3(rows, rhs, 1e-10, circle)) continue;
        double r2 = circle[2] + circle[0]*circle[0] + circle[1]*circle[1];
        if (!(r2 > 0 && circle[0] > std::sqrt(r2))) continue;
        auto a = Frame(Type::Torus, mean + (center + axis * circle[1]) * scale, axis);
        a.radius = std::sqrt(r2) * scale; a.majorRadius = circle[0] * scale;
        out.push_back(a);
    }
    return out;
}
static void ScoreInterior(GeneralScoringCandidate &candidate, SurfaceSampleContext &context,
                          const std::vector<size_t> &indices, GeneralClassificationDiagnostics &diagnostics)
{
    const auto &a = candidate.seed.parameters;
    auto &score = candidate.seed.score;
    const auto &points = context.Greville();
    for (size_t i : indices)
    {
        ++diagnostics.interiorPointEvaluations;
        const auto &s = points[i];
        V q = ProjectToMother(a, s.point);
        double d = Norm(s.point - q);
        if (!std::isfinite(d)) throw Failure(Status::IllConditioned, "Nonfinite interior distance");
        score.interiorMax = std::max(score.interiorMax, d);
        double weight = s.area > 0 ? s.area : 0;
        candidate.weightedSquaredSum += weight * d*d; candidate.weightSum += weight;
        if (s.normalValid)
        {
            auto f = FrameOf(a);
            V local = f.local(q), n;
            switch (a.type)
            {
            case Type::Plane: n = f.z; break;
            case Type::Sphere: n = Unit(q - f.o); break;
            case Type::Cylinder: n = Unit(f.vector({local.x, local.y, 0})); break;
            case Type::Cone: {
                V radial = f.vector({local.x, local.y, 0});
                double length = Norm(radial);
                // Projection onto the apex has no unique mother normal. Keep
                // its valid distance and report the worst unsigned angle;
                // normal diagnostics must not reject a geometrically valid seed.
                if (!(length > 0)) { score.normalMax = Pi / 2; continue; }
                n = radial / length * std::cos(a.halfAngle) - f.z * std::sin(a.halfAngle); break;
            }
            default: {
                V radial = Unit(f.vector({local.x, local.y, 0}));
                n = Unit(q - f.o - radial * a.majorRadius); break;
            }
            }
            score.normalMax = std::max(score.normalMax,
                std::acos(std::clamp(std::abs(Dot(s.normal, n)), 0.0, 1.0)));
        }
    }
    if (!std::isfinite(candidate.weightedSquaredSum) || !std::isfinite(candidate.weightSum))
        throw Failure(Status::IllConditioned, "Nonfinite score accumulation");
    // A coarse subset may have zero area; the remaining points can still make
    // the full score regular. Such a subset must not reject the candidate.
    score.interiorRms = candidate.weightSum > 0 ?
        std::sqrt(candidate.weightedSquaredSum / candidate.weightSum) : std::numeric_limits<double>::infinity();
}
static void CompleteScore(GeneralScoringCandidate &candidate, SurfaceSampleContext &context,
                          const std::vector<size_t> &remaining, GeneralClassificationDiagnostics &diagnostics)
{
    ++diagnostics.fullScoreCount;
    try
    {
        ScoreInterior(candidate, context, remaining, diagnostics);
        auto &score = candidate.seed.score;
        if (!(candidate.weightSum > 0)) throw Failure(Status::IllConditioned, "No regular independent score points");
        score.boundaryMax = 0;
        for (const auto &edge : context.Boundaries())
            for (const auto &s : edge)
            {
                ++diagnostics.boundaryPointEvaluations;
                double d = Norm(s.point - ProjectToMother(candidate.seed.parameters, s.point));
                if (!std::isfinite(d)) throw Failure(Status::IllConditioned, "Nonfinite boundary distance");
                score.boundaryMax = std::max(score.boundaryMax, d);
            }
        if (!std::isfinite(score.interiorRms) || !std::isfinite(score.normalMax))
            throw Failure(Status::IllConditioned, "Nonfinite general candidate score");
        candidate.stage = GeneralScoreStage::Complete;
        ++diagnostics.completedScoreCount;
    }
    catch (const std::runtime_error &) { candidate.stage = GeneralScoreStage::Failed; }
}
static bool Similar(double a, double b, const GeneralConversionOptions &o)
{
    return std::abs(a - b) <= o.similarScoreAbsoluteRatio * o.fitting.validation.distanceTolerance +
        o.similarScoreRelativeTolerance * std::max(a, b);
}
GeneralCandidatePool InitializeGeneralCandidates(SurfaceSampleContext &context, const GeneralConversionOptions &o,
                                                 GeneralClassificationDiagnostics &diagnostics)
{
    diagnostics = {};
    const auto &points = context.Greville();
    const auto &indices = context.Seeds();
    if (indices.size() < 4 || context.Scores().empty())
        throw Failure(Status::IllConditioned, "Insufficient regular Greville seed/score points");
    V center{};
    for (size_t i : indices) center = center + points[i].point / double(indices.size());
    double scale = 0;
    for (const auto &p : points) scale = std::max(scale, Norm(p.point - center));
    if (!(scale > Roundoff(Norm(center))))
        throw Failure(Status::Degenerate, "Surface samples collapse to a point");
    std::mt19937 random(0x42534132u);
    std::uniform_int_distribution<size_t> pick(0, indices.size() - 1);
    GeneralCandidatePool pool;
    const auto &coarse = context.CoarseScores(o.coarseScorePointCount);
    std::vector<bool> isCoarse(points.size(), false);
    for (size_t i : coarse) isCoarse[i] = true;
    for (size_t i : context.Scores()) if (!isCoarse[i]) pool.remainingScoreIndices.push_back(i);
    unsigned familyIndex = 0;
    for (Type type : {Type::Plane, Type::Sphere, Type::Cylinder, Type::Cone, Type::Torus})
    {
        auto &family = pool.families[familyIndex++];
        for (unsigned attempt = 0; attempt < o.seedAttemptsPerType; ++attempt)
        {
            ++diagnostics.seedAttemptCount;
            std::array<SurfaceSample, 4> p;
            std::array<size_t, 4> chosen{};
            bool complete = true;
            for (unsigned j = 0; j < 4; ++j)
            {
                bool found = false;
                for (unsigned retry = 0; retry < 32 && !found; ++retry)
                {
                    size_t k = indices[pick(random)];
                    found = true;
                    for (unsigned previous = 0; previous < j; ++previous)
                        if (chosen[previous] == k || Norm(points[k].point - p[previous].point) < 1e-6 * scale)
                            found = false;
                    if (found) { chosen[j] = k; p[j] = points[k]; }
                }
                if (!found) { complete = false; break; }
            }
            if (!complete) continue;
            try
            {
                std::vector<AnalyticParameters> seeds;
                switch (type)
                {
                case Type::Plane: seeds.push_back(PlaneSeed(p, scale)); break;
                case Type::Sphere: seeds.push_back(SphereSeed(p, scale)); break;
                case Type::Cylinder: seeds.push_back(CylinderSeed(p, scale)); break;
                case Type::Cone: seeds.push_back(ConeSeed(p, scale)); break;
                default: seeds = TorusSeeds(p, scale); break;
                }
                for (const auto &a : seeds)
                {
                    ++diagnostics.rawCandidateCount;
                    try
                    {
                        bool duplicate = false;
                        for (const auto &c : family)
                        {
                            const auto &existing = c.seed.parameters;
                            V delta = Get(a.origin) - Get(existing.origin);
                            if (type == Type::Plane || type == Type::Cylinder)
                            {
                                V axis = Get(a.axis);
                                delta = type == Type::Plane ? axis * Dot(delta, axis) : delta - axis * Dot(delta, axis);
                            }
                            duplicate |= Norm(delta) < 1e-6 * scale &&
                                std::abs(a.radius - existing.radius) < 1e-6 * scale &&
                                std::abs(a.majorRadius - existing.majorRadius) < 1e-6 * scale &&
                                std::abs(a.halfAngle - existing.halfAngle) < 1e-6 &&
                                (type == Type::Sphere || Parallel(Get(a.axis), Get(existing.axis), 1e-6));
                        }
                        // Deduplication is geometric, before point scoring. It
                        // does not consume extra sample groups to fill a quota.
                        if (duplicate)
                        {
                            ++diagnostics.duplicateCandidateCount;
                            continue;
                        }
                        GeneralScoringCandidate candidate;
                        candidate.seed.parameters = a;
                        candidate.seed.score.type = type;
                        candidate.seed.score.degreesOfFreedom = Dof(type);
                        // Undefined until CompleteScore evaluates the boundary.
                        candidate.seed.score.boundaryMax = std::numeric_limits<double>::infinity();
                        family.push_back(std::move(candidate));
                    }
                    catch (const std::runtime_error &) {}
                }
            }
            catch (const std::runtime_error &) {}
        }
    }
    // Build the same deduplicated mother pool in both modes before scoring.
    // Otherwise a representative failing only during full scoring can suppress
    // different near-duplicates depending on when that failure is discovered.
    for (auto &family : pool.families)
        for (auto &candidate : family)
        {
            ++diagnostics.scoreEvaluationCount;
            try
            {
                if (o.coarseScorePointCount)
                {
                    ++diagnostics.coarseScoreCount;
                    ScoreInterior(candidate, context, coarse, diagnostics);
                }
                else CompleteScore(candidate, context, pool.remainingScoreIndices, diagnostics);
            }
            catch (const std::runtime_error &) { candidate.stage = GeneralScoreStage::Failed; }
        }
    return pool;
}
static double PrimaryScore(const GeneralCandidateScore &s) { return std::max(s.interiorMax, s.boundaryMax); }
std::vector<GeneralSeed> NextGeneralCandidateBatch(GeneralCandidatePool &pool, SurfaceSampleContext &context,
    const GeneralConversionOptions &o, GeneralClassificationDiagnostics &diagnostics)
{
    std::vector<GeneralSeed> candidates;
    for (auto &family : pool.families)
    {
        std::vector<size_t> pending, best;
        for (size_t i = 0; i < family.size(); ++i)
            if (family[i].stage == GeneralScoreStage::Coarse || family[i].stage == GeneralScoreStage::Complete)
                pending.push_back(i);
        auto lower = [&](size_t i) {
            const auto &c = family[i];
            return c.stage == GeneralScoreStage::Complete ? PrimaryScore(c.seed.score) : c.seed.score.interiorMax;
        };
        std::stable_sort(pending.begin(), pending.end(), [&](size_t a, size_t b) {
            if (lower(a) != lower(b)) return lower(a) < lower(b);
            return family[a].seed.score.interiorRms < family[b].seed.score.interiorRms;
        });
        auto better = [&](size_t a, size_t b) {
            const auto &x = family[a].seed.score, &y = family[b].seed.score;
            if (PrimaryScore(x) != PrimaryScore(y)) return PrimaryScore(x) < PrimaryScore(y);
            if (x.interiorRms != y.interiorRms) return x.interiorRms < y.interiorRms;
            return a < b; // Original generation order on exact ties.
        };
        for (size_t i : pending)
        {
            // Only the maximum distance is a lower bound. A coarse RMS is NOT
            // a lower bound on the full weighted RMS. Equal maxima must finish.
            if (best.size() == o.candidatesPerType && lower(i) > PrimaryScore(family[best.back()].seed.score)) break;
            auto &c = family[i];
            if (c.stage == GeneralScoreStage::Coarse) CompleteScore(c, context, pool.remainingScoreIndices, diagnostics);
            if (c.stage != GeneralScoreStage::Complete) continue;
            best.push_back(i);
            std::sort(best.begin(), best.end(), better);
            if (best.size() > o.candidatesPerType) best.pop_back();
        }
        for (size_t i : best)
        {
            candidates.push_back(family[i].seed);
            family[i].stage = GeneralScoreStage::Consumed;
        }
    }
    if (!pool.batches)
        for (const auto &family : pool.families)
            for (const auto &c : family) diagnostics.deferredCandidateCount += c.stage == GeneralScoreStage::Coarse;
    else diagnostics.backupCandidateCount += static_cast<unsigned>(candidates.size());
    if (!candidates.empty()) { ++pool.batches; ++diagnostics.batchCount; }
    std::stable_sort(candidates.begin(), candidates.end(), [](const GeneralSeed &a, const GeneralSeed &b) {
        // Boundary error participates in priority, but never in final acceptance.
        double x = std::max(a.score.interiorMax, a.score.boundaryMax);
        double y = std::max(b.score.interiorMax, b.score.boundaryMax);
        if (x != y) return x < y;
        if (a.score.interiorRms != b.score.interiorRms) return a.score.interiorRms < b.score.interiorRms;
        return a.score.degreesOfFreedom < b.score.degreesOfFreedom;
    });
    // Apply preference after a strict sort: an approximate comparator would not
    // be transitive and would violate std::sort's ordering requirements.
    for (size_t i = 0; i < candidates.size(); ++i)
        if (candidates[i].parameters.type == Type::Torus)
            for (size_t j = i + 1; j < candidates.size(); ++j)
                if (candidates[j].parameters.type == Type::Sphere &&
                    Similar(candidates[i].score.interiorMax, candidates[j].score.interiorMax, o) &&
                    Similar(candidates[i].score.boundaryMax, candidates[j].score.boundaryMax, o))
                {
                    auto sphere = candidates[j];
                    candidates.erase(candidates.begin() + j);
                    candidates.insert(candidates.begin() + i, sphere);
                    break;
                }
    return candidates;
}
std::vector<GeneralSeed> ClassifyGeneral(SurfaceSampleContext &context, const GeneralConversionOptions &o,
                                       GeneralClassificationDiagnostics *diagnostics)
{
    GeneralClassificationDiagnostics local;
    auto &counts = diagnostics ? *diagnostics : local;
    auto pool = InitializeGeneralCandidates(context, o, counts);
    return NextGeneralCandidateBatch(pool, context, o, counts);
}
}
