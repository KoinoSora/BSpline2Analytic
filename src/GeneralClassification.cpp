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
using LineCoordinates = std::array<double, 6>; // direction, moment = position x direction
// Column-pivoted Householder QR of A^T (6 x 4). The last two columns
// of the orthogonal factor span null(A), without forming normal equations.
static bool AxisPencil(const std::array<V, 4> &points, const std::array<SurfaceSample, 4> &p,
                       LineCoordinates &u, LineCoordinates &v)
{
    double a[6][4]{}, q[6][6]{};
    for (unsigned i = 0; i < 6; ++i) q[i][i] = 1;
    for (unsigned j = 0; j < 4; ++j)
    {
        V moment = Cross(points[j], p[j].normal);
        a[0][j] = moment.x; a[1][j] = moment.y; a[2][j] = moment.z;
        a[3][j] = p[j].normal.x; a[4][j] = p[j].normal.y; a[5][j] = p[j].normal.z;
    }
    double first = 0;
    constexpr double eps = std::numeric_limits<double>::epsilon();
    for (unsigned k = 0; k < 4; ++k)
    {
        unsigned pivot = k;
        double best = 0;
        for (unsigned j = k; j < 4; ++j)
        {
            double length = 0;
            for (unsigned i = k; i < 6; ++i) length = std::hypot(length, a[i][j]);
            if (length > best) { best = length; pivot = j; }
        }
        if (!k) first = best;
        if (!(best > 512 * eps * first) || !std::isfinite(best)) return false;
        for (unsigned i = 0; i < 6; ++i) std::swap(a[i][k], a[i][pivot]);
        double alpha = -std::copysign(best, a[k][k]), length = 0;
        LineCoordinates h{};
        for (unsigned i = k; i < 6; ++i) h[i] = a[i][k];
        h[k] -= alpha;
        for (unsigned i = k; i < 6; ++i) length = std::hypot(length, h[i]);
        for (unsigned i = k; i < 6; ++i) h[i] /= length;
        for (unsigned j = k; j < 4; ++j)
        {
            double dot = 0;
            for (unsigned i = k; i < 6; ++i) dot += h[i] * a[i][j];
            for (unsigned i = k; i < 6; ++i) a[i][j] -= 2 * h[i] * dot;
        }
        for (unsigned i = 0; i < 6; ++i)
        {
            double dot = 0;
            for (unsigned j = k; j < 6; ++j) dot += q[i][j] * h[j];
            for (unsigned j = k; j < 6; ++j) q[i][j] -= 2 * dot * h[j];
        }
    }
    for (unsigned i = 0; i < 6; ++i) { u[i] = q[i][4]; v[i] = q[i][5]; }
    return true;
}
static V LinePart(const LineCoordinates &x, unsigned start) { return {x[start], x[start+1], x[start+2]}; }
struct TorusAxis { V direction, position; };
static std::vector<TorusAxis> TorusAxes(const std::array<V, 4> &points,
    const std::array<SurfaceSample, 4> &p, GeneralTorusSeedDiagnostics &diagnostics)
{
    LineCoordinates u, v;
    if (!AxisPencil(points, p, u, v)) { ++diagnostics.rankDeficientGroups; return {}; }
    V ud = LinePart(u, 0), um = LinePart(u, 3), vd = LinePart(v, 0), vm = LinePart(v, 3);
    // A line must satisfy direction . moment = 0. On the nullspace pencil
    // s*u+t*v this is C*s^2+B*s*t+A*t^2=0, including roots with s=0.
    double a = Dot(vd, vm), b = Dot(ud, vm) + Dot(vd, um), c = Dot(ud, um);
    double magnitude = std::max({std::abs(a), std::abs(b), std::abs(c)});
    constexpr double eps = std::numeric_limits<double>::epsilon();
    if (!(magnitude > 256 * eps)) { ++diagnostics.degenerateAxisGroups; return {}; }
    a /= magnitude; b /= magnitude; c /= magnitude;
    double discriminant = std::fma(b, b, -4*a*c);
    double roundoff = 64 * eps * (b*b + std::abs(4*a*c));
    if (discriminant < -roundoff) { ++diagnostics.noRealAxisGroups; return {}; }
    discriminant = std::max(0.0, discriminant);
    double root = -.5 * (b + std::copysign(std::sqrt(discriminant), b));
    // Homogeneous stable quadratic formula also handles linear and double roots.
    std::array<std::array<double, 2>, 2> roots{{{a, root}, {root, c}}};
    std::vector<TorusAxis> axes;
    for (const auto &st : roots)
    {
        double length = std::hypot(st[0], st[1]);
        if (!(length > 0)) continue;
        double s = st[0] / length, t = st[1] / length;
        V direction = ud*s + vd*t, moment = um*s + vm*t;
        double norm = Norm(direction);
        if (!(norm > 128 * eps)) { ++diagnostics.rejectedAxes; continue; }
        direction = direction / norm; moment = moment / norm;
        if (!Finite(direction) || !Finite(moment) ||
            std::abs(Dot(direction, moment)) > 4096 * eps * (1 + Norm(moment)))
        { ++diagnostics.rejectedAxes; continue; }
        V position = Cross(direction, moment);
        bool valid = Finite(position);
        for (unsigned i = 0; i < 4; ++i)
            valid &= std::abs(Dot(points[i] - position, Cross(p[i].normal, direction))) <=
                8192 * eps * (1 + Norm(position));
        if (!valid) { ++diagnostics.rejectedAxes; continue; }
        bool duplicate = false;
        for (const auto &existing : axes)
            duplicate |= Norm(position - existing.position) < 1e-10 * (1 + Norm(position)) &&
                         Parallel(direction, existing.direction, 1e-10);
        if (!duplicate) axes.push_back({direction, position});
    }
    diagnostics.axisCount += static_cast<unsigned>(axes.size());
    return axes;
}
// Algebraic four-normal axis enumeration, followed by a conditioned meridian
// circle fit. Supplemental points come only from the generation partition.
std::vector<AnalyticParameters> TorusSeeds(const std::array<SurfaceSample, 4> &p, double scale,
    const std::vector<SurfaceSample> &support, GeneralTorusSeedDiagnostics &diagnostics)
{
    ++diagnostics.sampleGroups;
    V mean{};
    for (const auto &s : p) mean = mean + s.point * .25;
    std::array<V, 4> points;
    for (unsigned i = 0; i < 4; ++i) points[i] = (p[i].point - mean) / scale;
    std::vector<AnalyticParameters> out;
    for (const auto &axis : TorusAxes(points, p, diagnostics))
    {
        std::vector<std::array<double, 2>> meridian;
        auto add = [&](V point) {
            V d = point - axis.position;
            double z = Dot(d, axis.direction);
            meridian.push_back({Norm(d - axis.direction*z), z});
        };
        for (V point : points) add(point);
        for (const auto &s : support)
        {
            bool duplicate = false;
            for (const auto &seed : p) duplicate |= s.u == seed.u && s.v == seed.v;
            if (!duplicate) add((s.point - mean) / scale);
        }
        double rhoMean = 0, zMean = 0;
        for (const auto &x : meridian)
        { rhoMean += x[0] / meridian.size(); zMean += x[1] / meridian.size(); }
        double spread = 0;
        for (const auto &x : meridian) spread = std::max(spread, std::hypot(x[0]-rhoMean, x[1]-zMean));
        if (!(spread > Roundoff(std::hypot(rhoMean, zMean))))
        { ++diagnostics.circleFitFailures; continue; }
        std::vector<std::array<double, 3>> rows;
        std::vector<double> rhs;
        for (const auto &x : meridian)
        {
            double rho = (x[0] - rhoMean) / spread, z = (x[1] - zMean) / spread;
            rows.push_back({2*rho, 2*z, 1}); rhs.push_back(rho*rho + z*z);
        }
        std::array<double, 3> circle;
        if (!LeastSquares3(rows, rhs, 1e-12, circle))
        { ++diagnostics.circleFitFailures; continue; }
        double r2 = circle[2] + circle[0]*circle[0] + circle[1]*circle[1];
        double major = rhoMean + spread * circle[0];
        double minor = r2 > 0 ? spread * std::sqrt(r2) : 0;
        V center = mean + (axis.position + axis.direction * (zMean + spread * circle[1])) * scale;
        if (!(minor > 0 && major > minor) || !std::isfinite(major*scale) ||
            !std::isfinite(minor*scale) || !Finite(center))
        { ++diagnostics.invalidRadii; continue; }
        auto a = Frame(Type::Torus, center, axis.direction);
        a.radius = minor * scale; a.majorRadius = major * scale;
        out.push_back(a);
    }
    diagnostics.candidateCount += static_cast<unsigned>(out.size());
    diagnostics.dualCandidateGroups += out.size() == 2;
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
        context.CheckDeadline("classification_scores");
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
    ++diagnostics.seedGroups[candidate.seed.score.seedGroup].fullyScored;
    try
    {
        ScoreInterior(candidate, context, remaining, diagnostics);
        auto &score = candidate.seed.score;
        if (!(candidate.weightSum > 0)) throw Failure(Status::IllConditioned, "No regular independent score points");
        score.boundaryMax = 0;
        for (const auto &edge : context.Boundaries())
            for (const auto &s : edge)
            {
                context.CheckDeadline("classification_boundary_scores");
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
// Spread a small, deterministic support set over the generation partition.
// This neither samples the source again nor consumes the random attempt budget.
static std::vector<SurfaceSample> SeedSupport(const std::vector<SurfaceSample> &points,
                                              const std::vector<size_t> &indices)
{
    double u0 = points[indices.front()].u, u1 = u0, v0 = points[indices.front()].v, v1 = v0;
    for (size_t i : indices)
    {
        u0 = std::min(u0, points[i].u); u1 = std::max(u1, points[i].u);
        v0 = std::min(v0, points[i].v); v1 = std::max(v1, points[i].v);
    }
    double du = u1 > u0 ? u1-u0 : 1, dv = v1 > v0 ? v1-v0 : 1;
    std::vector<double> nearest(indices.size(), std::numeric_limits<double>::infinity());
    std::vector<SurfaceSample> support;
    size_t next = 0;
    while (support.size() < std::min(size_t(12), indices.size()))
    {
        const auto &chosen = points[indices[next]];
        support.push_back(chosen);
        for (size_t j = 0; j < indices.size(); ++j)
        {
            double u = (points[indices[j]].u-chosen.u)/du, v = (points[indices[j]].v-chosen.v)/dv;
            nearest[j] = std::min(nearest[j], u*u+v*v);
        }
        nearest[next] = -1;
        next = static_cast<size_t>(std::max_element(nearest.begin(), nearest.end())-nearest.begin());
    }
    return support;
}
GeneralCandidatePool InitializeGeneralCandidates(SurfaceSampleContext &context, const GeneralConversionOptions &o,
                                                 GeneralClassificationDiagnostics &diagnostics)
{
    context.CheckDeadline("classification");
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
    std::vector<SurfaceSample> support;
    unsigned familyIndex = 0;
    for (Type type : {Type::Plane, Type::Sphere, Type::Cylinder, Type::Cone, Type::Torus})
    {
        auto &family = pool.families[familyIndex++];
        for (unsigned attempt = 0; attempt < o.seedAttemptsPerType; ++attempt)
        {
            context.CheckDeadline("classification_seeds");
            ++diagnostics.seedAttemptCount;
            unsigned groupIndex=static_cast<unsigned>(diagnostics.seedGroups.size());
            diagnostics.seedGroups.emplace_back();
            auto &group=diagnostics.seedGroups.back();
            group.type=type;group.attempt=attempt+1;
            group.multiPoint=attempt==0 && (type==Type::Sphere || type==Type::Cylinder || type==Type::Cone);
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
            if (!complete && !group.multiPoint)
            {group.failure=GeneralSeedFailure::InsufficientPoints;continue;}
            try
            {
                std::vector<AnalyticParameters> seeds;
                std::vector<SurfaceSample> fit;
                if(group.multiPoint || type==Type::Torus)
                    if(support.empty()) support=SeedSupport(points,indices);
                if(type==Type::Sphere || type==Type::Cylinder || type==Type::Cone)
                    fit=group.multiPoint?support:std::vector<SurfaceSample>(p.begin(),p.end());
                group.sampleCount=fit.empty()?4:static_cast<unsigned>(fit.size());
                switch (type)
                {
                case Type::Plane: seeds.push_back(PlaneSeed(p, scale)); break;
                case Type::Sphere: case Type::Cylinder: case Type::Cone:
                    seeds.push_back(ConditionedGeneralSeed(type,fit,scale,group));break;
                default:
                    seeds = TorusSeeds(p, scale, support, diagnostics.torus);
                    break;
                }
                group.generated=static_cast<unsigned>(seeds.size());
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
                            ++group.duplicates;
                            continue;
                        }
                        GeneralScoringCandidate candidate;
                        candidate.seed.parameters = a;
                        candidate.seed.score.type = type;
                        candidate.seed.score.degreesOfFreedom = Dof(type);
                        candidate.seed.score.seedGroup = groupIndex;
                        // Undefined until CompleteScore evaluates the boundary.
                        candidate.seed.score.boundaryMax = std::numeric_limits<double>::infinity();
                        family.push_back(std::move(candidate));
                    }
                    catch (const std::runtime_error &) {group.failure=GeneralSeedFailure::ConstructionFailed;}
                }
            }
            catch (const std::runtime_error &)
            {if(group.failure==GeneralSeedFailure::None) group.failure=GeneralSeedFailure::ConstructionFailed;}
        }
    }
    // Build the same deduplicated mother pool in both modes before scoring.
    // Otherwise a representative failing only during full scoring can suppress
    // different near-duplicates depending on when that failure is discovered.
    for (auto &family : pool.families)
        for (auto &candidate : family)
        {
            context.CheckDeadline("classification_scores");
            ++diagnostics.scoreEvaluationCount;
            ++diagnostics.seedGroups[candidate.seed.score.seedGroup].scored;
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
    context.CheckDeadline("classification_batch");
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
            context.CheckDeadline("classification_scores");
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
            ++diagnostics.seedGroups[family[i].seed.score.seedGroup].selected;
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
