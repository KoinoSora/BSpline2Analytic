#include "PairedValidationInternal.h"
#include "StageTimer.h"
#include "ValidationTiming.h"
#include <functional>

namespace bspline2analytic::detail::paired
{
namespace
{
using Clock = std::chrono::steady_clock;
double Milliseconds(Clock::time_point t)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}
struct Budget : std::runtime_error
{
    using std::runtime_error::runtime_error;
};
I MapValue(const PairedParameterMap &m, double t)
{
    auto it = std::lower_bound(m.source.begin(), m.source.end(), t);
    if (it != m.source.end() && *it == t)
        return I(m.target[it - m.source.begin()]);
    if (it == m.source.begin() || it == m.source.end())
        throw Failure(Status::InvalidInput, "Mapping outside source domain");
    size_t i = static_cast<size_t>(it - m.source.begin() - 1);
    return I(m.target[i]) + (I(t) - I(m.source[i])) / (I(m.source[i + 1]) - I(m.source[i])) *
                                (I(m.target[i + 1]) - I(m.target[i]));
}
void CheckMap(const PairedParameterMap &m)
{
    if (m.source.size() < 2 || m.source.size() != m.target.size() || m.target.front() != 0 ||
        m.target.back() != 1)
        throw Failure(Status::VerificationFailed, "Invalid mapping endpoints");
    for (size_t i = 1; i < m.source.size(); ++i)
        if (!(std::isfinite(m.source[i - 1]) && std::isfinite(m.source[i]) &&
              (I(m.source[i]) - I(m.source[i - 1])).lo > 0 &&
              (I(m.target[i]) - I(m.target[i - 1])).lo > 0))
            throw Failure(Status::VerificationFailed, "Mapping is not strictly monotone");
}
void AddChartBreaks(PairedParameterMap &m, const Chart &chart, bool u)
{
    std::vector<double> breaks;
    for (const auto &c : chart.cells)
    {
        breaks.push_back(u ? c.surface.u0 : c.surface.v0);
        breaks.push_back(u ? c.surface.u1 : c.surface.v1);
    }
    std::sort(breaks.begin(), breaks.end());
    breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
    for (double q : breaks)
    {
        auto it = std::lower_bound(m.target.begin(), m.target.end(), q);
        if (it != m.target.end() && *it == q)
            continue;
        if (it == m.target.begin() || it == m.target.end())
            throw Failure(Status::VerificationFailed, "Chart break outside mapping");
        size_t i = static_cast<size_t>(it - m.target.begin() - 1);
        // These doubles DEFINE the new map; no approximate inverse is treated
        // as an exact certificate. Subsequent interval composition uses it.
        double t = m.source[i] + (m.source[i + 1] - m.source[i]) * (q - m.target[i]) /
                                     (m.target[i + 1] - m.target[i]);
        m.source.insert(m.source.begin() + i + 1, t);
        m.target.insert(m.target.begin() + i + 1, q);
    }
    CheckMap(m);
}
PairedParameterMap AdaptiveMap(const sggk::BSplineSurface &source, const AnalyticPatch &patch,
                               const Chart &chart, bool alongU,
                               const PairedValidationOptions &options,
                               PairedValidationReport &report, SurfaceSampleContext &samples,
                               const ClassificationSamples &cache, const PatchRecoveryTrace &trace)
{
    unsigned axis = alongU ? 0 : 1, dim = patch.mapping.swappedUV ? 1 - axis : axis;
    bool reverse = alongU ? patch.mapping.reversedU : patch.mapping.reversedV;
    double start = cache.sourceDomain[2 * axis], end = cache.sourceDomain[2 * axis + 1];
    auto frame = FrameOf(patch.parameters);
    unsigned first = alongU ? 2 : 0;
    const auto &a = patch.parameters;
    double speed = std::max({1.0, std::abs(a.radius), std::abs(a.majorRadius),
                             std::abs(patch.domain[1] - patch.domain[0]),
                             std::abs(patch.domain[3] - patch.domain[2])}) *
                   2;
    const double targetError = options.validation.distanceTolerance * .2 / speed;
    std::map<double, double> values;
    values[start] = 0;
    values[end] = 1;
    auto value = [&](double t) {
        auto found = values.find(t);
        if (found != values.end())
            return found->second;
        double sum = 0;
        for (unsigned e = first; e < first + 2; ++e)
        {
            double lo = trace.paths[e].front()[dim], hi = trace.paths[e].back()[dim];
            if (!std::isfinite(hi - lo) || std::abs(hi - lo) < 1e-14)
                throw Failure(Status::Unsupported, "Collapsed parameter trace");
            double q;
            const auto &ts = cache.boundaryParameters[e];
            auto saved = std::lower_bound(ts.begin(), ts.end(), t);
            if (saved != ts.end() && *saved == t)
                q = trace.paths[e][saved - ts.begin()][dim];
            else
            {
                auto x = samples
                             .At(alongU ? t : (e == 0 ? source.MinParamU() : source.MaxParamU()),
                                 alongU ? (e == 2 ? source.MinParamV() : source.MaxParamV()) : t,
                                 sggk::DerivType::NONE, sggk::DerivType::NONE, false)
                             .point;
                ++report.mappingPointEvaluations;
                q = MotherParameters(a, frame.local(x))[dim];
            }
            if (chart.angular[dim] && (dim == 0 || a.type == sggk::SurfaceType::Torus))
                q = Unwrap(q, lo + (hi - lo) * (t - start) / (end - start));
            double f = (q - lo) / (hi - lo);
            double angle = reverse ? patch.domain[2 * dim + 1] -
                                         f * (patch.domain[2 * dim + 1] - patch.domain[2 * dim])
                                   : patch.domain[2 * dim] +
                                         f * (patch.domain[2 * dim + 1] - patch.domain[2 * dim]);
            double r = RationalCoordinate(chart, dim, angle);
            sum += reverse ? 1 - r : r;
        }
        double q = sum * .5;
        values[t] = q;
        return q;
    };
    PairedParameterMap map;
    map.source.push_back(start);
    map.target.push_back(0);
    std::function<void(double, double, double, double, unsigned)> refine;
    refine = [&](double l, double r, double a0, double a1, unsigned depth) {
        if (map.source.size() > options.maxMappingSegments)
            throw Budget("Mapping segment budget exhausted");
        double error = 0;
        for (double f : {.25, .5, .75})
        {
            double t = l + (r - l) * f;
            if (!(l < t && t < r))
                throw Failure(Status::VerificationFailed, "Mapping parameter resolution exhausted");
            error = std::max(error, std::abs(value(t) - (a0 + (a1 - a0) * f)));
        }
        if (error > targetError)
        {
            if (depth >= options.maxMappingDepth)
                throw Budget("Mapping depth exhausted");
            double mid = l + (r - l) * .5, q = value(mid);
            if (!(a0 < q && q < a1))
                throw Failure(Status::VerificationFailed, "Opposite-edge map is not monotone");
            refine(l, mid, a0, q, depth + 1);
            refine(mid, r, q, a1, depth + 1);
        }
        else
        {
            map.source.push_back(r);
            map.target.push_back(a1);
        }
    };
    refine(start, end, 0, 1, 0);
    CheckMap(map);
    return map;
}
std::vector<double> Breaks(const sggk::RealArray &knots, const PairedParameterMap &map)
{
    std::vector<double> out = map.source;
    for (double x : knots)
        if (x > out.front() && x < map.source.back())
            out.push_back(x);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}
TargetCell Compose(const TargetCell &target, const PairedParameterMap &u,
                   const PairedParameterMap &v, double u0, double u1, double v0, double v1)
{
    TargetCell out = target;
    auto axis = [&](bool along, double l, double r, double a, double b,
                    const PairedParameterMap &m) {
        I x = (MapValue(m, l) - I(a)) / (I(b) - I(a)), y = (MapValue(m, r) - I(a)) / (I(b) - I(a));
        // Preserve exact chart endpoints and avoid unnecessary interval growth.
        if (MapValue(m, l).lo == a && MapValue(m, l).hi == a)
            x = I(0);
        if (MapValue(m, r).lo == b && MapValue(m, r).hi == b)
            y = I(1);
        out.surface.p = Restrict(out.surface.p, along, x, y);
        out.surface.w = Restrict(out.surface.w, along, x, y);
        out.normal = Restrict(out.normal, along, x, y);
    };
    axis(true, u0, u1, target.surface.u0, target.surface.u1, u);
    axis(false, v0, v1, target.surface.v0, target.surface.v1, v);
    out.surface.u0 = u0;
    out.surface.u1 = u1;
    out.surface.v0 = v0;
    out.surface.v1 = v1;
    return out;
}
I Length(IV x)
{
    return Sqrt(Sqr(x.x) + Sqr(x.y) + Sqr(x.z));
}
double ErrorUpper(const BV &e, const BPoly &weight)
{
    unsigned u = std::max({e.x.u, e.y.u, e.z.u, weight.u}),
             v = std::max({e.x.v, e.y.v, e.z.v, weight.v});
    BPoly zero(u, v), x = e.x + zero, y = e.y + zero, z = e.z + zero, w = weight + zero;
    double upper = 0;
    for (size_t i = 0; i < w.c.size(); ++i)
    {
        if (w.c[i].lo <= 0)
            return INFINITY;
        upper = std::max(upper, Length({x.c[i] / w.c[i], y.c[i] / w.c[i], z.c[i] / w.c[i]}).hi);
    }
    return upper;
}
struct AngleBounds
{
    double ratio = INFINITY;
    bool exceeded = false;
};
AngleBounds Angle(const BV &a, const BV &b, I threshold, Products &products)
{
    BPoly dot = products.Dot(a, b);
    BV cross = products.Cross(a, b);
    I d = dot.range();
    double ratio = d.lo > 0 ? (Length(Range(cross)) / I(d.lo)).hi : INFINITY;
    if (ratio < threshold.lo)
        return {ratio, false};
    auto ap = Range(Evaluate(a, .5, 0)), bp = Range(Evaluate(b, .5, 0));
    I pointDot = detail::Dot(ap, bp), pointCross = Length(detail::Cross(ap, bp));
    bool exceeded =
        pointDot.hi < 0 || (pointDot.lo > 0 && (pointCross / pointDot).lo >= threshold.hi);
    return {ratio, exceeded};
}
BV Half(const BV &p, bool u, bool right)
{
    auto x = Split(p.x, u), y = Split(p.y, u), z = Split(p.z, u);
    return right ? BV{x.second, y.second, z.second} : BV{x.first, y.first, z.first};
}
BPoly Half(const BPoly &p, bool u, bool right)
{
    auto h = Split(p, u);
    return right ? h.second : h.first;
}
} // namespace

class Validator
{
    const sggk::BSplineSurface &source;
    const AnalyticPatch &patch;
    const PairedValidationOptions &options;
    PairedValidationReport &out;
    ValidationTimer &timer;
    Products products;
    std::vector<BezierCell> sourceCells;
    Chart chart;
    double distanceUpper = 0, tangentRatio = 0, normalRatio = 0;
    I tangentLimit, normalLimit;
    bool geometryExceeded = false, tangentExceeded = false, normalExceeded = false;
    void Visit(bool boundary)
    {
        if (++out.certificate.piecesVisited > options.validation.maxValidationPieces)
            throw Budget("Validation node budget exhausted");
        if (boundary)
            ++out.boundaryNodes;
        else
            ++out.surfaceNodes;
    }
    bool Geometry(const BV &error, const BPoly &weight, unsigned depth)
    {
        Visit(false);
        double upper;
        {
            ValidationTimer::Scope stage(timer, out.certificate.timing.surfaceBounds);
            upper = ErrorUpper(error, weight);
            if (upper < options.validation.distanceTolerance)
            {
                distanceUpper = std::max(distanceUpper, upper);
                return true;
            }
            IV e = Range(Evaluate(error, .5, .5));
            I w = Evaluate(weight, .5, .5).range();
            if (w.lo > 0 && (Length(e) / w).lo >= options.validation.distanceTolerance)
            {
                geometryExceeded = true;
                distanceUpper = std::max(distanceUpper, upper);
                return false;
            }
        }
        if (depth >= options.validation.maxSubdivisionDepth)
        {
            distanceUpper = INFINITY;
            return false;
        }
        ValidationTimer::Scope stage(timer, out.certificate.timing.subdivision);
        bool u = depth % 2 == 0;
        auto left = Half(error, u, false), right = Half(error, u, true);
        auto wl = Half(weight, u, false), wr = Half(weight, u, true);
        return Geometry(left, wl, depth + 1) && Geometry(right, wr, depth + 1);
    }
    bool Edge(const BV &t, const BV &normal, const BV &targetT, const BV &targetN, unsigned depth)
    {
        Visit(true);
        auto a = Angle(t, targetT, tangentLimit, products),
             b = Angle(normal, targetN, normalLimit, products);
        if (a.ratio < tangentLimit.lo && b.ratio < normalLimit.lo)
        {
            tangentRatio = std::max(tangentRatio, a.ratio);
            normalRatio = std::max(normalRatio, b.ratio);
            return true;
        }
        if (a.exceeded || b.exceeded)
        {
            tangentExceeded |= a.exceeded;
            normalExceeded |= b.exceeded;
            return false;
        }
        if (depth >= options.validation.maxSubdivisionDepth)
            return false;
        return Edge(Half(t, true, false), Half(normal, true, false), Half(targetT, true, false),
                    Half(targetN, true, false), depth + 1) &&
               Edge(Half(t, true, true), Half(normal, true, true), Half(targetT, true, true),
                    Half(targetN, true, true), depth + 1);
    }
    bool Attempt()
    {
        ++out.mappingAttempts;
        distanceUpper = tangentRatio = normalRatio = 0;
        geometryExceeded = tangentExceeded = normalExceeded = false;
        auto &report = out.certificate;
        report.geometry = {};
        report.boundaryTangent = {};
        report.boundaryNormal = {};
        report.correspondence = CheckState::Unresolved;
        AddChartBreaks(out.uMap, chart, true);
        AddChartBreaks(out.vMap, chart, false);
        if (out.uMap.source.size() - 1 > options.maxMappingSegments ||
            out.vMap.source.size() - 1 > options.maxMappingSegments)
            throw Budget("Mapping chart-break budget exhausted");
        auto us = Breaks(source.KnotsU(), out.uMap), vs = Breaks(source.KnotsV(), out.vMap);
        size_t nu = us.size() - 1, nv = vs.size() - 1;
        if (nu > options.maxCommonCells || nv > options.maxCommonCells / nu)
            throw Budget("Common grid budget exhausted");
        out.commonCells += static_cast<unsigned>(nu * nv);
        bool geometryOK = true, edgesOK = true;
        for (size_t i = 0; i < nu && geometryOK && edgesOK; ++i)
            for (size_t j = 0; j < nv && geometryOK && edgesOK; ++j)
            {
                double u0 = us[i], u1 = us[i + 1], v0 = vs[j], v1 = vs[j + 1];
                auto src =
                    std::find_if(sourceCells.begin(), sourceCells.end(), [&](const BezierCell &c) {
                        return u0 >= c.u0 && u1 <= c.u1 && v0 >= c.v0 && v1 <= c.v1;
                    });
                if (src == sourceCells.end())
                    throw Failure(Status::VerificationFailed, "Source partition coverage failure");
                I mu = MapValue(out.uMap, u0 + (u1 - u0) * .5),
                  mv = MapValue(out.vMap, v0 + (v1 - v0) * .5);
                auto dst =
                    std::find_if(chart.cells.begin(), chart.cells.end(), [&](const TargetCell &c) {
                        return mu.lo > c.surface.u0 && mu.hi < c.surface.u1 &&
                               mv.lo > c.surface.v0 && mv.hi < c.surface.v1;
                    });
                if (dst == chart.cells.end())
                    throw Failure(Status::VerificationFailed, "Target partition coverage failure");
                BezierCell s;
                TargetCell r;
                BV error;
                {
                    ValidationTimer::Scope stage(timer, report.timing.surfaceBounds);
                    s = Crop(*src, u0, u1, v0, v1);
                    r = Compose(*dst, out.uMap, out.vMap, u0, u1, v0, v1);
                    error = products.Mul(s.p, r.surface.w) - r.surface.p;
                }
                geometryOK = Geometry(error, r.surface.w, 0);
                if (!geometryOK)
                    break;
                for (unsigned edge = 0; edge < 4 && edgesOK; ++edge)
                {
                    bool outer = edge == 0   ? i == 0
                                 : edge == 1 ? i + 1 == nu
                                 : edge == 2 ? j == 0
                                             : j + 1 == nv;
                    if (!outer)
                        continue;
                    ValidationTimer::Scope stage(timer, report.timing.boundaryChecks);
                    auto id = static_cast<BoundaryId>(edge);
                    BV tangent = Boundary(Derivative(s.p, edge >= 2), id);
                    BV normal = products.Cross(Boundary(Derivative(s.p, true), id),
                                               Boundary(Derivative(s.p, false), id));
                    BV rp = Boundary(r.surface.p, id);
                    // Boundary packs the varying direction into the first axis.
                    BV rwVector{r.surface.w, r.surface.w, r.surface.w};
                    BPoly rw = Boundary(rwVector, id).x;
                    BV targetT = products.Mul(Derivative(rp, true), rw) -
                                 products.Mul(rp, Derivative(rw, true));
                    edgesOK = Edge(tangent, normal, targetT, Boundary(r.normal, id), 0);
                }
            }
        ValidationTimer::Scope stage(timer, report.timing.reporting);
        report.geometry.state = geometryOK         ? CheckState::Passed
                                : geometryExceeded ? CheckState::Exceeded
                                                   : CheckState::Unresolved;
        // Early exit means unvisited cells are not bounded by the partial maximum.
        report.geometry.upperBound = geometryOK && edgesOK ? distanceUpper : INFINITY;
        if (geometryOK && edgesOK)
        {
            report.boundaryTangent.state = report.boundaryNormal.state = report.correspondence =
                CheckState::Passed;
            report.boundaryTangent.upperBound = (Atan(I(tangentRatio)) * I(180) / IPi()).hi;
            report.boundaryNormal.upperBound = (Atan(I(normalRatio)) * I(180) / IPi()).hi;
            return report.boundaryTangent.upperBound <
                       options.validation.boundaryTangentToleranceDeg &&
                   report.boundaryNormal.upperBound < options.validation.boundaryNormalToleranceDeg;
        }
        report.boundaryTangent.state =
            tangentExceeded ? CheckState::Exceeded : CheckState::Unresolved;
        report.boundaryNormal.state =
            normalExceeded ? CheckState::Exceeded : CheckState::Unresolved;
        if (!edgesOK)
            report.geometry.state = CheckState::Unresolved;
        return false;
    }

  public:
    Validator(const sggk::BSplineSurface &s, const AnalyticPatch &p,
              const PairedValidationOptions &o, PairedValidationReport &r, ValidationTimer &t)
        : source(s), patch(p), options(o), out(r), timer(t)
    {
    }
    bool Run()
    {
        if (!IntervalEnvironmentSupported())
            throw Failure(Status::VerificationFailed,
                          "Unsupported interval floating-point environment");
        CheckApproximateInput(source, options.validation);
        CheckPatch(patch);
        if (!options.maxMappingSegments || options.maxMappingSegments > 4096 ||
            options.maxMappingDepth > 20 || !options.maxCommonCells ||
            options.maxCommonCells > 1000000)
            throw Failure(Status::InvalidInput, "Invalid paired validation budget");
        if (source.DegreeU() > 32 || source.DegreeV() > 32)
            throw Failure(Status::Unsupported, "Paired degree budget is 32 per direction");
        if (source.KnotsU().size() - 1 > options.maxCommonCells / (source.KnotsV().size() - 1))
            throw Budget("Source cell budget exhausted");
        sourceCells = ExtractCells(source);
        ++out.sourceExtractions;
        if (sourceCells.empty())
            throw Failure(Status::InvalidInput, "Empty source domain");
        for (auto &c : sourceCells)
            c.p = LocalSource(c.p, patch.parameters);
        auto start = Clock::now();
        chart = BuildChart(patch);
        out.nurbsMs += Milliseconds(start);
        out.targetPatches = static_cast<unsigned>(chart.cells.size());
        I t = I(options.validation.boundaryTangentToleranceDeg) * IPi() / I(180),
          n = I(options.validation.boundaryNormalToleranceDeg) * IPi() / I(180);
        tangentLimit = Sin(t) / Cos(t);
        normalLimit = Sin(n) / Cos(n);
        out.uMap = {{source.MinParamU(), source.MaxParamU()}, {0, 1}};
        out.vMap = {{source.MinParamV(), source.MaxParamV()}, {0, 1}};
        if (Attempt())
            return true;
        if (!options.enablePiecewiseMapping)
            return false;
        timer.Switch(out.certificate.timing.preparation);
        start = Clock::now();
        try
        {
            SurfaceSampleContext samples(source, 50000);
            auto cache = BuildClassificationSamples(source, nullptr, &samples);
            auto trace = TraceApproximateBoundaries(cache, patch.parameters);
            out.uMap = AdaptiveMap(source, patch, chart, true, options, out, samples, cache, trace);
            out.vMap =
                AdaptiveMap(source, patch, chart, false, options, out, samples, cache, trace);
        }
        catch (...)
        {
            out.mappingMs += Milliseconds(start);
            throw;
        }
        out.mappingMs += Milliseconds(start);
        return Attempt();
    }
};
} // namespace bspline2analytic::detail::paired

namespace bspline2analytic
{
const char *PairedValidationStatusName(PairedValidationStatus s) noexcept
{
    switch (s)
    {
    case PairedValidationStatus::Passed:
        return "Passed";
    case PairedValidationStatus::NotCertified:
        return "NotCertified";
    case PairedValidationStatus::Unsupported:
        return "Unsupported";
    case PairedValidationStatus::InvalidInput:
        return "InvalidInput";
    case PairedValidationStatus::BudgetExceeded:
        return "BudgetExceeded";
    default:
        return "NumericalFailure";
    }
}
PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &s,
                                                const AnalyticPatch &p,
                                                const PairedValidationOptions &o)
{
    PairedValidationReport out;
    detail::ValidationTimer timer(out.certificate.timing);
    try
    {
        detail::paired::Validator validator(s, p, o, out, timer);
        bool passed = validator.Run();
        out.status = passed ? PairedValidationStatus::Passed : PairedValidationStatus::NotCertified;
        out.certificate.message = passed ? "Positive-weight paired rational certificate: complete "
                                           "surface, oriented boundary tangents and normals"
                                         : "The chosen correspondence was not certified; geometric "
                                           "infeasibility is not established";
    }
    catch (const detail::paired::Budget &e)
    {
        out.status = PairedValidationStatus::BudgetExceeded;
        out.certificate.message = e.what();
    }
    catch (const detail::Failure &e)
    {
        out.status = e.status == Status::Unsupported    ? PairedValidationStatus::Unsupported
                     : e.status == Status::InvalidInput ? PairedValidationStatus::InvalidInput
                     : e.status == Status::Degenerate   ? PairedValidationStatus::Unsupported
                                                        : PairedValidationStatus::NotCertified;
        out.certificate.message = e.what();
    }
    catch (const std::exception &e)
    {
        out.status = PairedValidationStatus::NumericalFailure;
        out.certificate.message = e.what();
    }
    timer.Switch(out.certificate.timing.reporting);
    if (out.status != PairedValidationStatus::Passed)
    {
        out.certificate.correspondence = CheckState::Unresolved;
        // No partial domain may masquerade as a complete bound after early exit.
        out.certificate.geometry.upperBound = out.certificate.boundaryTangent.upperBound =
            out.certificate.boundaryNormal.upperBound = INFINITY;
    }
    timer.Finish();
    return out;
}
} // namespace bspline2analytic
