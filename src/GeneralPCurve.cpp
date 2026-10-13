#include "GeneralInternal.h"
#include <Geometry/2D/BSplineCurve2D.h>
#include <set>
#include <sstream>
namespace bspline2analytic
{
const char *PCurveFailureName(PCurveFailure f)
{
    switch (f)
    {
    case PCurveFailure::None: return "None";
    case PCurveFailure::NonlinearPCurve: return "NonlinearPCurve";
    case PCurveFailure::PCurveOutOfRange: return "PCurveOutOfRange";
    case PCurveFailure::NonMonotonePCurve: return "NonMonotonePCurve";
    case PCurveFailure::InconsistentPeriodicBranch: return "InconsistentPeriodicBranch";
    case PCurveFailure::UnresolvedDegeneratePCurve: return "UnresolvedDegeneratePCurve";
    case PCurveFailure::InvalidBoundaryLoop: return "InvalidBoundaryLoop";
    case PCurveFailure::BoundaryGeometryExceeded: return "BoundaryGeometryExceeded";
    case PCurveFailure::ValidationUnresolved: return "ValidationUnresolved";
    case PCurveFailure::ChartRecoveryUnresolved: return "ChartRecoveryUnresolved";
    }
    return "Unknown";
}
}
namespace bspline2analytic::detail
{
static constexpr unsigned Ends[4][2] = {{0,2}, {1,3}, {0,1}, {2,3}};
std::array<double, 2> EdgeUV(const LinearPCurve &e, double t)
{
    auto uv = e.CalcPoint(t); return {uv.X(), uv.Y()};
}
std::array<double, 2> TrimScale(const GeneralTrim &t)
{
    const auto &a = t.parameters;
    if (a.type == sggk::SurfaceType::Plane) return {1, 1};
    if (a.type == sggk::SurfaceType::Cylinder) return {a.radius, 1};
    if (a.type == sggk::SurfaceType::Sphere) return {a.radius, a.radius};
    if (a.type == sggk::SurfaceType::Torus) return {a.majorRadius + a.radius, a.radius};
    double height = 0;
    for (const auto &p : t.corners) height = std::max(height, std::abs(p.Y()));
    return {std::max(height * std::tan(a.halfAngle), 1e-12), 1 / std::cos(a.halfAngle)};
}
static double Cross2(std::array<double, 2> a, std::array<double, 2> b)
{ return a[0]*b[1] - a[1]*b[0]; }
std::vector<std::array<sggk::Point2D,3>> TrimTriangles(const GeneralTrim &t)
{
    auto area=[](const sggk::Point2D &a,const sggk::Point2D &b,const sggk::Point2D &c) {
        return Cross2({b.X()-a.X(),b.Y()-a.Y()},{c.X()-a.X(),c.Y()-a.Y()});
    };
    for (const auto &split: {std::array<unsigned,6>{0,1,3,0,3,2},std::array<unsigned,6>{0,1,2,1,3,2}})
    {
        double a=area(t.corners[split[0]],t.corners[split[1]],t.corners[split[2]]);
        double b=area(t.corners[split[3]],t.corners[split[4]],t.corners[split[5]]);
        if (a*b<0 || std::abs(a+b)<1e-14) continue;
        std::vector<std::array<sggk::Point2D,3>> result;
        if (std::abs(a)>1e-14) result.push_back({t.corners[split[0]],t.corners[split[1]],t.corners[split[2]]});
        if (std::abs(b)>1e-14) result.push_back({t.corners[split[3]],t.corners[split[4]],t.corners[split[5]]});
        return result;
    }
    return {};
}
bool InsideTrim(const GeneralTrim &t, std::array<double, 2> p, double tolerance)
{
    for (const auto &triangle:TrimTriangles(t))
    {
        double orientation=Cross2({triangle[1].X()-triangle[0].X(),triangle[1].Y()-triangle[0].Y()},
            {triangle[2].X()-triangle[0].X(),triangle[2].Y()-triangle[0].Y()});
        bool inside=true;
        for (unsigned i=0;i<3;++i)
        {
            const auto &a=triangle[i],&b=triangle[(i+1)%3];
            if ((orientation>0?1:-1)*Cross2({b.X()-a.X(),b.Y()-a.Y()},{p[0]-a.X(),p[1]-a.Y()})<-tolerance)
                inside=false;
        }
        if (inside) return true;
    }
    return false;
}
std::array<double, 2> NearestTrimUV(const GeneralTrim &t, std::array<double, 2> p)
{
    if (InsideTrim(t, p)) return p;
    auto scale = TrimScale(t);
    std::array<double, 2> best{};
    double error = INFINITY;
    for (const auto &edge : t.edges)
    {
        double x = (edge.end.X()-edge.start.X()) * scale[0];
        double y = (edge.end.Y()-edge.start.Y()) * scale[1];
        double px = (p[0]-edge.start.X())*scale[0], py = (p[1]-edge.start.Y())*scale[1];
        double d = x*x+y*y;
        double fraction = d > 0 ? std::clamp((px*x+py*y)/d, 0.0, 1.0) : 0;
        double e = std::hypot(px-x*fraction, py-y*fraction);
        if (e < error) { error = e; best = EdgeUV(edge, fraction); }
    }
    return best;
}
double EdgeReference(const GeneralTrim &t, unsigned edge, double fraction, unsigned coordinate)
{
    const auto &walk = t.walks[edge];
    auto next = std::lower_bound(walk.begin(), walk.end(), fraction,
        [](const auto &x, double f) { return x.first < f; });
    if (next == walk.begin()) return next->second[coordinate];
    if (next == walk.end()) return walk.back().second[coordinate];
    auto previous = next - 1;
    double k = (fraction - previous->first) / (next->first - previous->first);
    return previous->second[coordinate] + k * (next->second[coordinate] - previous->second[coordinate]);
}
static bool RecoverChart(SurfaceSampleContext &samples, const GeneralConversionOptions &o,
    const AnalyticParameters &a, GeneralTrim &trim, GeneralValidationReport &report)
{
    trim = {}; trim.parameters = a;
    auto frame = FrameOf(a);
    std::array<std::array<double, 2>, 4> cornerUV{};
    std::array<V, 4> cornerPoints;
    std::array<bool, 4> assigned{};
    const auto &s = samples.source;
    auto fail = [&](PCurveFailure reason, unsigned edge, double deviation = 0) {
        report.pcurveFailure = reason; report.failedBoundary = static_cast<BoundaryId>(edge);
        report.pcurveLinearity = CheckState::Exceeded;
        report.maxSampledPCurveDeviation = std::max(report.maxSampledPCurveDeviation, deviation);
        std::ostringstream message;
        message << PCurveFailureName(reason) << ": boundary " << edge << ", deviation=" << deviation;
        report.message = message.str(); return false;
    };
    for (unsigned i = 0; i < 4; ++i)
        cornerPoints[i] = ProjectToMother(a, samples.At(i%2 ? s.MaxParamU() : s.MinParamU(),
            i/2 ? s.MaxParamV() : s.MinParamV(),sggk::DerivType::NONE,
            sggk::DerivType::NONE,false).point);
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        bool fixedU = edge < 2;
        double lo = fixedU ? s.MinParamV() : s.MinParamU(), hi = fixedU ? s.MaxParamV() : s.MaxParamU();
        auto parameters = SampleParameters(fixedU ? s.KnotsV() : s.KnotsU(), lo, hi, o.boundarySamplesPerSpan);
        std::set<double> ts(parameters.begin(), parameters.end());
        for (const auto &x : samples.Boundaries()[edge]) ts.insert(fixedU ? x.v : x.u);
        std::vector<bool> regular;
        double world = 0, extent = 0;
        V first = cornerPoints[Ends[edge][0]];
        for (double t : ts)
        {
            const auto &x = samples.At(fixedU ? (edge ? s.MaxParamU() : s.MinParamU()) : t,
                fixedU ? t : (edge == 2 ? s.MinParamV() : s.MaxParamV()),
                sggk::DerivType::NONE,sggk::DerivType::NONE,false);
            V q = ProjectToMother(a, x.point), local = frame.local(q);
            world = std::max(world, Norm(q)); extent = std::max(extent, Norm(q-first));
            bool validU = a.type == sggk::SurfaceType::Plane ||
                std::hypot(local.x,local.y) > Roundoff(Norm(local));
            // A fitted pole/apex can miss an exact source corner by a few ulps
            // or by the fitting residual. Its azimuth is ill-conditioned. Use
            // the adjacent regular limit as a chart proposal; the continuous
            // 3D certificate, not this threshold, decides acceptance.
            // Include samples arbitrarily close to an endpoint (e.g. a knot
            // one ulp below it), not only parameters exactly equal to lo/hi.
            {
                double singularTolerance=std::max(Roundoff(world)*4,o.fitting.validation.distanceTolerance*.01);
                if((a.type==sggk::SurfaceType::Cone && Norm(local)<singularTolerance) ||
                   (a.type==sggk::SurfaceType::Sphere && std::hypot(local.x,local.y)<singularTolerance))
                    validU=false;
            }
            regular.push_back(validU);
            trim.walks[edge].push_back({(t-lo)/(hi-lo), MotherParameters(a,local)});
        }
        trim.edges[edge].collapsed3D = extent <= Roundoff(world) * 4;
        auto &walk = trim.walks[edge];
        size_t initial = 0;
        while (initial < regular.size() && !regular[initial]) ++initial;
        if (initial == regular.size())
        {
            if (!trim.edges[edge].collapsed3D)
                return fail(PCurveFailure::UnresolvedDegeneratePCurve, edge);
            continue;
        }
        double lastU = walk[initial].second[0], lastV = walk[initial].second[1];
        for (size_t j = 0; j < walk.size(); ++j)
        {
            if (regular[j])
            {
                if (a.type != sggk::SurfaceType::Plane) walk[j].second[0] = Unwrap(walk[j].second[0], lastU);
                lastU = walk[j].second[0];
            }
            else walk[j].second[0] = lastU;
            if (a.type == sggk::SurfaceType::Torus) walk[j].second[1] = Unwrap(walk[j].second[1], lastV);
            lastV = walk[j].second[1];
        }
        // At a pole only U is undefined. For a linear UV trajectory extrapolate
        // U against the defined polar/height coordinate using nearby regular UVs.
        for (size_t endpoint : {size_t(0),walk.size()-1})
            if (!regular[endpoint] && !trim.edges[edge].collapsed3D)
            {
                size_t nearby=regular.size(),second=regular.size();
                for (size_t k=0;k<regular.size();++k)
                {
                    size_t j=endpoint==0?k:regular.size()-1-k;
                    if (regular[j]) { if (nearby==regular.size()) nearby=j; else {second=j;break;} }
                }
                if (second==regular.size()) return fail(PCurveFailure::UnresolvedDegeneratePCurve,edge);
                double dv=walk[second].second[1]-walk[nearby].second[1];
                if (std::abs(dv)<1e-12) return fail(PCurveFailure::UnresolvedDegeneratePCurve,edge);
                walk[endpoint].second[0]=walk[nearby].second[0]+
                    (walk[second].second[0]-walk[nearby].second[0])*
                    (walk[endpoint].second[1]-walk[nearby].second[1])/dv;
            }
        if (!trim.edges[edge].collapsed3D)
            for (unsigned endpoint = 0; endpoint < 2; ++endpoint)
            {
                unsigned corner = Ends[edge][endpoint];
                auto uv = endpoint ? walk.back().second : walk.front().second;
                if (!assigned[corner]) { cornerUV[corner] = uv; assigned[corner] = true; }
            }
    }
    for (bool valid : assigned)
        if (!valid) return fail(PCurveFailure::UnresolvedDegeneratePCurve, 0);
    // Connect in perimeter order using whole-period translations only. Degenerate
    // sides use neighbouring limits and the opposite side's winding as a guide.
    constexpr unsigned order[] = {2,1,3,0};
    std::array<std::array<double,2>,4> connected{};
    connected[0] = cornerUV[0];
    for (unsigned k = 0; k < 4; ++k)
    {
        unsigned edge = order[k]; bool reversed = k >= 2;
        unsigned from = Ends[edge][reversed ? 1 : 0], to = Ends[edge][reversed ? 0 : 1];
        auto &walk = trim.walks[edge];
        auto first = reversed ? walk.back().second : walk.front().second;
        auto last = reversed ? walk.front().second : walk.back().second;
        std::array<double,2> destination{};
        for (unsigned d = 0; d < 2; ++d)
        {
            bool periodic = d == 0 ? a.type != sggk::SurfaceType::Plane : a.type == sggk::SurfaceType::Torus;
            if (trim.edges[edge].collapsed3D)
            {
                destination[d] = cornerUV[to][d];
                if (periodic)
                {
                    unsigned opposite = edge ^ 1;
                    const auto &other = trim.walks[opposite];
                    double delta = other.back().second[d] - other.front().second[d];
                    destination[d] = Unwrap(destination[d], connected[from][d] + (reversed ? -delta : delta));
                }
            }
            else
            {
                double shift = connected[from][d] - first[d];
                if (periodic && std::abs(shift - 2*Pi*std::round(shift/(2*Pi))) > 1e-8)
                    return fail(PCurveFailure::InconsistentPeriodicBranch, edge);
                if (!periodic && std::abs(shift) > 1e-8)
                    return fail(PCurveFailure::InconsistentPeriodicBranch, edge);
                destination[d] = last[d] + shift;
                for (auto &entry : walk) entry.second[d] += shift;
            }
        }
        if (to == 0)
        {
            if (std::hypot(destination[0]-connected[0][0], destination[1]-connected[0][1]) > 1e-8)
                return fail(PCurveFailure::InconsistentPeriodicBranch, edge);
            destination=connected[0];
        }
        else connected[to] = destination;
        auto start = reversed ? destination : connected[from];
        auto end = reversed ? connected[from] : destination;
        trim.edges[edge].start = {start[0],start[1]}; trim.edges[edge].end = {end[0],end[1]};
        trim.edges[edge].projectedStart = Point(cornerPoints[Ends[edge][0]]);
        trim.edges[edge].projectedEnd = Point(cornerPoints[Ends[edge][1]]);
        if (trim.edges[edge].collapsed3D)
            for (auto &entry : walk) entry.second = EdgeUV(trim.edges[edge],entry.first);
    }
    for (unsigned i = 0; i < 4; ++i) trim.corners[i] = {connected[i][0],connected[i][1]};
    constexpr unsigned loop[] = {0,1,3,2};
    double sign = 0, area = 0;
    for (unsigned i = 0; i < 4; ++i)
    {
        auto p = connected[loop[i]], q = connected[loop[(i+1)%4]], r = connected[loop[(i+2)%4]];
        double turn = Cross2({q[0]-p[0],q[1]-p[1]}, {r[0]-q[0],r[1]-q[1]});
        if (std::abs(turn)<=1e-14 &&
            (q[0]-p[0])*(r[0]-q[0])+(q[1]-p[1])*(r[1]-q[1])<-1e-14)
            return fail(PCurveFailure::InvalidBoundaryLoop,0);
        if (std::abs(turn) > 1e-14)
        {
            sign = turn;
        }
        area += p[0]*q[1]-p[1]*q[0];
    }
    if (std::abs(area) < 1e-14 || sign == 0 || TrimTriangles(trim).empty())
        return fail(PCurveFailure::InvalidBoundaryLoop,0);
    for (unsigned d=0;d<2;++d)
    {
        double lo=INFINITY,hi=-INFINITY;
        for (const auto &uv:connected) {lo=std::min(lo,uv[d]);hi=std::max(hi,uv[d]);}
        bool periodic=d==0?a.type!=sggk::SurfaceType::Plane:a.type==sggk::SurfaceType::Torus;
        if (periodic && hi-lo>2*Pi+1e-10) return fail(PCurveFailure::InvalidBoundaryLoop,0);
        if (d==1 && ((a.type==sggk::SurfaceType::Sphere && (lo<-1e-10 || hi>Pi+1e-10)) ||
                     (a.type==sggk::SurfaceType::Cone && lo<-1e-10)))
            return fail(PCurveFailure::InvalidBoundaryLoop,0);
    }
    auto scale = TrimScale(trim);
    AnalyticPatch mother; mother.parameters = a; MakeNative(mother);
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        auto &line = trim.edges[edge];
        double dx = (line.end.X()-line.start.X())*scale[0], dy = (line.end.Y()-line.start.Y())*scale[1];
        double length2 = dx*dx+dy*dy, last = -INFINITY;
        for (const auto &entry : trim.walks[edge])
        {
            auto uv = entry.second;
            double x = (uv[0]-line.start.X())*scale[0], y = (uv[1]-line.start.Y())*scale[1];
            double t = length2 > 0 ? (x*dx+y*dy)/length2 : entry.first;
            double margin = o.fitting.validation.distanceTolerance / std::max(std::sqrt(length2),1e-12);
            if (t < -margin || t > 1+margin) return fail(PCurveFailure::PCurveOutOfRange,edge);
            if (t < last-margin) return fail(PCurveFailure::NonMonotonePCurve,edge);
            last = t;
            auto onLine = EdgeUV(line,std::clamp(t,0.0,1.0));
            double deviation = Norm(CanonicalPoint(a,uv[0],uv[1])-CanonicalPoint(a,onLine[0],onLine[1]));
            report.maxSampledPCurveDeviation = std::max(report.maxSampledPCurveDeviation,deviation);
            if (deviation >= o.fitting.validation.distanceTolerance)
                return fail(PCurveFailure::NonlinearPCurve,edge,deviation);
        }
        auto sdkAt = [&](double t, std::array<double,2> reference) {
            auto uv = EdgeUV(line,t);
            V point = frame.global(CanonicalPoint(a,uv[0],uv[1]));
            auto native = mother.surface->CalcParam(Point(point));
            // Singular angle is fixed by an interior chart probe, not arbitrary
            // CalcParam output at a sphere pole or cone apex.
            V local=frame.local(point);
            if (std::hypot(local.x,local.y)<=std::max(Roundoff(Norm(point))*4,
                    o.fitting.validation.distanceTolerance*.01) &&
                (a.type == sggk::SurfaceType::Sphere || a.type == sggk::SurfaceType::Cone))
            {
                auto probe = mother.surface->CalcParam(Point(frame.global(CanonicalPoint(a,uv[0],.8))));
                native = {probe.X(),native.Y()};
            }
            return std::array<double,2>{a.type == sggk::SurfaceType::Plane ? native.X() : Unwrap(native.X(),reference[0]),
                a.type == sggk::SurfaceType::Torus ? Unwrap(native.Y(),reference[1]) : native.Y()};
        };
        auto start = sdkAt(0,{line.start.X(),line.start.Y()}), end = start;
        for (unsigned step = 1; step <= 64; ++step) end = sdkAt(double(step)/64,end);
        for (double t : {0.0,.25,.5,.75,1.0})
        {
            auto uv = EdgeUV(line,t);
            V expected = frame.global(CanonicalPoint(a,uv[0],uv[1]));
            V native = Get(mother.surface->CalcPoint(start[0]+(end[0]-start[0])*t,start[1]+(end[1]-start[1])*t));
            if (Norm(expected-native) > std::max(Roundoff(Norm(expected))*16,1e-9))
                return fail(PCurveFailure::InconsistentPeriodicBranch,edge,Norm(expected-native));
        }
        line.sdkCurve = std::make_shared<sggk::BSplineCurve2D>(1,
            sggk::Point2DArray{{start[0],start[1]},{end[0],end[1]}}, sggk::RealArray{0,1}, sggk::UIntArray{2,2});
    }
    // Share SDK-UV vertex records as well as canonical vertices. Translate an
    // entire PCurve by whole periods before assigning its endpoint, then check
    // the closing edge; never independently snap unrelated 3D-equal corners.
    std::array<sggk::Point2D,4> nativeCorners;
    nativeCorners[0]=trim.edges[2].sdkCurve->CalcPoint(0);
    for (unsigned k=0;k<4;++k)
    {
        unsigned edge=order[k];bool reversed=k>=2;
        unsigned from=Ends[edge][reversed?1:0],to=Ends[edge][reversed?0:1];
        auto first=trim.edges[edge].sdkCurve->CalcPoint(reversed?1:0);
        auto last=trim.edges[edge].sdkCurve->CalcPoint(reversed?0:1);
        std::array<double,2> shift{nativeCorners[from].X()-first.X(),nativeCorners[from].Y()-first.Y()};
        for (unsigned d=0;d<2;++d)
        {
            bool periodic=d==0?a.type!=sggk::SurfaceType::Plane:a.type==sggk::SurfaceType::Torus;
            double residual=periodic?shift[d]-2*Pi*std::round(shift[d]/(2*Pi)):shift[d];
            if (std::abs(residual)>1e-8) return fail(PCurveFailure::InconsistentPeriodicBranch,edge);
        }
        sggk::Point2D destination(last.X()+shift[0],last.Y()+shift[1]);
        if (to==0)
        {
            if (std::hypot(destination.X()-nativeCorners[0].X(),destination.Y()-nativeCorners[0].Y())>1e-8)
                return fail(PCurveFailure::InconsistentPeriodicBranch,edge);
            destination=nativeCorners[0];
        }
        else nativeCorners[to]=destination;
        sggk::Point2D start=reversed?destination:nativeCorners[from];
        sggk::Point2D end=reversed?nativeCorners[from]:destination;
        trim.edges[edge].sdkCurve=std::make_shared<sggk::BSplineCurve2D>(1,
            sggk::Point2DArray{start,end},sggk::RealArray{0,1},sggk::UIntArray{2,2});
    }
    report.pcurveLinearity = CheckState::Passed; report.pcurveFailure = PCurveFailure::None;
    return true;
}
bool RecoverLinearPCurves(SurfaceSampleContext &samples, const GeneralConversionOptions &o,
    const AnalyticParameters &a, GeneralTrim &trim, GeneralValidationReport &report,
    const std::function<bool(const GeneralTrim &,GeneralValidationReport &)> &accept)
{
    std::vector<AnalyticParameters> charts{a};
    if (a.type == sggk::SurfaceType::Sphere)
    {
        V center = Get(a.origin);
        auto add = [&](V axis) {
            if (Norm(axis) <= Roundoff(Norm(center))) return;
            AnalyticParameters chart = a; chart.axis = Direction(axis);
            chart.xDirection = Direction(Perpendicular(Unit(axis))); charts.push_back(chart);
        };
        add({1,0,0}); add({0,1,0});
        for (const auto &edge : samples.Boundaries())
        {
            if (edge.size() < 3) continue;
            add(Cross(edge[edge.size()/2].point-edge.front().point, edge.back().point-edge.front().point));
            if (Norm(edge.back().point-edge.front().point) <= Roundoff(Norm(edge.front().point)))
                add(edge.front().point-center);
        }
    }
    bool recovered=false;
    unsigned attempts=0,iterations=0;
    GeneralValidationReport continuous;
    auto attemptChart=[&](const AnalyticParameters &chart)
    {
        ++attempts;
        GeneralValidationReport attempt;
        try
        {
            if (RecoverChart(samples,o,chart,trim,attempt))
            {
                recovered=true;
                if(!accept || accept(trim,attempt))
                {report=attempt;report.chartAttempts=attempts;report.chartIterations=iterations;return true;}
                continuous=attempt;
            }
        }
        catch (const std::exception &e)
        {
            attempt.pcurveFailure = PCurveFailure::UnresolvedDegeneratePCurve;
            attempt.pcurveLinearity = CheckState::Unresolved; attempt.message = e.what();
        }
        report = attempt;
        return false;
    };
    for(const auto &chart:charts) if(attemptChart(chart)) return true;
    if(a.type==sggk::SurfaceType::Sphere)
    {
        for(const auto &chart:RefinedSphereCharts(samples,o,a,charts,iterations)) if(attemptChart(chart)) return true;
        if(!recovered)
        {
            report.pcurveLinearity=CheckState::Unresolved;
            report.pcurveFailure=PCurveFailure::ChartRecoveryUnresolved;
            report.message="Sphere chart recovery exhausted; no global nonlinearity claim";
        }
    }
    if(recovered) report=continuous;
    report.chartAttempts=attempts;report.chartIterations=iterations;
    return false;
}
}
