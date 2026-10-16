#include "GeneralRationalBoundaryInternal.h"
#include <set>
#include <sstream>

namespace bspline2analytic::detail::rational_boundary
{
namespace
{
bool Collapsed(const sggk::BSplineSurface &s,unsigned e)
{
    const auto &cp=s.ControlPoints();bool fixedU=e<2,high=e%2;
    const auto &k=fixedU?s.KnotsU():s.KnotsV();const auto &m=fixedU?s.MultsU():s.MultsV();
    unsigned degree=fixedU?s.DegreeU():s.DegreeV();
    double value=fixedU?(high?s.MaxParamU():s.MinParamU()):(high?s.MaxParamV():s.MinParamV());
    if(value!=(high?k.back():k.front())||(high?m.back():m.front())!=degree+1)return false;
    unsigned fixed=high?unsigned((fixedU?cp.size():cp.front().size())-1):0;
    V first=Get(fixedU?cp[fixed][0]:cp[0][fixed]);
    for(unsigned j=0;j<(fixedU?cp.front().size():cp.size());++j)
    {V p=Get(fixedU?cp[fixed][j]:cp[j][fixed]);if(p.x!=first.x||p.y!=first.y||p.z!=first.z)return false;}
    return true;
}
std::array<double,2> UV(const Trim &t,V p,std::array<double,2> reference)
{
    auto a=MotherParameters(t.mother,FrameOf(t.mother).local(ProjectToMother(t.mother,p)));
    for(unsigned d=0;d<2;++d)
        if((d==0&&t.mother.type!=sggk::SurfaceType::Plane)||(d==1&&t.mother.type==sggk::SurfaceType::Torus))
            a[d]=Unwrap(a[d],reference[d]);
    return a;
}
}
bool Recover(SurfaceSampleContext &samples,const AnalyticParameters &a,
             const GeneralRationalBoundaryOptions &options,Trim &out,RecoveryDiagnostics *diagnostics)
{
    if(diagnostics)*diagnostics={};
    auto reject=[&](const char *reason,int edge=-1,int sample=-1,double parameter=NAN,double value=NAN,double limit=NAN)
    {
        if(diagnostics){diagnostics->reason=reason;diagnostics->edge=edge;diagnostics->sample=sample;
            diagnostics->parameter=parameter;diagnostics->value=value;diagnostics->limit=limit;}
        return false;
    };
    // Reuse the established periodic-walk proposals, but never its analytic-UV
    // linearity gate or native PCurve construction. Acceptance is rational-domain.
    GeneralTrim trim;const auto &o=options.general;
    trim = {}; trim.parameters = a;
    auto frame = FrameOf(a);
    std::array<std::array<double, 2>, 4> cornerUV{};
    std::array<V, 4> cornerPoints;
    std::array<bool, 4> assigned{};
    const auto &s = samples.source;
    auto fail = [&](PCurveFailure reason, unsigned edge, double deviation = 0)->bool {
        reject(PCurveFailureName(reason),int(edge),-1,NAN,deviation);
        std::ostringstream message;
        message << PCurveFailureName(reason) << ": boundary " << edge << ", deviation=" << deviation;
        throw Failure(Status::VerificationFailed,message.str());
    };
    for (unsigned i = 0; i < 4; ++i)
        cornerPoints[i] = ProjectToMother(a, samples.At(i%2 ? s.MaxParamU() : s.MinParamU(),
            i/2 ? s.MaxParamV() : s.MinParamV(),sggk::DerivType::NONE,
            sggk::DerivType::NONE,false).point);
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        samples.CheckDeadline("chart_recovery");
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
        // At a pole only U is undefined. Extrapolate a LIMIT PROPOSAL using
        // nearby regular values; it does not assert analytic-UV linearity.
        // The rational-domain gate and exported-curve certificate decide acceptance.
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
    out={};out.mother=a;out.domain={INFINITY,-INFINITY,INFINITY,-INFINITY};
    for(const auto &walk:trim.walks)for(const auto &p:walk)for(unsigned d=0;d<2;++d)
    {
        out.domain[d*2]=std::min(out.domain[d*2],p.second[d]);
        out.domain[d*2+1]=std::max(out.domain[d*2+1],p.second[d]);
    }
    if(options.chartDomain)out.domain=*options.chartDomain;
    for(unsigned d=0;d<2;++d)
    {
        double lo=out.domain[2*d],hi=out.domain[2*d+1];
        if(!(std::isfinite(lo)&&std::isfinite(hi)&&lo<hi))return reject("InvalidChartExtent",int(d),-1,NAN,hi-lo);
        bool angular=d==0?a.type!=sggk::SurfaceType::Plane:
            (a.type==sggk::SurfaceType::Sphere||a.type==sggk::SurfaceType::Torus);
        if(angular&&hi-lo>2*Pi+1e-12)return reject("AngularExtentExceeded",int(d),-1,NAN,hi-lo,2*Pi+1e-12);
    }
    if(a.type==sggk::SurfaceType::Sphere&&(out.domain[2]<0||out.domain[3]>Pi))return reject("SphereLatitudeOutOfRange");
    if(a.type==sggk::SurfaceType::Cone&&out.domain[2]<0)return reject("NegativeConeHeight");
    AnalyticPatch patch;patch.parameters=a;patch.domain=out.domain;
    // This validator certifies boundary position/tangent, not surface normals.
    // Omit unused normal coefficients and their subsequent frame transforms.
    out.chart=paired::BuildChart(patch,false);paired::PlaceChart(out.chart,a,Get(a.origin));
    auto map=[&](std::array<double,2> uv)
    {
        // Extrapolation/clamping is only a recovery proposal. Gross violations
        // must not disappear into RationalCoordinate's endpoint clamps.
        for(unsigned d=0;d<2;++d)
            if(uv[d]<out.domain[2*d]-1e-12||uv[d]>out.domain[2*d+1]+1e-12)
            {
                reject("WalkOutsideChart",int(d),-1,NAN,uv[d]);
                throw Failure(Status::VerificationFailed,"Walk outside selected rational chart");
            }
        return std::array<double,2>{paired::RationalCoordinate(out.chart,0,uv[0]),
                                   paired::RationalCoordinate(out.chart,1,uv[1])};
    };
    for(unsigned i=0;i<4;++i)
    {
        auto q=map(connected[i]);out.polygon.corners[i]={q[0],q[1]};
    }
    if(diagnostics)for(unsigned e=0;e<4;++e)
        for(const auto &entry:trim.walks[e])diagnostics->walks[e].push_back({entry.first,map(entry.second)});
    // Certify a strictly convex parameter quadrilateral. A spatially collapsed
    // pole/apex edge can still have two distinct rational-chart endpoints.
    constexpr unsigned loop[]={0,1,3,2};int orientation=0;
    for(unsigned i=0;i<4;++i)
    {
        auto p=out.polygon.corners[loop[i]],q=out.polygon.corners[loop[(i+1)%4]],
             r=out.polygon.corners[loop[(i+2)%4]];
        I turn=(I(q.X())-I(p.X()))*(I(r.Y())-I(q.Y()))-
               (I(q.Y())-I(p.Y()))*(I(r.X())-I(q.X()));
        int sign=turn.lo>0?1:turn.hi<0?-1:0;
        if(!sign)return reject("UnresolvedConvexTurn",int(i),-1,NAN,turn.lo,turn.hi);
        if(orientation&&sign!=orientation)return reject("NonConvexCorners",int(i),-1,NAN,double(sign),double(orientation));
        orientation=sign;
    }
    if(TrimTriangles(out.polygon).empty())return reject("EmptyTrimTriangles");
    for(unsigned e=0;e<4;++e)
    {
        samples.CheckDeadline("chart_recovery");
        out.collapsed[e]=Collapsed(samples.source,e);
        for(const auto &entry:trim.walks[e])out.walks[e].push_back({entry.first,map(entry.second)});
        auto start=out.polygon.corners[Ends[e][0]],end=out.polygon.corners[Ends[e][1]];
        double dx=end.X()-start.X(),dy=end.Y()-start.Y(),length=dx*dx+dy*dy,last=-INFINITY;
        if(!(length>0))return reject("ZeroParameterEdge",int(e),-1,NAN,length);
        if(out.collapsed[e])continue;
        for(size_t i=0;i<out.walks[e].size();++i)
        {
            auto q=out.walks[e][i].second;
            double f=((q[0]-start.X())*dx+(q[1]-start.Y())*dy)/length;
            if(f<last-1e-10)return reject("NonMonotoneParameter",int(e),int(i),out.walks[e][i].first,f,last-1e-10);
            last=f;f=std::clamp(f,0.,1.);
            auto uv=trim.walks[e][i].second;
            V projected=FrameOf(a).global(CanonicalPoint(a,uv[0],uv[1]));
            double distance=Norm(projected-ChartPoint(out,{start.X()+dx*f,start.Y()+dy*f}));
            if(!(distance<o.fitting.validation.distanceTolerance))
                return reject("RationalLineWitnessExceeded",int(e),int(i),out.walks[e][i].first,distance,o.fitting.validation.distanceTolerance);
        }
    }
    if(diagnostics)diagnostics->reason="Recovered";
    return true;
}
std::array<double,2> Coordinates(const Trim &t,V p,std::array<double,2> reference)
{
    auto uv=UV(t,p,reference);
    return {paired::RationalCoordinate(t.chart,0,uv[0]),paired::RationalCoordinate(t.chart,1,uv[1])};
}
double Fraction(const Trim &t,unsigned e,V p,double reference)
{
    auto a=t.polygon.corners[Ends[e][0]],b=t.polygon.corners[Ends[e][1]];
    std::array<double,2> ref{};
    for(unsigned d=0;d<2;++d)ref[d]=(t.domain[2*d]+t.domain[2*d+1])*.5;
    auto q=Coordinates(t,p,ref);
    double dx=b.X()-a.X(),dy=b.Y()-a.Y(),length=dx*dx+dy*dy;
    return length>0?((q[0]-a.X())*dx+(q[1]-a.Y())*dy)/length:reference;
}
double FiniteDistance(const Trim &t,V p)
{
    std::array<double,2> ref{(t.domain[0]+t.domain[1])*.5,(t.domain[2]+t.domain[3])*.5};
    auto q=Coordinates(t,p,ref);
    if(InsideTrim(t.polygon,q,0))return Norm(p-ChartPoint(t,q));
    double best=INFINITY;
    for(const auto &tri:TrimTriangles(t.polygon))
    {
        double x=tri[1].X()-tri[0].X(),y=tri[1].Y()-tri[0].Y();
        double z=tri[2].X()-tri[0].X(),w=tri[2].Y()-tri[0].Y(),det=x*w-y*z;
        if(det==0)continue;
        double px=q[0]-tri[0].X(),py=q[1]-tri[0].Y();
        double b=(px*w-py*z)/det,c=(x*py-y*px)/det;
        double weight[]={std::max(0.,1-b-c),std::max(0.,b),std::max(0.,c)};
        double sum=weight[0]+weight[1]+weight[2];std::array<double,2> on{};
        for(unsigned i=0;i<3;++i){on[0]+=tri[i].X()*weight[i]/sum;on[1]+=tri[i].Y()*weight[i]/sum;}
        best=std::min(best,Norm(p-ChartPoint(t,on)));
    }
    return best;
}
}
