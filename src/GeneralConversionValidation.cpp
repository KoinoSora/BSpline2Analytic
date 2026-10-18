/**
 * @file GeneralConversionValidation.cpp
 * @brief 历史解析 UV 有限裁剪域的投影与边界验证。
 * 内部位置、边界覆盖和 PCurve 线性约束均属于该历史合同，不可套用到当前有理边界输出。
 */
#include "GeneralInternal.h"
#include <functional>
namespace bspline2analytic::detail
{
static I Abs(I a)
{
    return {a.lo <= 0 && a.hi >= 0 ? 0 : std::min(std::abs(a.lo), std::abs(a.hi)), a.absMax()};
}
struct MotherBound
{
    I distance = I::Whole();
    bool regularProjection = false;
};
static MotherBound BoundMother(const BezierCell &cell, const AnalyticParameters &a)
{
    I weight = cell.w.range();
    if (!(weight.lo > 0)) return {};
    IV axis = Normalize(Exact(Get(a.axis)));
    V origin = Get(a.origin);
    BV constant{BPoly(0,0,I(origin.x)),BPoly(0,0,I(origin.y)),BPoly(0,0,I(origin.z))};
    BV d = cell.p - constant * cell.w;
    BPoly z = Dot(d, axis);
    BV radial = d - BV{z * axis.x, z * axis.y, z * axis.z};
    BPoly radial2 = Dot(radial, radial), norm2 = Dot(d, d), ww = cell.w * cell.w;
    I rho = Sqrt(radial2.range()) / weight;
    I height = z.range() / weight;
    I radius(a.radius);
    auto residualBound = [&](const BPoly &f, I denominator) {
        I error = Abs(f.range()) / denominator;
        return I(0, error.hi);
    };
    switch (a.type)
    {
    case sggk::SurfaceType::Plane:
        return {Abs(height), true};
    case sggk::SurfaceType::Sphere:
        return {residualBound(norm2 - ww * Sqr(radius), Sqr(I(weight.lo)) * radius),
                norm2.range().lo > 0};
    case sggk::SurfaceType::Cylinder:
        return {residualBound(radial2 - ww * Sqr(radius), Sqr(I(weight.lo)) * radius),
                radial2.range().lo > 0};
    case sggk::SurfaceType::Cone: {
        I c = Cos(I(a.halfAngle)), s = Sin(I(a.halfAngle));
        I along = rho * s + height * c;
        if (along.hi <= 0)
            return {Sqrt(norm2.range()) / weight, true};
        I direct = Abs(rho * c - height * s);
        I denominator = rho * c + height * s;
        if (height.lo >= 0 && denominator.lo > 0)
            direct = residualBound(radial2 * Sqr(c) - (z * z) * Sqr(s),
                                   Sqr(I(weight.lo)) * I(denominator.lo));
        if (along.lo < 0)
            direct = {0, std::max(direct.hi, (Sqrt(norm2.range()) / weight).hi)};
        return {direct, radial2.range().lo > 0 || along.hi <= 0};
    }
    case sggk::SurfaceType::Torus: {
        I R(a.majorRadius);
        BPoly sum = norm2 + ww * (Sqr(R) - Sqr(radius));
        BPoly f = sum * sum - (ww * radial2) * (I(4) * Sqr(R));
        I factor = Sqr(rho + R) + Sqr(height) - Sqr(radius);
        I tube = Sqr(rho - R) + Sqr(height);
        if (!(factor.lo > 0)) return {};
        return {residualBound(f, Sqr(Sqr(I(weight.lo))) * I(factor.lo) * radius),
                radial2.range().lo > 0 && tube.lo > 0};
    }
    default: return {};
    }
}
// 有理 Bézier 单元及其对单元归一化参数的导数所需的区间几何数据。
struct Field { IV value, derivative; };
static IV Divide(IV a, I b) { return a * (I(1)/b); }
static Field SourceField(const BezierCell &cell, bool alongU, bool midpoint)
{
    BV derivative = Derivative(cell.p,alongU);
    BPoly dw = Derivative(cell.w,alongU);
    BV numerator = midpoint ? Evaluate(cell.p,.5,.5) : cell.p;
    BPoly w = midpoint ? Evaluate(cell.w,.5,.5) : cell.w;
    I weight = w.range();
    return {Divide(Range(numerator),weight),
        Divide(Range(derivative*cell.w-cell.p*dw),Sqr(cell.w.range()))};
}
static Field LocalField(Field field, const AnalyticParameters &a)
{
    IV z = Normalize(Exact(Get(a.axis))), x = Exact(Get(a.xDirection));
    x = Normalize(x-z*Dot(x,z)); IV y = Cross(z,x);
    IV p = field.value-Exact(Get(a.origin));
    return {{Dot(p,x),Dot(p,y),Dot(p,z)},
            {Dot(field.derivative,x),Dot(field.derivative,y),Dot(field.derivative,z)}};
}
static I Hypot(IV p) { return Sqrt(Sqr(p.x)+Sqr(p.y)+Sqr(p.z)); }
static I Clamp01(I t) { return {std::clamp(t.lo,0.0,1.0),std::clamp(t.hi,0.0,1.0)}; }
static I Angle(I y, I x, double reference)
{
    if (!x.finite() || !y.finite() || (x.lo <= 0 && x.hi >= 0 && y.lo <= 0 && y.hi >= 0))
        return I::Whole();
    I out;
    out.lo=INFINITY;out.hi=-INFINITY;
    for (double xx : {x.lo,x.hi}) for (double yy : {y.lo,y.hi})
    {
        I angle;
        if (xx == 0) angle = IPi()*I(yy >= 0 ? .5 : -.5);
        else
        {
            angle = Atan(I(yy)/I(xx));
            if (xx < 0) angle = angle + IPi()*I(yy >= 0 ? 1 : -1);
        }
        double raw = std::atan2(yy,xx);
        angle = angle + IPi()*I(2*std::round((reference-raw)/(2*Pi)));
        out.lo = std::min(out.lo,angle.lo); out.hi = std::max(out.hi,angle.hi);
    }
    return out.hi-out.lo < Pi ? out : I::Whole();
}
static std::array<I,2> Chart(Field p, const AnalyticParameters &a, std::array<double,2> reference)
{
    if (a.type == sggk::SurfaceType::Plane) return {p.value.x,p.value.y};
    I u = Angle(p.value.y,p.value.x,reference[0]);
    I rho = Sqrt(Sqr(p.value.x)+Sqr(p.value.y));
    if (a.type == sggk::SurfaceType::Cylinder) return {u,p.value.z};
    if (a.type == sggk::SurfaceType::Sphere)
        return {u,Acos(-p.value.z/Hypot(p.value))};
    if (a.type == sggk::SurfaceType::Torus)
        return {u,Angle(p.value.z,rho-I(a.majorRadius),reference[1])};
    I c = Cos(I(a.halfAngle)), sn = Sin(I(a.halfAngle));
    I along = rho*sn+p.value.z*c;
    return {u,I(std::max(0.0,along.lo),std::max(0.0,along.hi))*c};
}
static IV MotherPoint(const AnalyticParameters &a, I u, I v)
{
    I c = Cos(u), s = Sin(u), R(a.radius);
    switch (a.type)
    {
    case sggk::SurfaceType::Plane: return {u,v,I(0)};
    case sggk::SurfaceType::Cylinder: return {R*c,R*s,v};
    case sggk::SurfaceType::Sphere: return {R*Sin(v)*c,R*Sin(v)*s,-R*Cos(v)};
    case sggk::SurfaceType::Cone: {
        I tangent = Sin(I(a.halfAngle))/Cos(I(a.halfAngle));
        return {v*tangent*c,v*tangent*s,v};
    }
    default: {
        I r = I(a.majorRadius)+R*Cos(v); return {r*c,r*s,R*Sin(v)};
    }
    }
}
static IV MotherTangent(const AnalyticParameters &a, I u, I v, I du, I dv)
{
    I c=Cos(u),s=Sin(u),R(a.radius);
    IV pu,pv;
    switch (a.type)
    {
    case sggk::SurfaceType::Plane: pu={I(1),I(0),I(0)};pv={I(0),I(1),I(0)};break;
    case sggk::SurfaceType::Cylinder: pu={-R*s,R*c,I(0)};pv={I(0),I(0),I(1)};break;
    case sggk::SurfaceType::Sphere:
        pu={-R*Sin(v)*s,R*Sin(v)*c,I(0)};pv={R*Cos(v)*c,R*Cos(v)*s,R*Sin(v)};break;
    case sggk::SurfaceType::Cone: {
        I t=Sin(I(a.halfAngle))/Cos(I(a.halfAngle));
        pu={-v*t*s,v*t*c,I(0)};pv={t*c,t*s,I(1)};break;
    }
    default: {
        I r=I(a.majorRadius)+R*Cos(v);
        pu={-r*s,r*c,I(0)};pv={-R*Sin(v)*c,-R*Sin(v)*s,R*Cos(v)};break;
    }
    }
    return pu*du+pv*dv;
}
static Field ProjectField(Field p, const AnalyticParameters &a)
{
    if (a.type == sggk::SurfaceType::Plane)
        return {{p.value.x,p.value.y,I(0)},{p.derivative.x,p.derivative.y,I(0)}};
    if (a.type == sggk::SurfaceType::Sphere)
    {
        I r=Hypot(p.value), dr=Dot(p.value,p.derivative)/r;
        return {p.value*(I(a.radius)/r),
            (Divide(p.derivative,r)-p.value*(dr/Sqr(r)))*I(a.radius)};
    }
    I rho=Sqrt(Sqr(p.value.x)+Sqr(p.value.y));
    I dr=(p.value.x*p.derivative.x+p.value.y*p.derivative.y)/rho;
    I radial, radialDerivative, z=p.value.z, dz=p.derivative.z;
    if (a.type == sggk::SurfaceType::Cylinder)
    { radial=I(a.radius);radialDerivative=I(0); }
    else if (a.type == sggk::SurfaceType::Cone)
    {
        I c=Cos(I(a.halfAngle)),s=Sin(I(a.halfAngle));
        I t=rho*s+z*c,dt=dr*s+dz*c;
        if (t.hi <= 0) return {{I(0),I(0),I(0)},{I(0),I(0),I(0)}};
        if (t.lo < 0) dt=I(std::min(0.0,dt.lo),std::max(0.0,dt.hi));
        t=I(std::max(0.0,t.lo),std::max(0.0,t.hi));
        radial=t*s;radialDerivative=dt*s;z=t*c;dz=dt*c;
    }
    else
    {
        I x=rho-I(a.majorRadius),r=Sqrt(Sqr(x)+Sqr(z));
        I derivative=(x*dr+z*dz)/r;
        radial=I(a.majorRadius)+I(a.radius)*x/r;
        radialDerivative=I(a.radius)*(dr/r-x*derivative/Sqr(r));
        dz=I(a.radius)*(dz/r-z*derivative/Sqr(r));z=I(a.radius)*z/r;
    }
    I factor=radial/rho,df=radialDerivative/rho-radial*dr/Sqr(rho);
    return {{p.value.x*factor,p.value.y*factor,z},
            {p.derivative.x*factor+p.value.x*df,p.derivative.y*factor+p.value.y*df,dz}};
}
struct EdgeParameter { I value, derivative; };
static EdgeParameter MapEdge(Field p, const GeneralTrim &trim, unsigned edge,
                             std::array<double,2> reference, I sourceFraction, I fractionDerivative)
{
    const auto &e=trim.edges[edge];const auto &a=trim.parameters;
    double du=e.end.X()-e.start.X(),dv=e.end.Y()-e.start.Y();
    if (e.collapsed3D) return {sourceFraction,fractionDerivative};
    auto uv=Chart(p,a,reference);
    I rho=Sqrt(Sqr(p.value.x)+Sqr(p.value.y));
    I dr=(p.value.x*p.derivative.x+p.value.y*p.derivative.y)/rho;
    I ud=(p.value.x*p.derivative.y-p.value.y*p.derivative.x)/Sqr(rho),vd;
    if (a.type == sggk::SurfaceType::Plane)
    { ud=p.derivative.x;vd=p.derivative.y; }
    else if (a.type == sggk::SurfaceType::Cylinder) vd=p.derivative.z;
    else if (a.type == sggk::SurfaceType::Sphere)
        vd=(rho*p.derivative.z-p.value.z*dr)/(Sqr(rho)+Sqr(p.value.z));
    else if (a.type == sggk::SurfaceType::Torus)
        vd=((rho-I(a.majorRadius))*p.derivative.z-p.value.z*dr)/
            (Sqr(rho-I(a.majorRadius))+Sqr(p.value.z));
    else
    {
        I c=Cos(I(a.halfAngle)),s=Sin(I(a.halfAngle));vd=Sqr(c)*p.derivative.z+s*c*dr;
    }
    // 固定 U 的经线在球极点处仍有正则的圆参数；使用带符号径向坐标，避免在极点除以 rho。
    if (std::abs(du) < 1e-12 && std::abs(dv) > 1e-12 &&
        a.type != sggk::SurfaceType::Plane && a.type != sggk::SurfaceType::Cylinder)
    {
        I c=Cos(I(e.start.X())),s=Sin(I(e.start.X()));
        I w=p.value.x*c+p.value.y*s,dw=p.derivative.x*c+p.derivative.y*s;
        if (a.type == sggk::SurfaceType::Sphere)
        { uv[1]=Angle(w,-p.value.z,reference[1]);vd=(w*p.derivative.z-p.value.z*dw)/(Sqr(w)+Sqr(p.value.z)); }
        else if (a.type == sggk::SurfaceType::Torus)
        { uv[1]=Angle(p.value.z,w-I(a.majorRadius),reference[1]);vd=((w-I(a.majorRadius))*p.derivative.z-p.value.z*dw)/(Sqr(w-I(a.majorRadius))+Sqr(p.value.z)); }
        else
        { I ca=Cos(I(a.halfAngle)),sa=Sin(I(a.halfAngle));uv[1]=Sqr(ca)*p.value.z+sa*ca*w;vd=Sqr(ca)*p.derivative.z+sa*ca*dw; }
        return {(uv[1]-I(e.start.Y()))/I(dv),vd/I(dv)};
    }
    if (std::abs(dv) < 1e-12 && std::abs(du) > 1e-12)
        return {(uv[0]-I(e.start.X()))/I(du),ud/I(du)};
    if (std::abs(du)<1e-12 && std::abs(dv)>1e-12)
        return {(uv[1]-I(e.start.Y()))/I(dv),vd/I(dv)};
    auto scale=TrimScale(trim);
    I x(du*scale[0]*scale[0]),y(dv*scale[1]*scale[1]);
    I denominator=Sqr(I(du*scale[0]))+Sqr(I(dv*scale[1]));
    return {((uv[0]-I(e.start.X()))*x+(uv[1]-I(e.start.Y()))*y)/denominator,
            (ud*x+vd*y)/denominator};
}
static bool ValidateBoundaryCell(const BezierCell &cell, unsigned edge, const GeneralTrim &trim,
    const sggk::BSplineSurface &s, GeneralValidationReport &report, double tolerance, unsigned depth,
    const ValidationOptions &options)
{
    if (++report.piecesVisited > options.maxValidationPieces) return false;
    bool alongU=edge>=2;
    double lo=alongU?s.MinParamU():s.MinParamV(),hi=alongU?s.MaxParamU():s.MaxParamV();
    double begin=alongU?cell.u0:cell.v0,end=alongU?cell.u1:cell.v1;
    I fraction((begin-lo)/(hi-lo),(end-lo)/(hi-lo));
    double middle=(begin+end)*.5;
    std::array<double,2> reference{EdgeReference(trim,edge,(middle-lo)/(hi-lo),0),
                                 EdgeReference(trim,edge,(middle-lo)/(hi-lo),1)};
    auto source=LocalField(SourceField(cell,alongU,false),trim.parameters);
    auto mid=LocalField(SourceField(cell,alongU,true),trim.parameters);
    auto map=MapEdge(source,trim,edge,reference,fraction,I((end-begin)/(hi-lo)));
    auto centerMap=MapEdge(mid,trim,edge,reference,I((middle-lo)/(hi-lo)),I(0));
    const auto &e=trim.edges[edge];
    I du(e.end.X()-e.start.X()),dv(e.end.Y()-e.start.Y());
    I t=Clamp01(map.value),t0=Clamp01(centerMap.value);
    I dt=map.derivative;
    if (map.value.hi<=0 || map.value.lo>=1) dt=I(0);
    else if (map.value.lo<0 || map.value.hi>1) dt=I(std::min(0.0,dt.lo),std::max(0.0,dt.hi));
    I u=I(e.start.X())+du*t,v=I(e.start.Y())+dv*t;
    IV target=MotherPoint(trim.parameters,I(e.start.X())+du*t0,I(e.start.Y())+dv*t0);
    IV derivative=MotherTangent(trim.parameters,u,v,du*dt,dv*dt);
    I distance=Hypot(mid.value-target)+Hypot(source.derivative-derivative)*I(.5);
    auto projection=ProjectField(source,trim.parameters),midProjection=ProjectField(mid,trim.parameters);
    I deviation=Hypot(midProjection.value-target)+Hypot(projection.derivative-derivative)*I(.5);
    I midDistance=Hypot(mid.value-target),midDeviation=Hypot(midProjection.value-target);
    if (midDistance.lo>=tolerance || midDeviation.lo>=tolerance)
    {
        report.boundaryHausdorff=CheckState::Exceeded;
        report.pcurveFailure=midDeviation.lo>=tolerance ? PCurveFailure::NonlinearPCurve : PCurveFailure::BoundaryGeometryExceeded;
        report.pcurveLinearity=midDeviation.lo>=tolerance ? CheckState::Exceeded : report.pcurveLinearity;
        report.failedBoundary=static_cast<BoundaryId>(edge);
        report.hasWitness=true;report.witness={(cell.u0+cell.u1)*.5,(cell.v0+cell.v1)*.5};
        report.message=std::string(PCurveFailureName(report.pcurveFailure))+": boundary "+std::to_string(edge);
        return false;
    }
    // 奇异端点处可改用直接区间包络求界。
    I direct=Hypot(source.value-MotherPoint(trim.parameters,u,v));
    if (direct.finite() && (!distance.finite() || direct.hi<distance.hi)) distance=direct;
    I projectionFallback=distance+BoundMother(cell,trim.parameters).distance;
    if (projectionFallback.finite() && (!deviation.finite() || projectionFallback.hi<deviation.hi))
        deviation=projectionFallback;
    if (e.collapsed3D && trim.parameters.type==sggk::SurfaceType::Cone)
        deviation=Hypot(source.value)+Hypot(MotherPoint(trim.parameters,u,v));
    bool mapRegular=map.value.finite() && centerMap.value.finite();
    // 若任意源/目标点对距离均足够小，则微小边即使参数图坐标无定义也能保证覆盖；否则需正则连续坐标，以介值定理论证覆盖。
    bool allPairs=direct.finite() && direct.hi<tolerance;
    if ((mapRegular || allPairs) && distance.finite() && deviation.finite() &&
        distance.hi<tolerance && deviation.hi<tolerance)
    {
        report.boundaryHausdorffUpperBound=std::max(report.boundaryHausdorffUpperBound,distance.hi);
        return true;
    }
    report.hasWitness=true;report.witness={(cell.u0+cell.u1)*.5,(cell.v0+cell.v1)*.5};
    report.failedBoundary=static_cast<BoundaryId>(edge);
    if (depth>=options.maxSubdivisionDepth) return false;
    auto halves=Split(cell,alongU);
    return ValidateBoundaryCell(halves.first,edge,trim,s,report,tolerance,depth+1,options) &&
           ValidateBoundaryCell(halves.second,edge,trim,s,report,tolerance,depth+1,options);
}
static bool CertifiedInside(const std::array<sggk::Point2D,3> &triangle, std::array<I,2> uv)
{
    double area=(triangle[1].X()-triangle[0].X())*(triangle[2].Y()-triangle[0].Y())-
        (triangle[1].Y()-triangle[0].Y())*(triangle[2].X()-triangle[0].X());
    for (unsigned i=0;i<3;++i)
    {
        const auto &a=triangle[i],&b=triangle[(i+1)%3];
        I cross=(I(b.X())-I(a.X()))*(uv[1]-I(a.Y()))-(I(b.Y())-I(a.Y()))*(uv[0]-I(a.X()));
        if (area>0 ? cross.lo<0 : cross.hi>0) return false;
    }
    return true;
}
static I EndpointCoordinate(const BezierCell &cell,unsigned edge,bool start,const GeneralTrim &trim)
{
    bool alongU=edge>=2;
    BezierCell point=cell;
    double t=start?0:1;
    point.p=Evaluate(cell.p,alongU?t:.5,alongU?.5:t);
    point.w=Evaluate(cell.w,alongU?t:.5,alongU?.5:t);
    auto uv=EdgeUV(trim.edges[edge],t);
    return Clamp01(MapEdge(LocalField(SourceField(point,alongU,true),trim.parameters),
        trim,edge,uv,I(t),I(0)).value);
}
static double FiniteDistance(const BezierCell &cell, const GeneralTrim &trim, I motherDistance)
{
    double uRef=0,vRef=0;
    for (const auto &p:trim.corners) { uRef+=p.X()*.25;vRef+=p.Y()*.25; }
    auto uv=Chart(LocalField(SourceField(cell,true,false),trim.parameters),trim.parameters,{uRef,vRef});
    auto fallback=[&] {
        Field mid=LocalField(SourceField(cell,true,true),trim.parameters);
        V p{(mid.value.x.lo+mid.value.x.hi)*.5,(mid.value.y.lo+mid.value.y.hi)*.5,
            (mid.value.z.lo+mid.value.z.hi)*.5};
        auto raw=MotherParameters(trim.parameters,p);
        if (trim.parameters.type!=sggk::SurfaceType::Plane) raw[0]=Unwrap(raw[0],uRef);
        if (trim.parameters.type==sggk::SurfaceType::Torus) raw[1]=Unwrap(raw[1],vRef);
        auto target=NearestTrimUV(trim,raw);
        return Hypot(LocalField(SourceField(cell,true,false),trim.parameters).value-
            MotherPoint(trim.parameters,I(target[0]),I(target[1]))).hi;
    };
    if (!uv[0].finite() || !uv[1].finite()) return fallback();
    auto scale=TrimScale(trim);
    if (trim.parameters.type==sggk::SurfaceType::Cone)
    {
        I angle(trim.parameters.halfAngle);
        scale[0]=std::max(scale[0],(I(uv[1].absMax())*Sin(angle)/Cos(angle)).hi);
        scale[1]=(I(1)/Cos(angle)).hi;
    }
    for (auto &x:scale) x=(I(x)+I(Roundoff(x))).hi;
    double outside=INFINITY;
    // 每个三角形均为凸集。
    // 分别对整个 UV 矩形求到两个三角形的距离界，再取较小值作为到其并集的上界；因此也能处理简单凹四边形，无需误当作凸多边形。
    for (const auto &triangle:TrimTriangles(trim))
    {
        double triangleBound=0;
        for (double u:{uv[0].lo,uv[0].hi}) for (double v:{uv[1].lo,uv[1].hi})
        {
            if (CertifiedInside(triangle,{I(u),I(v)})) continue;
            double best=INFINITY;
            for (unsigned edge=0;edge<3;++edge)
            {
                const auto &a=triangle[edge],&b=triangle[(edge+1)%3];
                double dx=(b.X()-a.X())*scale[0],dy=(b.Y()-a.Y())*scale[1];
                double denominator=dx*dx+dy*dy;
                double t=denominator>0 ? std::clamp(((u-a.X())*scale[0]*dx+(v-a.Y())*scale[1]*dy)/denominator,0.0,1.0):0;
                I x=I(a.X())+(I(b.X())-I(a.X()))*I(t);
                I y=I(a.Y())+(I(b.Y())-I(a.Y()))*I(t);
                I d=Sqrt(Sqr((I(u)-x)*I(scale[0]))+Sqr((I(v)-y)*I(scale[1])));
                best=std::min(best,d.hi);
            }
            triangleBound=std::max(triangleBound,best);
        }
        outside=std::min(outside,triangleBound);
    }
    return std::min((motherDistance+Sqrt(I(2))*I(outside)).hi,fallback());
}
static GeneralValidationReport ValidateTrim(const sggk::BSplineSurface &s, const GeneralTrim &trim,
                                        const GeneralConversionOptions &options, SurfaceSampleContext &samples,
                                        unsigned consumed)
{
    GeneralValidationReport report;
    report.pcurveLinearity=CheckState::Passed;
    report.piecesVisited=consumed;
    const auto &o=options.fitting.validation;
    if (!(o.distanceTolerance>0 && std::isfinite(o.distanceTolerance)) ||
        o.maxSubdivisionDepth>30 || !o.maxValidationPieces ||
        options.boundarySamplesPerSpan<3 || options.boundarySamplesPerSpan>64)
        throw Failure(Status::InvalidInput,"Invalid general PCurve validation options");
    const auto &parameters=trim.parameters;
    if (!IntervalEnvironmentSupported())
    {
        report.geometry = report.boundaryHausdorff = CheckState::Unresolved;
        report.message = "Directed interval environment is unavailable";
        return report;
    }
    const auto cells = ExtractCells(s);
    if(TryPairedTrimBound(s,cells,trim,o,report))
        return report;
    report.boundaryHausdorffUpperBound=0;
    for (unsigned edge=0;edge<4;++edge)
    {
        double previousBound=report.boundaryHausdorffUpperBound;
        report.boundaryHausdorffUpperBound=0;
        I startCoordinate=I::Whole(),endCoordinate=I::Whole();
        for (const auto &cell:cells)
        {
            bool onEdge=edge==0 ? cell.u0==s.MinParamU() : edge==1 ? cell.u1==s.MaxParamU() :
                        edge==2 ? cell.v0==s.MinParamV() : cell.v1==s.MaxParamV();
            if (!onEdge) continue;
            BezierCell boundary=BoundaryCell(cell,static_cast<BoundaryId>(edge));
            double begin=edge<2?cell.v0:cell.u0,end=edge<2?cell.v1:cell.u1;
            if (begin==(edge<2?s.MinParamV():s.MinParamU()))
                startCoordinate=EndpointCoordinate(boundary,edge,true,trim);
            if (end==(edge<2?s.MaxParamV():s.MaxParamU()))
                endCoordinate=EndpointCoordinate(boundary,edge,false,trim);
            if (!ValidateBoundaryCell(boundary,edge,trim,s,report,o.distanceTolerance,0,o))
            {
                if (report.boundaryHausdorff!=CheckState::Exceeded)
                {
                    report.boundaryHausdorff=CheckState::Unresolved;report.pcurveLinearity=CheckState::Unresolved;
                    report.pcurveFailure=PCurveFailure::ValidationUnresolved;
                    report.message="Linear PCurve continuous boundary bound unresolved: boundary "+std::to_string(edge);
                }
                report.boundaryHausdorffUpperBound=INFINITY;
                return report;
            }
        }
        const auto &line=trim.edges[edge];
        I u(std::min(line.start.X(),line.end.X()),std::max(line.start.X(),line.end.X()));
        I v(std::min(line.start.Y(),line.end.Y()),std::max(line.start.Y(),line.end.Y()));
        I speed=Hypot(MotherTangent(parameters,u,v,I(line.end.X()-line.start.X()),I(line.end.Y()-line.start.Y())));
        double gap=std::max(startCoordinate.hi,1-endCoordinate.lo);
        I bound=I(report.boundaryHausdorffUpperBound)+speed*I(gap);
        if (!startCoordinate.finite() || !endCoordinate.finite() || !bound.finite() || bound.hi>=o.distanceTolerance)
        {
            report.boundaryHausdorff=CheckState::Unresolved;report.pcurveLinearity=CheckState::Unresolved;
            report.pcurveFailure=PCurveFailure::ValidationUnresolved;report.failedBoundary=static_cast<BoundaryId>(edge);
            report.message="PCurve endpoint coverage bound unresolved: boundary "+std::to_string(edge);
            report.boundaryHausdorffUpperBound=INFINITY;return report;
        }
        report.boundaryHausdorffUpperBound=std::max(previousBound,bound.hi);
    }
    report.boundaryHausdorff=CheckState::Passed;
    double bound = 0;
    bool exceeded = false, unresolved = false;
    std::function<void(const BezierCell &, unsigned)> visit = [&](const BezierCell &cell, unsigned depth) {
        if (exceeded || unresolved) return;
        if (++report.piecesVisited > o.maxValidationPieces)
        { unresolved = true; return; }
        auto certified = BoundMother(cell, parameters);
        double finiteBound=FiniteDistance(cell,trim,certified.distance);
        if (certified.distance.finite() && finiteBound < o.distanceTolerance)
        {
            bound = std::max(bound, finiteBound);
            return;
        }
        double u = cell.u0 + (cell.u1 - cell.u0) * .5;
        double v = cell.v0 + (cell.v1 - cell.v0) * .5;
        const auto &p = samples.At(u, v,sggk::DerivType::NONE,sggk::DerivType::NONE,false);
        double distance = Norm(p.point - ProjectToMother(parameters, p.point));
        report.hasWitness = true; report.witness = {u, v};
        // 采样只用于寻找超差见证；单元通过必须依赖向外舍入的连续误差上界。
        if (certified.distance.lo > o.distanceTolerance ||
            distance > o.distanceTolerance * (1 + 1e-8))
        { exceeded = true; return; }
        if (depth >= o.maxSubdivisionDepth)
        { unresolved = true; return; }
        auto children = Split(cell, depth % 2 == 0);
        visit(children.first, depth + 1); visit(children.second, depth + 1);
    };
    for (const auto &cell : cells) visit(cell, 0);
    if (cells.empty()) unresolved = true;
    report.geometry = exceeded ? CheckState::Exceeded : unresolved ? CheckState::Unresolved : CheckState::Passed;
    if (report.Passed())
    {
        report.distanceUpperBound = bound;
        report.hasWitness = false;
        // 每条边使用端点为 0 和 1 的连续截断参数坐标，其像覆盖整条线性 PCurve；同一位移上界因此同时控制两个方向的边界距离。
        report.message = "Continuous source-to-finite-trim distance and bidirectional boundary Hausdorff bounds; four linear PCurves";
    }
    else report.message = exceeded ? "General geometric constraint exceeded" :
        "General continuous bound or projection regularity unresolved within budget";
    return report;
}
GeneralValidationReport ValidateGeneral(const sggk::BSplineSurface &s,const AnalyticParameters &a,
    const GeneralConversionOptions &options,SurfaceSampleContext &samples,GeneralTrim *outputTrim)
{
    GeneralValidationReport report;
    const auto &o=options.fitting.validation;
    if(!(o.distanceTolerance>0 && std::isfinite(o.distanceTolerance)) ||
       o.maxSubdivisionDepth>30 || !o.maxValidationPieces ||
       options.boundarySamplesPerSpan<3 || options.boundarySamplesPerSpan>64)
        throw Failure(Status::InvalidInput,"Invalid general PCurve validation options");
    GeneralTrim trim;
    unsigned consumed=0;
    auto accept=[&](const GeneralTrim &candidate,GeneralValidationReport &attempt) {
        double sampled=attempt.maxSampledPCurveDeviation;
        attempt=ValidateTrim(s,candidate,options,samples,consumed);
        attempt.maxSampledPCurveDeviation=sampled;
        consumed=attempt.piecesVisited;
        return attempt.Passed();
    };
    if(RecoverLinearPCurves(samples,options,a,trim,report,accept) && outputTrim) *outputTrim=trim;
    report.piecesVisited=consumed;
    return report;
}
}
