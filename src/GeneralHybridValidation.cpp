/**
 * @file GeneralHybridValidation.cpp
 * @brief 已废弃的解析 UV 混合验证器与共享一维提取实现。
 * 完整验证接口禁测；ExtractHybridBoundaries 等提取辅助函数仍可被当前有理边界验证器复用。
 */
#include "GeneralHybridInternal.h"
#include "PairedValidationInternal.h"
#include <chrono>

namespace bspline2analytic::detail
{
namespace
{
using Clock = std::chrono::steady_clock;
double Milliseconds(Clock::time_point start)
{ return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }
struct Measure
{
    double &value; Clock::time_point start=Clock::now();
    explicit Measure(double &v):value(v){}
    ~Measure(){value+=Milliseconds(start);}
};
struct Budget : std::runtime_error { using std::runtime_error::runtime_error; };
std::vector<double> Expand(const sggk::RealArray &k,const sggk::UIntArray &m)
{
    std::vector<double> out;
    for(size_t i=0;i<k.size();++i)out.insert(out.end(),m[i],k[i]);
    return out;
}
// 在次数个参数处计算多仿射形式（blossom）：参数全为 t 时求值，重复端点时提取 Bézier 系数；所有运算向外舍入。
std::vector<I> Coefficients(const std::vector<double> &k,unsigned p,unsigned span,
                            double lo,double hi,unsigned left)
{
    std::vector<std::vector<I>> d(p+1,std::vector<I>(p+1));
    for(unsigned i=0;i<=p;++i)d[i][i]=I(1);
    for(unsigned r=1;r<=p;++r)
    {
        double t=r<=left?lo:hi;
        for(unsigned i=p;i>=r;--i)
        {
            I alpha=(I(t)-I(k[span-p+i]))/(I(k[span+1+i-r])-I(k[span-p+i]));
            for(unsigned j=0;j<=p;++j)d[i][j]=d[i-1][j]*(I(1)-alpha)+d[i][j]*alpha;
        }
    }
    return d[p];
}
bool Zero(I x){return x.lo==0&&x.hi==0;}
BPoly Constant(I x){return BPoly(0,0,x);}
BPoly Linear(I a,I b){BPoly p(1,0);p.c={a,b};return p;}
struct HybridThresholds
{
    double distanceTolerance,tangentRatioLimit;
    HybridThresholds(double distance,double degrees):distanceTolerance(distance)
    {
        I angle=I(degrees)*IPi()/I(180);
        tangentRatioLimit=(Sin(angle)/Cos(angle)).lo;
    }
};
// 恢复得到的参数图独占坐标架和 PCurve，不能跨尝试共享该上下文，因为球面恢复可能改变坐标架和边方向。
struct TargetBoundaryContext
{
    struct Edge {I u,v,du,dv;};
    sggk::SurfaceType type;
    paired::LocalFrame frame;
    I radius,majorRadius,halfRadius,coneSlope;
    std::array<Edge,4> edges;
    explicit TargetBoundaryContext(const GeneralTrim &trim)
        :type(trim.parameters.type),frame(paired::MakeLocalFrame(trim.parameters)),
         radius(trim.parameters.radius),majorRadius(trim.parameters.majorRadius),
         halfRadius(radius/I(2))
    {
        if(type==sggk::SurfaceType::Cone)
            coneSlope=Sin(I(trim.parameters.halfAngle))/Cos(I(trim.parameters.halfAngle));
        for(unsigned e=0;e<4;++e)
        {
            const auto &line=trim.edges[e];
            edges[e]={I(line.start.X()),I(line.start.Y()),
                I(line.end.X())-I(line.start.X()),I(line.end.Y())-I(line.start.Y())};
        }
    }
};
// 采用含一致十阶余项界的十项 Taylor 包络；后续切向使用独立包络的解析导数，不能对该区间位置近似直接求导代替。
struct TrigPair {BPoly sin,cos;};
TrigPair EncloseSinCos(const BPoly &angle,paired::Products &products)
{
    I range=angle.range();
    if(!range.finite())return {Constant(I::Whole()),Constant(I::Whole())};
    double center=range.lo+(range.hi-range.lo)*.5;
    BPoly d=angle-Constant(I(center)),power=Constant(I(1));TrigPair out;
    I sn=Sin(I(center)),cs=Cos(I(center)),factorial(1);
    const std::array<I,4> derivatives{sn,cs,-sn,-cs};
    for(unsigned k=0;k<10;++k)
    {
        // 复用角度幂，但分别保留 sin/cos 的误差包络。
        out.sin=out.sin+power*(derivatives[k%4]/factorial);
        out.cos=out.cos+power*(derivatives[(k+1)%4]/factorial);
        if(k<9)power=products.Mul(power,d);
        factorial=factorial*I(k+1);
    }
    I remainder(1),radius(d.range().absMax());
    for(unsigned k=1;k<=10;++k)remainder=remainder*radius/I(k);
    auto tail=Constant(I(-remainder.hi,remainder.hi));
    return {out.sin+tail,out.cos+tail};
}
struct Target {BV point,tangent;};
Target TargetBoundary(const TargetBoundaryContext &a,unsigned edge,double lo,double hi,
                       paired::Products &products)
{
    const auto &line=a.edges[edge];I du=line.du,dv=line.dv;
    BPoly u=Linear(line.u+du*I(lo),line.u+du*I(hi));
    BPoly v=Linear(line.v+dv*I(lo),line.v+dv*I(hi));
    auto trig=[&](const BPoly &p){return EncloseSinCos(p,products);};
    BPoly zero=Constant(I(0));BV p,t;
    using T=sggk::SurfaceType;
    if(a.type==T::Plane){p={u,v,zero};t={Constant(du),Constant(dv),zero};}
    else if(a.type==T::Cylinder)
    {auto q=trig(u);p={q.cos*a.radius,q.sin*a.radius,v};t={q.sin*(-du*a.radius),q.cos*(du*a.radius),Constant(dv)};}
    else if(a.type==T::Cone)
    {
        auto q=trig(u);I h=a.coneSlope;
        auto vc=products.Mul(v,q.cos),vs=products.Mul(v,q.sin);
        p={vc*h,vs*h,v};
        t={(q.cos*dv-vs*du)*h,(q.sin*dv+vc*du)*h,Constant(dv)};
    }
    else if(a.type==T::Sphere)
    {
        auto plus=trig(v+u),minus=trig(v-u),latitude=trig(v);I r=a.halfRadius;
        p={(plus.sin+minus.sin)*r,(minus.cos-plus.cos)*r,latitude.cos*(-a.radius)};
        t={(plus.cos*(dv+du)+minus.cos*(dv-du))*r,
           (plus.sin*(dv+du)-minus.sin*(dv-du))*r,latitude.sin*(a.radius*dv)};
    }
    else
    {
        auto plus=trig(u+v),minus=trig(u-v),longitude=trig(u),latitude=trig(v);I r=a.halfRadius;
        p={longitude.cos*a.majorRadius+(plus.cos+minus.cos)*r,
           longitude.sin*a.majorRadius+(plus.sin+minus.sin)*r,latitude.sin*a.radius};
        t={longitude.sin*(-a.majorRadius*du)-(plus.sin*(du+dv)+minus.sin*(du-dv))*r,
           longitude.cos*(a.majorRadius*du)+(plus.cos*(du+dv)+minus.cos*(du-dv))*r,latitude.cos*(a.radius*dv)};
    }
    const auto &f=a.frame;
    auto rotate=[&](const BV &q){return BV{q.x*f.x.x+q.y*f.y.x+q.z*f.z.x,
        q.x*f.x.y+q.y*f.y.y+q.z*f.z.y,q.x*f.x.z+q.y*f.y.z+q.z*f.z.z};};
    return {rotate(p),rotate(t)}; // 源几何也以同一个原点中心化。
}
double LineFraction(const GeneralTrim &trim,unsigned edge,V world,double reference)
{
    const auto &a=trim.parameters;const auto &line=trim.edges[edge];
    V local=FrameOf(a).local(ProjectToMother(a,world));
    auto uv=MotherParameters(a,local);auto ref=EdgeUV(line,reference);
    for(unsigned d=0;d<2;++d)
        if((d==0&&a.type!=sggk::SurfaceType::Plane)||(d==1&&a.type==sggk::SurfaceType::Torus))
            uv[d]+=2*Pi*std::round((ref[d]-uv[d])/(2*Pi));
    auto scale=TrimScale(trim);
    double dx=(line.end.X()-line.start.X())*scale[0],dy=(line.end.Y()-line.start.Y())*scale[1];
    double length=dx*dx+dy*dy;
    if(!(length>0))return reference;
    return ((uv[0]-line.start.X())*scale[0]*dx+(uv[1]-line.start.Y())*scale[1]*dy)/length;
}
// 此处构造的有限面片点只提供距离上界见证，不声称加权 UV 截断能解出精确最近点。
double FiniteSampleDistance(const GeneralTrim &trim,V p)
{
    const auto &a=trim.parameters;auto frame=FrameOf(a);
    V projected=ProjectToMother(a,p);auto uv=MotherParameters(a,frame.local(projected));
    std::array<double,2> center{};
    for(const auto &c:trim.corners){center[0]+=c.X()*.25;center[1]+=c.Y()*.25;}
    for(unsigned d=0;d<2;++d)
        if((d==0&&a.type!=sggk::SurfaceType::Plane)||(d==1&&a.type==sggk::SurfaceType::Torus))
            uv[d]+=2*Pi*std::round((center[d]-uv[d])/(2*Pi));
    if(InsideTrim(trim,uv,0))return Norm(p-projected);
    // 内部三角剖分对角线上的点可能同时未通过两个严格浮点符号测试。
    // 应在两个闭三角形内分别构造实际点，避免错误地改求到外边界的距离。
    double best=INFINITY;
    for(const auto &triangle:TrimTriangles(trim))
    {
        double x=triangle[1].X()-triangle[0].X(),y=triangle[1].Y()-triangle[0].Y();
        double z=triangle[2].X()-triangle[0].X(),w=triangle[2].Y()-triangle[0].Y();
        double determinant=x*w-y*z;
        if(determinant==0)continue;
        double px=uv[0]-triangle[0].X(),py=uv[1]-triangle[0].Y();
        double b=(px*w-py*z)/determinant,c=(x*py-y*px)/determinant;
        double weights[3]={std::max(0.,1-b-c),std::max(0.,b),std::max(0.,c)};
        double sum=weights[0]+weights[1]+weights[2];
        std::array<double,2> on{};
        for(unsigned i=0;i<3;++i){on[0]+=triangle[i].X()*(weights[i]/sum);on[1]+=triangle[i].Y()*(weights[i]/sum);}
        best=std::min(best,Norm(p-frame.global(CanonicalPoint(a,on[0],on[1]))));
    }
    return best;
}
bool SourceCollapsed(const sggk::BSplineSurface &s,unsigned edge)
{
    bool u=edge<2,high=edge%2!=0;const auto &cp=s.ControlPoints();
    const auto &knots=u?s.KnotsU():s.KnotsV();const auto &mults=u?s.MultsU():s.MultsV();
    unsigned degree=u?s.DegreeU():s.DegreeV();double value=u?(high?s.MaxParamU():s.MinParamU()):(high?s.MaxParamV():s.MinParamV());
    if(value!=(high?knots.back():knots.front())||(high?mults.back():mults.front())!=degree+1)return false;
    size_t fixed=high?(u?cp.size():cp.front().size())-1:0;
    V first=Get(u?cp[fixed][0]:cp[0][fixed]);
    for(size_t i=0;i<(u?cp.front().size():cp.size());++i)
    {V p=Get(u?cp[fixed][i]:cp[i][fixed]);if(p.x!=first.x||p.y!=first.y||p.z!=first.z)return false;}
    return true;
}
class BoundaryValidator
{
    const sggk::BSplineSurface &source;const GeneralTrim &trim;
    const GeneralHybridOptions &options;SurfaceSampleContext &samples;GeneralHybridReport &out;
    const HybridThresholds &thresholds;
    const TargetBoundaryContext targetContext;
    paired::Products products;
    double maxDistance=0,maxTangentRatio=0;
    bool Visit(const BezierCell &cell,unsigned edge,double a,double b,unsigned depth,bool collapsed)
    {
        if(++out.boundaryNodes>options.maxBoundaryPieces)throw Budget("Boundary node budget exhausted");
        auto target=TargetBoundary(targetContext,edge,a,b,products);
        auto error=products.Subtract(cell.p,products.MulShared(target.point,cell.w));
        double distance=paired::PositionBoxUpper(error,cell.w);
        if(distance>=thresholds.distanceTolerance)
            distance=std::min(distance,paired::PositionControlUpper(error,cell.w,products));
        double ratio=INFINITY;
        if(collapsed)ratio=0;
        else
        {
            BV tangent=Derivative(cell.p,true);
            if(cell.w.u||cell.w.v)
                tangent=products.Subtract(products.MulShared(tangent,cell.w),
                                          products.MulShared(cell.p,Derivative(cell.w,true)));
            I dot=products.Dot(tangent,target.tangent).range();
            IV cross=Range(products.Cross(tangent,target.tangent));
            if(dot.lo>0)ratio=(Sqrt(Sqr(cross.x)+Sqr(cross.y)+Sqr(cross.z))/I(dot.lo)).hi;
        }
        if(distance<thresholds.distanceTolerance&&ratio<thresholds.tangentRatioLimit)
        {
            maxDistance=std::max(maxDistance,distance);
            maxTangentRatio=std::max(maxTangentRatio,ratio);
            ++out.boundarySegments;return true;
        }
        if(depth>=options.maxBoundaryDepth)return false;
        double mid=cell.u0+(cell.u1-cell.u0)*.5;
        if(!(mid>cell.u0&&mid<cell.u1))return false;
        double mapped=a+(b-a)*.5;
        if(!collapsed)
        {
            bool u=edge<2;
            V p=samples.At(u?(edge==0?source.MinParamU():source.MaxParamU()):mid,
                u?mid:(edge==2?source.MinParamV():source.MaxParamV()),sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
            double proposed=LineFraction(trim,edge,p,mapped);
            if(std::isfinite(proposed)&&proposed>a&&proposed<b)mapped=proposed;
        }
        if(!(mapped>a&&mapped<b))return false;
        ++out.mappingRefinements;
        auto halves=Split(cell,true);
        return Visit(halves.first,edge,a,mapped,depth+1,collapsed)&&Visit(halves.second,edge,mapped,b,depth+1,collapsed);
    }
public:
    BoundaryValidator(const sggk::BSplineSurface &s,const GeneralTrim &t,const GeneralHybridOptions &o,
        SurfaceSampleContext &c,GeneralHybridReport &r,const HybridThresholds &limits)
        :source(s),trim(t),options(o),samples(c),out(r),thresholds(limits),targetContext(t){}
    bool Run(const std::array<std::vector<BezierCell>,4> &edges)
    {
        unsigned collapsedCount=0;
        for(unsigned e=0;e<4;++e)
        {
            out.failedBoundary=e;
            if(edges[e].empty())return false;
            bool collapsed=SourceCollapsed(source,e);
            out.tangentNotApplicable[e]=collapsed;
            if(collapsed)++collapsedCount;
            double lo=edges[e].front().u0,hi=edges[e].back().u1,previous=0;
            for(size_t i=0;i<edges[e].size();++i)
            {
                const auto &cell=edges[e][i];double end=i+1==edges[e].size()?1:(cell.u1-lo)/(hi-lo);
                if(i+1<edges[e].size()&&!collapsed)
                {
                    bool u=e<2;V p=samples.At(u?(e==0?source.MinParamU():source.MaxParamU()):cell.u1,
                        u?cell.u1:(e==2?source.MinParamV():source.MaxParamV()),sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
                    end=LineFraction(trim,e,p,end);
                }
                if(!(end>previous&&end<=1)||!Visit(cell,e,previous,end,0,collapsed))return false;
                previous=end;
            }
        }
        out.collapsedEdges=collapsedCount;out.boundaryDistanceUpperBound=maxDistance;
        // atan 单调，验收不需要换算为度；仅组装报告时转换角度单位。
        out.boundaryTangentUpperBoundDeg=(Atan(I(maxTangentRatio))*I(180)/IPi()).hi;
        out.boundaryPosition=out.boundaryTangent=CheckState::Passed;return true;
    }
};
}
std::array<std::vector<BezierCell>,4> ExtractHybridBoundaries(const sggk::BSplineSurface &s,V origin)
{
    std::array<std::vector<BezierCell>,4> result;
    const auto &cp=s.ControlPoints();
    auto ku=Expand(s.KnotsU(),s.MultsU()),kv=Expand(s.KnotsV(),s.MultsV());
    for(unsigned edge=0;edge<4;++edge)
    {
        bool fixedU=edge<2,high=edge%2!=0;
        const auto &fk=fixedU?ku:kv;const auto &vk=fixedU?kv:ku;
        unsigned fp=fixedU?s.DegreeU():s.DegreeV(),vp=fixedU?s.DegreeV():s.DegreeU();
        unsigned fn=unsigned(fixedU?cp.size():cp.front().size()),vn=unsigned(fixedU?cp.front().size():cp.size());
        double value=fixedU?(high?s.MaxParamU():s.MinParamU()):(high?s.MaxParamV():s.MinParamV());
        double lo=fixedU?s.MinParamV():s.MinParamU(),hi=fixedU?s.MaxParamV():s.MaxParamU();
        unsigned span=fp;while(span+1<fn&&fk[span+1]<=value)++span;
        auto fixed=Coefficients(fk,fp,span,value,value,fp);
        std::vector<IV> p(vn);std::vector<I>w(vn);
        for(unsigned j=0;j<vn;++j)for(unsigned i=0;i<=fp;++i)
        {
            if(Zero(fixed[i]))continue;
            unsigned u=fixedU?span-fp+i:j,v=fixedU?j:span-fp+i;
            I weight=fixed[i]*(s.IsRational()?I(s.Weights()[u][v]):I(1));
            p[j]=p[j]+(Exact(Get(cp[u][v]))-Exact(origin))*weight;w[j]=w[j]+weight;
        }
        for(unsigned k=vp;k<vn;++k)
        {
            if(vk[k]>=vk[k+1]||vk[k]>=hi||vk[k+1]<=lo)continue;
            BezierCell cell{{BPoly(vp,0),BPoly(vp,0),BPoly(vp,0)},vk[k],vk[k+1],0,1,BPoly(vp,0)};
            for(unsigned j=0;j<=vp;++j)
            {
                auto coefficients=Coefficients(vk,vp,k,vk[k],vk[k+1],vp-j);
                IV sum{};I weight;
                for(unsigned i=0;i<=vp;++i)if(!Zero(coefficients[i]))
                {sum=sum+p[k-vp+i]*coefficients[i];weight=weight+w[k-vp+i]*coefficients[i];}
                cell.p.x.c[j]=sum.x;cell.p.y.c[j]=sum.y;cell.p.z.c[j]=sum.z;cell.w.c[j]=weight;
            }
            // 多项式源面的基函数精确满足单位分解；保留 v 次的常数 1 表示会无谓地使乘积次数翻倍。
            if(!s.IsRational())cell.w=Constant(I(1));
            if(cell.u0<lo||cell.u1>hi)cell=paired::Crop(cell,std::max(lo,cell.u0),std::min(hi,cell.u1),0,1);
            result[edge].push_back(std::move(cell));
        }
    }
    return result;
}
GeneralHybridReport ValidateGeneralHybridWithSamples(const sggk::BSplineSurface &s,const AnalyticParameters &a,
    const GeneralHybridOptions &options,SurfaceSampleContext &samples)
{
    auto start=Clock::now();GeneralHybridReport out;
    try
    {
        ValidateSurface(s);
        const auto &v=options.general.fitting.validation;
        if(!(v.distanceTolerance>0&&std::isfinite(v.distanceTolerance))||
            !(v.boundaryTangentToleranceDeg>0&&v.boundaryTangentToleranceDeg<90)||
            !options.maxBoundaryPieces||options.maxBoundaryDepth>20||
            options.general.boundarySamplesPerSpan<3||options.general.boundarySamplesPerSpan>64||
            !options.general.fitting.sampling.maxSamples)
            throw Failure(Status::InvalidInput,"Invalid hybrid validation options");
        if(!IntervalEnvironmentSupported())throw Failure(Status::Unsupported,"Directed interval environment unavailable");
        if(s.IsRational())for(const auto &row:s.Weights())for(double w:row)
            if(!(w>0))throw Failure(Status::Unsupported,"Positive weights required");
        for(auto axis:{false,true})
        {
            const auto &m=axis?s.MultsU():s.MultsV();unsigned p=axis?s.DegreeU():s.DegreeV();
            for(size_t i=1;i+1<m.size();++i)if(m[i]>p)throw Failure(Status::Unsupported,"Discontinuous source knot");
        }
        auto frame=FrameOf(a);
        if(!Finite(frame.o)||!Finite(frame.x)||!Finite(frame.z))throw Failure(Status::InvalidInput,"Invalid mother frame");
        {
            Measure timer(out.samplingMs);
            const auto &points=samples.Greville(false);
            if(points.empty())throw Failure(Status::Unsupported,"No interior samples");
            for(const auto &p:points)
            {
                double d=Norm(p.point-ProjectToMother(a,p.point));++out.interiorSampleCount;
                out.interiorSampleMaxDistance=std::max(out.interiorSampleMaxDistance,d);
            }
        }
        if(!(out.interiorSampleMaxDistance<v.distanceTolerance))
        {out.interiorSamples=CheckState::Exceeded;out.message="Interior sample distance exceeded";}
        else
        {
            const auto thresholds=[&]
            {
                Measure timer(out.boundaryMs);
                return HybridThresholds(v.distanceTolerance,v.boundaryTangentToleranceDeg);
            }();
            std::optional<std::array<std::vector<BezierCell>,4>> edges;
            GeneralTrim trim;GeneralValidationReport recovery;
            auto recoveryStart=Clock::now();double before=out.samplingMs+out.preparationMs+out.boundaryMs;
            auto accept=[&](const GeneralTrim &candidate,GeneralValidationReport &)
            {
                {
                    Measure timer(out.samplingMs);double maximum=0;
                    for(const auto &p:samples.Greville(false))maximum=std::max(maximum,FiniteSampleDistance(candidate,p.point));
                    out.interiorSampleMaxFiniteDistance=maximum;
                    if(!(maximum<v.distanceTolerance)){out.message="Interior finite-patch sample witness exceeds tolerance";return false;}
                }
                if(!edges){Measure timer(out.preparationMs);edges=ExtractHybridBoundaries(s,Get(a.origin));}
                Measure timer(out.boundaryMs);
                BoundaryValidator validator(s,candidate,options,samples,out,thresholds);
                bool passed=validator.Run(*edges);
                if(!passed)out.message="Boundary correspondence unresolved on edge "+std::to_string(out.failedBoundary);
                if(passed)out.interiorSamples=CheckState::Passed;
                return passed;
            };
            bool passed=RecoverLinearPCurves(samples,options.general,a,trim,recovery,accept);
            out.recoveryMs=Milliseconds(recoveryStart)-(out.samplingMs+out.preparationMs+out.boundaryMs-before);
            out.chartAttempts=recovery.chartAttempts;out.chartIterations=recovery.chartIterations;
            if(passed)
            {
                out.patch.emplace(trim.parameters,std::array<double,4>{s.MinParamU(),s.MaxParamU(),s.MinParamV(),s.MaxParamV()},trim.corners,trim.edges);
                out.status=HybridValidationStatus::Passed;
                out.message="Sampled interior (augmented Greville grid); continuous boundary position and tangent; collapsed-edge tangent not applicable";
            }
            else if(out.message.empty())out.message="Hybrid validation not certified: "+recovery.message;
            if(out.boundaryNodes>options.maxBoundaryPieces)out.status=HybridValidationStatus::BudgetExceeded;
        }
    }
    catch(const Budget &e){out.status=HybridValidationStatus::BudgetExceeded;out.message=e.what();}
    catch(const Failure &e)
    {out.status=e.status==Status::InvalidInput?HybridValidationStatus::InvalidInput:e.status==Status::Unsupported?HybridValidationStatus::Unsupported:HybridValidationStatus::NotCertified;out.message=e.what();}
    catch(const std::exception &e){out.status=HybridValidationStatus::NotCertified;out.message=e.what();}
    if(out.status!=HybridValidationStatus::Passed)
    {
        out.patch.reset();out.boundaryDistanceUpperBound=out.boundaryTangentUpperBoundDeg=INFINITY;
        out.tangentNotApplicable={};
        if(out.interiorSamples!=CheckState::Exceeded)out.interiorSamples=CheckState::Unresolved;
        out.boundaryPosition=out.boundaryTangent=CheckState::Unresolved;
    }
    out.totalMs=Milliseconds(start);return out;
}
}
namespace bspline2analytic
{
GeneralHybridReport ValidateGeneralHybrid(const sggk::BSplineSurface &s,const AnalyticParameters &a,const GeneralHybridOptions &o)
{
    auto start=std::chrono::steady_clock::now();detail::SurfaceSampleContext samples(s,o.general.fitting.sampling.maxSamples);
    auto out=detail::ValidateGeneralHybridWithSamples(s,a,o,samples);
    out.totalMs=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();return out;
}
const char *HybridValidationStatusName(HybridValidationStatus s) noexcept
{
    switch(s){case HybridValidationStatus::Passed:return "Passed";case HybridValidationStatus::NotCertified:return "NotCertified";
    case HybridValidationStatus::InvalidInput:return "InvalidInput";case HybridValidationStatus::Unsupported:return "Unsupported";
    default:return "BudgetExceeded";}
}
}
