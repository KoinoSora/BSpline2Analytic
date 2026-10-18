/**
 * @file GeneralRationalBoundaryCurves.cpp
 * @brief 将有理参数域直线复合为可导出的三维 NURBS 边界。
 * 按参数图接缝分段并保持 C0 端点一致；从 SDK 读回真实输出系数后认证导出误差。
 */
#include "GeneralRationalBoundaryInternal.h"

namespace bspline2analytic::detail::rational_boundary
{
namespace
{
BPoly C(I x){return BPoly(0,0,x);}
BPoly L(I a,I b){BPoly p(1,0);p.c={a,b};return p;}
double Mid(I x)
{
    if(!x.finite())throw Failure(Status::VerificationFailed,"Nonfinite rational coefficient");
    return x.lo+(x.hi-x.lo)*.5;
}
std::vector<BPoly> Basis(unsigned n,const BPoly &x,paired::Products &products)
{
    std::vector<BPoly> b{C(I(1))};auto one=C(I(1))-x;
    for(unsigned k=1;k<=n;++k)
    {
        std::vector<BPoly> next(k+1);
        for(unsigned i=0;i<=k;++i)
        {
            if(i<k)next[i]=products.Mul(b[i],one);
            if(i)next[i]=next[i]+products.Mul(b[i-1],x);
        }
        b=std::move(next);
    }
    return b;
}
BPoly Compose(const BPoly &p,const BPoly &u,const BPoly &v,paired::Products &products)
{
    auto bu=Basis(p.u,u,products),bv=Basis(p.v,v,products);BPoly out;
    for(unsigned i=0;i<=p.u;++i)for(unsigned j=0;j<=p.v;++j)
        out=out+products.Mul(bu[i],bv[j])*p.at(i,j);
    return out;
}
const paired::TargetCell &CellAt(const paired::Chart &chart,std::array<double,2> q)
{
    for(const auto &c:chart.cells)
        if(q[0]>=c.surface.u0&&q[0]<=c.surface.u1&&q[1]>=c.surface.v0&&q[1]<=c.surface.v1)return c;
    throw Failure(Status::VerificationFailed,"Point outside rational chart");
}
}
// 将参数直线 u(t)、v(t) 代入张量积齐次 Bernstein 形式，分别复合三维分子和权重。
// 必须在同一参数图单元内处理；若跨接缝，ComposeEdge 会先分割再调用，避免跨公式混算。
BezierCell ComposeLine(const paired::TargetCell &cell,sggk::Point2D a,sggk::Point2D b,
                       double lo,double hi,paired::Products &products)
{
    const auto &s=cell.surface;
    auto axis=[&](double start,double end,double lower,double upper)
    {
        I origin=(I(start)-I(lower))/(I(upper)-I(lower));
        if(start==end)return C(origin);
        I slope=(I(end)-I(start))/(I(upper)-I(lower));
        return L(origin+slope*I(lo),origin+slope*I(hi));
    };
    auto u=axis(a.X(),b.X(),s.u0,s.u1),v=axis(a.Y(),b.Y(),s.v0,s.v1);
    return {{Compose(s.p.x,u,v,products),Compose(s.p.y,u,v,products),Compose(s.p.z,u,v,products)},
            lo,hi,0,1,Compose(s.w,u,v,products)};
}
std::vector<BezierCell> ComposeEdge(const Trim &t,unsigned edge,paired::Products &products)
{
    auto a=t.polygon.corners[Ends[edge][0]],b=t.polygon.corners[Ends[edge][1]];
    std::vector<double> breaks{0,1};
    for(unsigned d=0;d<2;++d)
    {
        double x=d?a.Y():a.X(),y=d?b.Y():b.X();
        if(x==y)continue;
        unsigned n=unsigned(t.chart.angles[d].size()-1);
        for(unsigned i=1;i<n;++i)
        {
            double f=(double(i)/n-x)/(y-x);
            if(f>0&&f<1)breaks.push_back(f);
        }
    }
    std::sort(breaks.begin(),breaks.end());breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
    std::vector<BezierCell> out;
    for(size_t i=0;i+1<breaks.size();++i)
    {
        double f=breaks[i]+(breaks[i+1]-breaks[i])*.5;
        auto &cell=CellAt(t.chart,{a.X()+(b.X()-a.X())*f,a.Y()+(b.Y()-a.Y())*f});
        // 舍入后的分割参数精确定义输出分段。
        // 接缝附近公式可能延伸若干 ulp，但仍为母面上的有理轨迹；所有端点及后续 C0 拼接误差均计入输出误差界。
        out.push_back(ComposeLine(cell,a,b,breaks[i],breaks[i+1],products));
    }
    return out;
}
V ChartPoint(const Trim &t,std::array<double,2> uv,std::array<V,2> *derivatives)
{
    const auto &s=CellAt(t.chart,uv).surface;
    double u=(uv[0]-s.u0)/(s.u1-s.u0),v=(uv[1]-s.v0)/(s.v1-s.v0);
    // 这里只做采样和候选求值，不构成连续上界；输出系数及边界判据另行认证。
    auto value=[&](const BPoly &p,unsigned derivative=0)
    {
        if(p.u>2||p.v>2)throw Failure(Status::Unsupported,"Unexpected rational chart degree");
        if((derivative==1&&!p.u)||(derivative==2&&!p.v))return 0.;
        unsigned degreeV=p.v-(derivative==2?1:0);
        double rows[3]{};
        for(unsigned i=0;i<=p.u;++i)
        {
            double a[3]{};for(unsigned j=0;j<=degreeV;++j)
                a[j]=derivative==2?p.v*(Mid(p.at(i,j+1))-Mid(p.at(i,j))):Mid(p.at(i,j));
            for(unsigned n=degreeV;n>0;--n)for(unsigned j=0;j<n;++j)a[j]=a[j]*(1-v)+a[j+1]*v;
            rows[i]=a[0];
        }
        unsigned degreeU=p.u;
        if(derivative==1){for(unsigned i=0;i<p.u;++i)rows[i]=p.u*(rows[i+1]-rows[i]);--degreeU;}
        for(unsigned n=degreeU;n>0;--n)for(unsigned i=0;i<n;++i)rows[i]=rows[i]*(1-u)+rows[i+1]*u;
        return rows[0];
    };
    double w=value(s.w);
    if(!(w>0&&std::isfinite(w)))throw Failure(Status::VerificationFailed,"Invalid chart sample weight");
    V local{value(s.p.x)/w,value(s.p.y)/w,value(s.p.z)/w};
    V point=Get(t.mother.origin)+local;
    if(!Finite(point))throw Failure(Status::VerificationFailed,"Nonfinite chart sample point");
    if(derivatives)for(unsigned d=0;d<2;++d)
    {
        double span=d?s.v1-s.v0:s.u1-s.u0,dw=value(s.w,d+1);
        (*derivatives)[d]=(V{value(s.p.x,d+1),value(s.p.y,d+1),value(s.p.z,d+1)}-local*dw)*(1/(w*span));
    }
    return point;
}
double RefineFraction(const Trim &t,unsigned edge,V point,double initial)
{
    auto a=t.polygon.corners[Ends[edge][0]],b=t.polygon.corners[Ends[edge][1]];
    double dx=b.X()-a.X(),dy=b.Y()-a.Y(),f=std::clamp(initial,0.,1.);
    if(!std::isfinite(f))return initial;
    // 这里只提出恢复候选；最多六步、每步最多四次试探求值以限制开销。
    // 不调用 SDK 最近点搜索或无约束求根器。
    for(unsigned iteration=0;iteration<6;++iteration)
    {
        std::array<V,2> derivative;
        V residual=ChartPoint(t,{a.X()+dx*f,a.Y()+dy*f},&derivative)-point;
        V tangent=derivative[0]*dx+derivative[1]*dy;
        double speed2=Dot(tangent,tangent),error2=Dot(residual,residual);
        if(!(speed2>0&&std::isfinite(speed2)&&std::isfinite(error2)))break;
        double step=Dot(residual,tangent)/speed2;bool improved=false;
        for(unsigned backtrack=0;backtrack<4;++backtrack,step*=.5)
        {
            double next=std::clamp(f-step,0.,1.);if(next==f)break;
            V nextResidual=ChartPoint(t,{a.X()+dx*next,a.Y()+dy*next})-point;
            if(Dot(nextResidual,nextResidual)<error2){f=next;improved=true;break;}
        }
        if(!improved)break;
    }
    return f;
}
// 构造理想有理轨迹后导出 SDK 曲线，再读回其真实控制点、权重和结点进行认证。
// 端点统一和 C0 拼接可能引入舍入扰动，必须纳入 outputTolerance，而不能仅验证导出前系数。
Target BuildTarget(const Trim &t,double outputTolerance,GeneralRationalBoundaryReport &report,
                   const ExecutionDeadline *deadline)
{
    Target result;auto &patch=result.patch;patch.mother=t.mother;patch.chartDomain=t.domain;
    patch.rationalCorners=t.polygon.corners;patch.tangentNotApplicable=t.collapsed;
    for(unsigned i=0;i<4;++i)
    {
        auto p=t.polygon.corners[i];patch.vertices[i]=Point(ChartPoint(t,{p.X(),p.Y()}));
    }
    auto triangles=TrimTriangles(t.polygon);
    if(triangles.empty())throw Failure(Status::VerificationFailed,"Empty rational trim");
    auto &tri=triangles.front();
    patch.interiorReference=Point(ChartPoint(t,{(tri[0].X()+tri[1].X()+tri[2].X())/3,
                                                (tri[0].Y()+tri[1].Y()+tri[2].Y())/3}));
    paired::Products products;
    for(unsigned e=0;e<4;++e)
    {
        if(deadline)deadline->Check("target_preparation");
        auto ideal=ComposeEdge(t,e,products);unsigned degree=1;
        for(const auto &c:ideal)degree=std::max({degree,c.p.x.u,c.p.y.u,c.p.z.u,c.w.u});
        report.maxTargetDegree=std::max(report.maxTargetDegree,degree);
        sggk::Point3DArray cp;sggk::RealArray weights,knots{0};sggk::UIntArray mults{degree+1};
        for(size_t k=0;k<ideal.size();++k)
        {
            if(deadline)deadline->Check("target_preparation");
            auto &c=ideal[k];c.p={products.Elevate(c.p.x,degree,0),products.Elevate(c.p.y,degree,0),
                                  products.Elevate(c.p.z,degree,0)};c.w=products.Elevate(c.w,degree,0);
            if(c.w.range().lo<=0)throw Failure(Status::VerificationFailed,"Nonpositive target weight");
            std::vector<double> segmentWeights(degree+1);sggk::Point3DArray segmentPoints(degree+1);
            double scale=k?weights.back()/Mid(c.w.c[0]):1;
            for(unsigned i=0;i<=degree;++i)
            {
                segmentWeights[i]=Mid(c.w.c[i])*scale;
                segmentPoints[i]=Point(Get(t.mother.origin)+V{Mid(c.p.x.c[i]/c.w.c[i]),
                    Mid(c.p.y.c[i]/c.w.c[i]),Mid(c.p.z.c[i]/c.w.c[i])});
            }
            segmentPoints[0]=k?cp.back():patch.vertices[Ends[e][0]];
            if(k)segmentWeights[0]=weights.back();
            if(k+1==ideal.size())segmentPoints.back()=patch.vertices[Ends[e][1]];
            for(unsigned i=0;i<=degree;++i)
            {
                double w=segmentWeights[i];
                if(!(w>0&&std::isfinite(w)))throw Failure(Status::VerificationFailed,"Invalid exported weight");
                if(!k||i){cp.push_back(segmentPoints[i]);weights.push_back(w);}
            }
            knots.push_back(c.u1);mults.push_back(k+1==ideal.size()?degree+1:degree);
        }
        patch.boundaries[e]=std::make_shared<sggk::BSplineCurve3D>(degree,cp,weights,knots,mults);
        // 认证前必须从 SDK 对象读回实际系数；尤其要处理 SDK 将等权有理曲线规范化为多项式曲线的情况。
        const auto &curve=*patch.boundaries[e];
        if(curve.Degree()!=degree||curve.ControlPoints().size()!=cp.size()||curve.Knots()!=knots||curve.Mults()!=mults)
            throw Failure(Status::VerificationFailed,"SDK changed exported rational curve");
        for(size_t i=0;i<cp.size();++i)
        {
            auto a=Get(curve.ControlPoints()[i]),b=Get(cp[i]);
            if(a.x!=b.x||a.y!=b.y||a.z!=b.z)throw Failure(Status::VerificationFailed,"SDK changed curve controls");
        }
        for(size_t k=0;k<ideal.size();++k)
        {
            if(deadline)deadline->Check("target_export_certificate");
            const auto &c=ideal[k];
            BezierCell actual{{BPoly(degree,0),BPoly(degree,0),BPoly(degree,0)},c.u0,c.u1,0,1,BPoly(degree,0)};
            for(unsigned i=0;i<=degree;++i)
            {
                size_t index=k*degree+i;double w=curve.IsRational()?curve.Weights()[index]:1;
                if(!(w>0&&std::isfinite(w)))throw Failure(Status::VerificationFailed,"Invalid SDK weight");
                IV p=(Exact(Get(curve.ControlPoints()[index]))-Exact(Get(t.mother.origin)))*I(w);
                actual.p.x.c[i]=p.x;actual.p.y.c[i]=p.y;actual.p.z.c[i]=p.z;actual.w.c[i]=I(w);
            }
            if(!curve.IsRational())actual.w=C(I(1));
            auto g=products.Subtract(products.MulShared(actual.p,c.w),products.MulShared(c.p,actual.w));
            auto h=products.Mul(actual.w,c.w);
            double bound=paired::PositionBoxUpper(g,h);
            if(!(bound<outputTolerance))bound=paired::PositionControlUpper(g,h,products);
            if(!(bound<outputTolerance))throw Failure(Status::VerificationFailed,"Export-to-mother bound unresolved");
            result.motherBound=std::max(result.motherBound,bound);
            result.edges[e].push_back(std::move(actual));
        }
        report.targetSegments+=unsigned(ideal.size());
    }
    return result;
}
}
