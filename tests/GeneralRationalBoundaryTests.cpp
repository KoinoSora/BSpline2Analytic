/**
 * @file GeneralRationalBoundaryTests.cpp
 * @brief 当前有理边界验证器的单元、数据与性能驱动。
 * 独立审计导出系数和边界；区分给定母面的验证模式与包含分类的算法模式，旧混合验证器禁用。
 */
#include "GeneralRationalBoundaryInternal.h"
#include "GeneralDataset.h"
#include "GeneralExtremeDataset.h"
#include "GeneralRationalDataset.h"
#include <Foundation/init.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <sstream>

using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace rb=bspline2analytic::detail::rational_boundary;
namespace
{
unsigned assertions=0;
void Require(bool p,const std::string &m){++assertions;if(!p)throw std::runtime_error(m);}
double Ms(std::chrono::steady_clock::time_point t){return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count();}
struct StageTimer
{
    double &out;std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
    explicit StageTimer(double &value):out(value){}~StageTimer(){out+=Ms(start);}
};
double Mid(I a){return a.lo+(a.hi-a.lo)*.5;}
std::string Quote(std::string s){std::string o="\"";for(char c:s){if(c=='"')o+='"';o+=c=='\n'?' ':c;}return o+'"';}
AnalyticParameters Mother(unsigned type)
{
    AnalyticParameters a;a.type=GeneralDatasetType(type);a.radius=2;a.majorRadius=5;a.halfAngle=.4;
    a.origin={1,-2,3};a.axis=Direction({.2,.3,1});a.xDirection=Direction(Perpendicular(Get(a.axis)));return a;
}
std::vector<BPoly> Basis(unsigned n,const BPoly &p)
{
    BPoly one(0,0,I(1));std::vector<BPoly> b{one};
    for(unsigned k=1;k<=n;++k)
    {std::vector<BPoly> q(k+1);for(unsigned i=0;i<=k;++i){if(i<k)q[i]=b[i]*(one-p);if(i)q[i]=q[i]+b[i-1]*p;}b=std::move(q);}
    return b;
}
BPoly ComposeSurface(const BPoly &p,const BPoly &u,const BPoly &v)
{
    auto a=Basis(p.u,u),b=Basis(p.v,v);BPoly q;
    for(unsigned i=0;i<=p.u;++i)for(unsigned j=0;j<=p.v;++j)q=q+(a[i]*b[j])*p.at(i,j);
    return q;
}
sggk::BSplineSurfacePtr Fixture(const AnalyticParameters &a,const std::array<double,4> &domain,bool skew,bool heavy=false)
{
    AnalyticPatch patch;patch.parameters=a;patch.domain=domain;
    auto chart=paired::BuildChart(patch);Require(chart.cells.size()==1,"Fixture expects one chart cell");
    paired::PlaceChart(chart,a,{0,0,0});const auto &c=chart.cells.front().surface;
    BPoly u(1,1),v(1,1);
    u.c=skew?std::vector<I>{I(.05),I(.2),I(.85),I(.95)}:std::vector<I>{I(0),I(0),I(1),I(1)};
    v.c=skew?std::vector<I>{I(.1),I(.8),I(.25),I(.95)}:std::vector<I>{I(0),I(1),I(0),I(1)};
    auto x=ComposeSurface(c.p.x,u,v),y=ComposeSurface(c.p.y,u,v),z=ComposeSurface(c.p.z,u,v),w=ComposeSurface(c.w,u,v);
    unsigned p=std::max({x.u,y.u,z.u,w.u}),q=std::max({x.v,y.v,z.v,w.v});paired::Products products;
    x=products.Elevate(x,p,q);y=products.Elevate(y,p,q);z=products.Elevate(z,p,q);w=products.Elevate(w,p,q);
    sggk::Point3DMatrix cp(p+1,std::vector<sggk::Point3D>(q+1));sggk::RealMatrix weights(p+1,std::vector<double>(q+1));
    for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
    {double d=Mid(w.at(i,j));weights[i][j]=d;cp[i][j]={Mid(x.at(i,j))/d,Mid(y.at(i,j))/d,Mid(z.at(i,j))/d};}
    auto s=std::make_shared<sggk::BSplineSurface>(p,q,cp,weights,sggk::RealArray{0,1},sggk::RealArray{0,1},
        sggk::UIntArray{p+1,p+1},sggk::UIntArray{q+1,q+1});
    if(heavy)
    {
        if(p<5)s->UDegreeElevation(5-p);if(q<5)s->VDegreeElevation(5-q);
        for(unsigned i=1;s->ControlPoints().size()<10;++i)s->InsertUKnots(std::fmod(i*.6180339887498949,1.),1);
        for(unsigned i=1;s->ControlPoints().front().size()<20;++i)s->InsertVKnots(std::fmod(i*.6180339887498949,1.),1);
    }
    return s;
}
double MotherDistance(const AnalyticParameters &a,V point)
{
    V p=FrameOf(a).local(point);double rho=std::hypot(p.x,p.y);
    using T=sggk::SurfaceType;
    switch(a.type)
    {
    case T::Plane:return std::abs(p.z);
    case T::Cylinder:return std::abs(rho-a.radius);
    case T::Cone:return std::abs(rho*std::cos(a.halfAngle)-p.z*std::sin(a.halfAngle));
    case T::Sphere:return std::abs(Norm(p)-a.radius);
    default:return std::abs(std::hypot(rho-a.majorRadius,p.z)-a.radius);
    }
}
double Audit(const sggk::BSplineSurface &s,const GeneralRationalBoundaryReport &r,double motherTolerance=1e-8)
{
    Require(r.Passed(),"Cannot audit failed output");double maximum=0;const auto &patch=*r.patch;
    SurfaceSampleContext samples(s);
    for(const auto &piece:patch.tangentExemptions)
    {
        Require(piece.edge<4&&piece.sourceFirst<=piece.sourceLast&&piece.targetFirst<=piece.targetLast&&
            piece.sourceDiameterUpper<=1e-8&&piece.targetDiameterUpper<=1e-8&&piece.distanceUpper<1e-5,
            "Invalid tangent exemption certificate");
        for(unsigned i=0;i<=2;++i)for(unsigned j=0;j<=2;++j)
        {
            unsigned e=piece.edge;double t=piece.sourceFirst+(piece.sourceLast-piece.sourceFirst)*(.5*i);
            auto source=s.CalcPoint(e<2?(e==0?s.MinParamU():s.MaxParamU()):t,
                                   e<2?t:(e==2?s.MinParamV():s.MaxParamV()));
            auto target=patch.boundaries[e]->CalcPoint(piece.targetFirst+(piece.targetLast-piece.targetFirst)*(.5*j));
            Require(Norm(Get(source)-Get(target))<1.001e-5,"Exempt piece failed independent position audit");
        }
    }
    for(unsigned e=0;e<4;++e)
    {
        const auto &curve=*patch.boundaries[e];
        for(unsigned endpoint=0;endpoint<2;++endpoint)
        {
            auto p=Get(curve.ControlPoints()[endpoint?curve.ControlPoints().size()-1:0]);
            auto v=Get(patch.vertices[rb::Ends[e][endpoint]]);
            Require(p.x==v.x&&p.y==v.y&&p.z==v.z,"Output wire endpoint mismatch");
        }
        for(unsigned i=0;i<=32;++i)
        {
            double f=double(i)/32;
            Require(MotherDistance(patch.mother,Get(curve.CalcPoint(f)))<motherTolerance,"Exported boundary off mother");
            double u=e<2?(e==0?s.MinParamU():s.MaxParamU()):s.MinParamU()+(s.MaxParamU()-s.MinParamU())*f;
            double v=e<2?s.MinParamV()+(s.MaxParamV()-s.MinParamV())*f:(e==2?s.MinParamV():s.MaxParamV());
            auto point=samples.At(u,v);double parameter=f;
            auto nearest=curve.CalcNearestPoint(Point(point.point),parameter);
            double distance=Norm(point.point-Get(nearest));maximum=std::max(maximum,distance);
            Require(distance<1.001e-5,"Independent SDK boundary distance audit failed");
            bool exempt=false;
            for(const auto &piece:patch.tangentExemptions)
                if(piece.edge==e&&(e<2?v:u)>=piece.sourceFirst&&(e<2?v:u)<=piece.sourceLast)exempt=true;
            if(!patch.tangentNotApplicable[e]&&!exempt&&i>0&&i<32)
            {
                V a=e<2?point.dv:point.du,b=Get(curve.CalcDeriv1(parameter));
                Require(Dot(a,b)>0&&Norm(Cross(a,b))/Dot(a,b)<std::tan(.1*Pi/180),"Independent tangent audit failed");
            }
        }
    }
    return maximum;
}
void ErrorKernelUnit()
{
    rb::CurveErrors errors;paired::Products products;
    auto value=[](const BPoly &p,long double x)
    {
        std::vector<long double> c;for(I a:p.c)c.push_back(a.lo);
        for(unsigned n=p.u;n>0;--n)for(unsigned i=0;i<n;++i)c[i]=(1-x)*c[i]+x*c[i+1];
        return c[0];
    };
    auto contains=[](I box,long double x){return box.lo<=x&&x<=box.hi;};
    // 对精确二进制系数独立作 de Casteljau 求值，覆盖不同次数以及常权重、变权重两条路径。
    for(unsigned p=1;p<=6;++p)for(unsigned q=1;q<=4;++q)for(unsigned mode=0;mode<4;++mode)
    {
        auto cell=[](unsigned n,bool rational,int sign)
        {
            BezierCell c;c.p={BPoly(n,0),BPoly(n,0),BPoly(n,0)};c.w=BPoly(rational?n:0,0,I(1));
            for(unsigned i=0;i<=n;++i)
            {
                c.p.x.c[i]=I(sign*(int(i)-2)*.25);c.p.y.c[i]=I((i*i+1)*.125);
                c.p.z.c[i]=I(sign*(int(i%3)-1)*.5);
                if(rational)c.w.c[i]=I(1+i*.25);
            }
            return c;
        };
        auto s=cell(p,(mode&1)!=0,1),t=cell(q,(mode&2)!=0,-1);bool fused=false;
        auto e=errors.Build(s,t,products,fused);Require(fused==(mode==3),"Wrong fused error path");
        for(unsigned i=0;i<=8;++i)
        {
            double x=i*.125;IV g=Range(Evaluate(e.numerator,x,0));I h=Evaluate(e.denominator,x,0).range();
            long double sw=value(s.w,x),tw=value(t.w,x);
            Require(contains(h,sw*tw),"Error denominator failed to enclose exact value");
            for(auto axis:{std::make_tuple(&s.p.x,&t.p.x,g.x),std::make_tuple(&s.p.y,&t.p.y,g.y),std::make_tuple(&s.p.z,&t.p.z,g.z)})
                Require(contains(std::get<2>(axis)/h,value(*std::get<1>(axis),x)/tw-value(*std::get<0>(axis),x)/sw),
                        "Fused/constant-weight residual failed independent evaluation");
        }
    }
    auto line=[](double speed,double slope,double offset)
    {BezierCell c;c.p={BPoly(1,0),BPoly(1,0),BPoly(1,0)};c.p.x.c={I(0),I(speed)};c.p.y.c={I(offset),I(offset+slope)};return c;};
    auto ratio=[&](const BezierCell &s,const BezierCell &t)
    {bool fused=false;auto e=errors.Build(s,t,products,fused);return rb::ResidualTangentRatio(e,s,rb::TangentNumerator(s,products),paired::PositionBoxUpper(e.numerator,e.denominator));};
    auto s=line(1,0,0);
    double close=ratio(s,line(1,1e-4,0));
    Require(close>=1e-4&&close<std::tan(.1*Pi/180),"Close tangent bound invalid");
    Require(ratio(s,line(1,.01,0))>=.01,"Angular violation understated");
    Require(ratio(s,line(1,0,1e-6))<1e-12,"Translation should have near-zero derivative error");
    auto fast=line(2,0,0);Require(!std::isfinite(ratio(s,fast)),"Speed mismatch must be inconclusive");
    auto st=rb::TangentNumerator(s,products),tt=rb::TangentNumerator(fast,products);
    auto speedScale=rb::ChooseTangentScale(st,tt);Require(speedScale.has_value(),"Missing positive speed scale");
    Require(rb::ScaledTangentSineSquaredUpper(st,rb::BuildScaledTangentResidual(st,tt,*speedScale,products),*speedScale)<1e-24,
            "Scaled residual should accept parallel speed mismatch");
    Require(!std::isfinite(ratio(s,line(-1,0,0))),"Opposite direction accepted by residual criterion");
    Require(!std::isfinite(ratio(line(0,0,0),s)),"Zero source speed accepted");
    // 用权重 W=1+t 的有理二次形式表示同一直线。
    BezierCell rational;rational.p={BPoly(2,0),BPoly(2,0),BPoly(2,0)};
    rational.p.x.c={I(0),I(.5),I(2)};rational.w=BPoly(2,0);rational.w.c={I(1),I(1.5),I(2)};
    auto tilted=rational;tilted.p.y=rational.p.x*I(1e-5);
    double rationalRatio=ratio(rational,tilted);
    Require(rationalRatio>=1e-5&&rationalRatio<std::tan(.1*Pi/180),"Rational derivative bound invalid");
    for(double scale:{1./1024,1024.})
    {
        auto scaled=rational;scaled.p=scaled.p*I(scale);scaled.w=scaled.w*I(scale);
        double scaledRatio=ratio(scaled,tilted);
        Require(scaledRatio>=1e-5&&scaledRatio<std::tan(.1*Pi/180),"Homogeneous scaling changed tangent certificate");
    }
    auto invalid=s;invalid.w.c[0]=I(-1,1);
    Require(!std::isfinite(ratio(invalid,s)),"Nonpositive denominator accepted");
    rb::CurveEnvelopes envelopes;
    auto scaledBound=[&](const BV &a,const BV &b,bool envelope)->double
    {
        auto lambda=rb::ChooseTangentScale(a,b);if(!lambda)return INFINITY;
        auto residual=envelope?envelopes.TangentResidual(a,b,*lambda):
            std::optional<BV>(rb::BuildScaledTangentResidual(a,b,*lambda,products));
        return residual?rb::ScaledTangentSineSquaredUpper(a,*residual,*lambda):INFINITY;
    };
    auto sine2=[](double tangent){return tangent*tangent/(1+tangent*tangent);};
    auto parallel=envelopes.Tangent(s),scaled=envelopes.Tangent(fast),opposite=envelopes.Tangent(line(-1,0,0));
    Require(parallel&&scaled&&opposite,"Missing tangent envelopes");
    Require(scaledBound(*parallel,*scaled,true)<1e-24,"Parallel envelope tangent failed");
    Require(!std::isfinite(scaledBound(*parallel,*opposite,true)),"Opposite envelope tangent accepted");
    Require(envelopes.TangentExceeded(*parallel,*opposite,sine2(.001))&&
            !envelopes.TangentExceeded(*parallel,*scaled,sine2(.001)),"Wrong tangent midpoint witness");
    auto tiltedTangent=envelopes.Tangent(line(1,1e-4,0));
    double envelopeSquared=scaledBound(*parallel,*tiltedTangent,true);
    Require(envelopeSquared>=sine2(1e-4)&&envelopeSquared<sine2(std::tan(.1*Pi/180)),"Envelope tangent bound understated");
    Require(envelopes.TangentExceeded(*parallel,*tiltedTangent,sine2(.99e-4))&&
            !envelopes.TangentExceeded(*parallel,*tiltedTangent,sine2(1.01e-4)),"Tangent witness crossed tolerance");
    BV direction{BPoly(1,0,I(1)),BPoly(1,0),BPoly(1,0)};direction.y.c={I(0),I(.2)};
    BPoly varyingScale(2,0);varyingScale.c={I(1),I(2),I(4)};
    auto variable=products.MulShared(direction,varyingScale);
    for(bool envelope:{false,true})
    {
        Require(scaledBound(variable,direction,envelope)<1e-20,"Quadratic positive speed scale rejected");
        BV zero{BPoly(2,0),BPoly(2,0),BPoly(2,0)};zero.x.c={I(1),I(-1),I(1)};
        Require(!std::isfinite(scaledBound(zero,*parallel,envelope)),"Vanishing interior source tangent accepted");
        Require(!std::isfinite(scaledBound(*parallel,BV{},envelope)),"Zero target tangent accepted");
        auto gap=*parallel;gap.y=BPoly(4,0);gap.y.c={I(0),I(.04),I(-.08/1.5),I(.04),I(0)};
        Require(!(scaledBound(*parallel,gap,envelope)<sine2(.001)),"Between-sample angular defect accepted");
        Require(!envelopes.TangentExceeded(*parallel,gap,sine2(.001)),"Residual failure mislabeled as midpoint witness");
    }
    Require(!envelopes.TangentResidual(*parallel,*parallel,BPoly(0,0,I(-1))),"Negative auxiliary scale accepted");
    auto near=envelopes.Build(s,line(1,0,.99e-5)),far=envelopes.Build(s,line(1,0,1.01e-5));
    Require(near&&far&&!envelopes.PositionExceeded(*near,1e-5)&&envelopes.PositionExceeded(*far,1e-5),
            "Envelope midpoint witness crossed tolerance");
    auto uncertain=*near;uncertain.numerator.y.c[0]=I(0,2e-5);
    Require(!envelopes.PositionExceeded(uncertain,1e-5),"Uncertain midpoint incorrectly excluded");
    auto invalidWeight=*far;invalidWeight.denominator.c[0]=I(-1,1);
    Require(!envelopes.PositionExceeded(invalidWeight,1e-5),"Nonpositive witness denominator accepted");
}
void EnvelopeExport(const std::filesystem::path &path)
{
    std::ofstream out(path);Require(bool(out),"Cannot open envelope audit output");
    out<<std::setprecision(17)<<'[';rb::CurveEnvelopes envelope;
    auto poly=[&](const BPoly &p)
    {
        out<<'[';for(size_t i=0;i<p.c.size();++i){if(i)out<<',';out<<'['<<p.c[i].lo<<','<<p.c[i].hi<<']';}out<<']';
    };
    auto cell=[&](const BezierCell &c)
    {out<<'[';poly(c.p.x);out<<',';poly(c.p.y);out<<',';poly(c.p.z);out<<',';poly(c.w);out<<']';};
    unsigned count=0;
    for(unsigned p:{1u,3u,5u,7u,12u})for(unsigned q:{1u,2u,4u,8u})for(unsigned mode=0;mode<8;++mode)
    {
        auto make=[&](unsigned n,int sign,bool rational)
        {
            BezierCell c;c.p={BPoly(n,0),BPoly(n,0),BPoly(n,0)};c.w=BPoly(rational?n:0,0,I(1));
            double scale=mode==4?0x1p-500:mode==5?0x1p400:mode==6?0x1p-1030:mode==7?0x1p-600:1;
            for(unsigned i=0;i<=n;++i)
            {
                unsigned k=0;for(auto axis:{&c.p.x,&c.p.y,&c.p.z})
                {
                    double x=scale*(sign*(int((i*17+3*k)%23)-11)+.125*i);
                    axis->c[i]=mode==3?I(std::nextafter(x,-INFINITY),std::nextafter(x,INFINITY)):I(x);++k;
                }
                if(rational)c.w.c[i]=I((1+(i%5)*.25)*(mode==7?0x1p-500:1));
            }
            return c;
        };
        auto s=make(p,1,(mode&1)!=0),t=make(q,-1,(mode&2)!=0);
        auto error=envelope.Build(s,t);auto tangent=envelope.Tangent(s);
        Require(error.has_value()&&tangent.has_value(),"Finite envelope unavailable");
        if(count++)out<<',';out<<"{\"source\":";cell(s);out<<",\"target\":";cell(t);
        out<<",\"error\":";cell({error->numerator,0,1,0,1,error->denominator});
        out<<",\"tangent\":";cell({*tangent});
        auto direction=[](unsigned n,double slope)
        {
            BV v{BPoly(n,0),BPoly(n,0),BPoly(n,0)};
            for(unsigned i=0;i<=n;++i){v.x.c[i]=I(1+i*.125);v.y.c[i]=I(slope*(1+i*.125));v.z.c[i]=I((i%3)*1e-5);}
            return v;
        };
        auto a=direction(p,.001),b=direction(q,.0011);
        auto lambda=rb::ChooseTangentScale(a,b);Require(lambda.has_value(),"Missing residual scale");
        auto residual=envelope.TangentResidual(a,b,*lambda);Require(residual.has_value(),"Missing residual envelope");
        double q2=rb::ScaledTangentSineSquaredUpper(a,*residual,*lambda);Require(std::isfinite(q2),"Missing residual angle bound");
        V e=rb::TangentResidualReference(a);
        out<<",\"angle_source\":";cell({a});out<<",\"angle_target\":";cell({b});
        out<<",\"angle_scale\":";poly(*lambda);out<<",\"angle_residual\":";cell({*residual});
        out<<",\"angle_reference\":["<<e.x<<','<<e.y<<','<<e.z<<']';
        out<<",\"angle_sine_squared\":"<<q2<<'}';
    }
    // 非法或非有限算术结果必须回退到原区间路径，不能视为通过。
    BezierCell invalid;invalid.p.x=BPoly(1,1);
    Require(!envelope.Build(invalid,invalid),"Tensor input accepted by 1D envelope");
    invalid={};invalid.p.x.c[0]=I::Whole();
    Require(!envelope.Build(invalid,invalid),"Nonfinite input accepted by envelope");
    out<<"]\n";std::cout<<"Envelope export cases="<<count<<'\n';
}
void RecoveryRepairUnit()
{
    using P=sggk::Point2D;
    std::array<P,4> quad{P{0,0},P{1,0},P{0,1},P{1,1}};
    std::array<bool,4> collapsed{};
    Require(rb::ValidRationalPolygon(quad,collapsed),"Regular rational quadrilateral rejected");
    for(unsigned edge=0;edge<4;++edge)
    {
        auto triangle=quad;triangle[rb::Ends[edge][1]]=triangle[rb::Ends[edge][0]];
        collapsed={};collapsed[edge]=true;
        Require(rb::ValidRationalPolygon(triangle,collapsed),"Exact adjacent triangle rejected");
        Require(!rb::ValidRationalPolygon(triangle,{}),"UV collapse without point source edge accepted");
    }
    collapsed.fill(true);
    Require(!rb::ValidRationalPolygon({P{0,0},P{0,0},P{0,0},P{0,0}},collapsed),"Point domain accepted");
    Require(!rb::ValidRationalPolygon({P{0,0},P{1,0},P{0,0},P{2,0}},collapsed),"Collinear triangle accepted");
    Require(!rb::ValidRationalPolygon({P{0,0},P{1,1},P{0,1},P{1,0}},collapsed),"Folded quadrilateral accepted");
    Require(!rb::ValidRationalPolygon({P{0,0},P{1,0},P{0,1},P{0,0}},collapsed),"Nonadjacent duplicate accepted");
    auto line=[](double length,double offset,double lo=0,double hi=1)
    {
        BezierCell c;c.p={BPoly(1,0),BPoly(1,0,I(offset)),BPoly(1,0)};
        c.p.x.c={I(0),I(length)};c.u0=lo;c.u1=hi;return c;
    };
    auto small=line(4e-9,0),near=line(4e-9,2e-6);
    auto merge=rb::MicroSegmentBound(small,{near},0,1,1e-8,1e-5);
    Require(merge&&merge->distanceUpper>=2e-6&&merge->sourceDiameterUpper<=1e-8,"Micro spatial certificate missing");
    Require(!rb::MicroSegmentBound(small,{line(4e-9,2e-5)},0,1,1e-8,1e-5),"Micro merge waived position violation");
    Require(!rb::MicroSegmentBound(line(1e-4,0),{line(1e-4,0)},0,1,1e-8,1e-5),"Long coincident curves exempted");
    Require(!rb::MicroSegmentBound(small,{near},0,1,0,1e-5),"Disabled micro merge accepted");
    Require(!rb::MicroSegmentBound(small,{line(2e-9,0,0,.4),line(2e-9,0,.6,1)},0,1,1e-8,1e-5),"Micro merge omitted target gap");
    auto point=paired::Crop(small,0,0,0,1);
    Require(rb::MicroSegmentBound(point,{near},0,1,1e-8,1e-5).has_value(),"Endpoint bridge not certified");
    rb::Trim trim;trim.mother=Mother(3);trim.domain={.2,.3,1,100};
    AnalyticPatch patch;patch.parameters=trim.mother;patch.domain=trim.domain;
    trim.chart=paired::BuildChart(patch,false);paired::PlaceChart(trim.chart,trim.mother,Get(trim.mother.origin));
    trim.polygon.corners={P{0,0},P{1,.3},P{0,1},P{1,1}};
    V expected=rb::ChartPoint(trim,{.4,.12});
    double f=rb::RefineFraction(trim,2,expected,.401);
    Require(f>=0&&f<=1&&Norm(expected-rb::ChartPoint(trim,{f,.3*f}))<1e-10,"Bounded 3D correspondence failed");
}
void Unit(const std::filesystem::path &data)
{
    ErrorKernelUnit();
    RecoveryRepairUnit();
    const std::array<double,4> domain{.2,.9,.6,1.2};
    for(unsigned type=0;type<5;++type)
    {
        auto a=Mother(type);auto s=Fixture(a,domain,true);GeneralRationalBoundaryOptions o;o.chartDomain=domain;
        auto r=ValidateGeneralRationalBoundaries(*s,a,o);
        Require(r.Passed(),"Rational skew fixture rejected type "+std::to_string(type)+": "+r.message);Audit(*s,r);
        Require((type==0||r.envelopeAccepted>0)&&r.envelopeChecks==r.envelopeAccepted+r.envelopeFallbacks+r.envelopePositionWitnesses+r.envelopeTangentWitnesses,"Envelope filter not exercised");
        Require(r.maxTargetDegree<=4&&!r.certifiesInteriorContinuum,"Wrong certificate contract");
        s->UReverse();auto reverse=ValidateGeneralRationalBoundaries(*s,a,o);
        Require(reverse.Passed(),"Reversed rational boundary rejected");Audit(*s,reverse);
        o.maxBoundaryNodes=1;auto budget=ValidateGeneralRationalBoundaries(*s,a,o);
        Require(budget.status==HybridValidationStatus::BudgetExceeded&&!budget.patch,"Budget leaked output");
    }
    auto plane=Mother(0);auto flat=Fixture(plane,domain,false);
    {
        GeneralRationalBoundaryOptions deadlineOptions;deadlineOptions.chartDomain=domain;
        SurfaceSampleContext samples(*flat);
        ExecutionDeadline expired(50,ExecutionDeadline::Clock::now()-std::chrono::seconds(1));
        samples.deadline=&expired;
        auto timed=rb::Validate(samples,plane,deadlineOptions);
        Require(timed.timeBudgetExceeded&&timed.status==HybridValidationStatus::BudgetExceeded&&!timed.patch,
                "Expired shared budget reset or leaked output");
        Require(samples.deadline==&expired,"Validator replaced the caller deadline");
        GeneralClassificationDiagnostics d;bool stopped=false;
        try{InitializeGeneralCandidates(samples,deadlineOptions.general,d);}
        catch(const DeadlineExceeded &){stopped=true;}
        Require(stopped,"Classification swallowed deadline cancellation");
        deadlineOptions.maxTotalMs=-1;
        Require(ValidateGeneralRationalBoundaries(*flat,plane,deadlineOptions).status==HybridValidationStatus::InvalidInput,
                "Invalid time budget accepted");
        deadlineOptions.maxTotalMs=0;
        Require(ValidateGeneralRationalBoundaries(*flat,plane,deadlineOptions).Passed(),"Explicit diagnostic no-deadline mode failed");
    }
    auto variableSpeed=Fixture(plane,domain,false);variableSpeed->UDegreeElevation(3-variableSpeed->DegreeU());
    for(unsigned j=0;j<variableSpeed->ControlPoints().front().size();++j)
        variableSpeed->SetControlPoint(Point(Get(variableSpeed->ControlPoints()[1][j])+Get(plane.xDirection)*1e-6),1,j);
    GeneralRationalBoundaryOptions speedOptions;speedOptions.chartDomain=domain;
    speedOptions.general.fitting.validation.boundaryTangentToleranceDeg=1e-6;
    auto speedReport=ValidateGeneralRationalBoundaries(*variableSpeed,plane,speedOptions);
    Require(speedReport.Passed()&&speedReport.tangentFallbacks>0,"Parallel variable-speed boundary must pass via full tangent fallback");
    GeneralRationalBoundaryOptions bad;bad.general.fitting.validation.distanceTolerance=-1;
    Require(ValidateGeneralRationalBoundaries(*flat,plane,bad).status==HybridValidationStatus::InvalidInput,"Invalid tolerance accepted");
    bad={};bad.microSegmentDiameter=1e-7;
    Require(ValidateGeneralRationalBoundaries(*flat,plane,bad).status==HybridValidationStatus::InvalidInput,"Oversized tangent exemption accepted");
    auto wrong=plane;wrong.origin=Point(Get(plane.origin)+Get(plane.axis)*.01);
    Require(!ValidateGeneralRationalBoundaries(*flat,wrong).Passed(),"Wrong mother accepted");
    auto bump=Fixture(plane,domain,false);bump->UDegreeElevation(3-bump->DegreeU());bump->VDegreeElevation(3-bump->DegreeV());
    for(unsigned i=0;i<4;++i)for(unsigned j=1;j<3;++j)
        bump->SetControlPoint(Point(Get(bump->ControlPoints()[i][j])+Get(plane.axis)*1e-7),i,j);
    GeneralRationalBoundaryOptions tight;tight.general.fitting.validation.boundaryTangentToleranceDeg=1e-8;
    auto tangentFailure=ValidateGeneralRationalBoundaries(*bump,plane,tight);
    Require(!tangentFailure.Passed()&&tangentFailure.tangentFallbacks>0,"Tangent violation must remain rejected through full fallback");
    auto outside=Fixture(plane,domain,false);outside->UDegreeElevation(3-outside->DegreeU());outside->VDegreeElevation(3-outside->DegreeV());
    for(unsigned i=1;i<3;++i)for(unsigned j=1;j<3;++j)
        outside->SetControlPoint(Point(Get(outside->ControlPoints()[i][j])+Get(plane.xDirection)*10),i,j);
    Require(!ValidateGeneralRationalBoundaries(*outside,plane).Passed(),"Interior finite-region excursion accepted");
    for(unsigned id:{0u,218u,318u})
    {
        auto f=ReadGeneralFixture(data/"main"/RationalFixtureName("main",id));
        Require(f.rationalDomain&&CertifyRationalDomain(f).passed,"Unit fixture requires a certified V2 source");
        GeneralRationalBoundaryOptions fixed;fixed.chartDomain=f.chartDomain;auto r=ValidateGeneralRationalBoundaries(*f.input,f.truth,fixed);
        Require(r.Passed(),"Known rectangular/collapsed fixture failed "+std::to_string(id)+": "+r.message);Audit(*f.input,r);
        if(id)Require(r.patch->tangentNotApplicable[0]||r.patch->tangentNotApplicable[1]||
            r.patch->tangentNotApplicable[2]||r.patch->tangentNotApplicable[3],"Collapsed source tangent not marked");
    }
    // 检查跨内结点的分段参数图复合，参考测试独立实现的圆参数逆公式，并覆盖周期接缝。
    for(unsigned type:{1u,2u,3u,4u})
    {
        rb::Trim t;t.mother=Mother(type);t.domain={5.8,8.2,.3,2.2};
        AnalyticPatch a;a.parameters=t.mother;a.domain=t.domain;t.chart=paired::BuildChart(a);
        paired::PlaceChart(t.chart,t.mother,Get(t.mother.origin));
        t.polygon.corners={sggk::Point2D{0,0},sggk::Point2D{1,.2},sggk::Point2D{.1,1},sggk::Point2D{.9,.8}};
        paired::Products products;auto edges=rb::ComposeEdge(t,2,products);Require(edges.size()>1,"Missing chart cuts");
        for(const auto &c:edges)for(unsigned i=0;i<=20;++i)
        {
            double x=double(i)/20,f=c.u0+(c.u1-c.u0)*x;std::array<double,2> uv{f,.2*f},angle{};
            for(unsigned d=0;d<2;++d)
            {
                const auto &knots=t.chart.angles[d];unsigned n=unsigned(knots.size()-1);
                unsigned k=std::min(n-1,unsigned(uv[d]*n));double q=uv[d]*n-k;
                if(t.chart.angular[d]){double h=(knots[k+1]-knots[k])*.5;angle[d]=knots[k]+2*std::atan2(q*std::sin(h),1-q+q*std::cos(h));}
                else angle[d]=knots[k]+q*(knots[k+1]-knots[k]);
            }
            IV p=Range(Evaluate(c.p,x,0));I w=Evaluate(c.w,x,0).range();
            V actual{Mid(p.x/w),Mid(p.y/w),Mid(p.z/w)};
            Require(Norm(actual-FrameOf(t.mother).vector(CanonicalPoint(t.mother,angle[0],angle[1])))<1e-10,"Rational line composition mismatch");
        }
        GeneralRationalBoundaryReport report;auto target=rb::BuildTarget(t,1e-7,report);
        Require(target.motherBound<1e-7,"Multi-cell export bound failed");
        for(const auto &c:target.patch.boundaries)for(unsigned i=0;i<=64;++i)
            Require(MotherDistance(t.mother,Get(c->CalcPoint(double(i)/64)))<1e-8,"Multi-cell exported curve off mother");
    }
    std::cout<<"Rational boundary unit assertions="<<assertions<<'\n';
}
void FullCase(const std::filesystem::path &data,const std::filesystem::path &csv,
              const std::string &suite,unsigned id,unsigned repeats,bool forceFixed=false,double maxCaseMs=50)
{
    bool synthetic=suite=="synthetic",extra=suite=="extra",extreme=suite=="extreme";
    Require(synthetic||extra||extreme||suite=="main","Unknown rational suite");
    Require(repeats>0&&id<(synthetic?5u:extra?60u:extreme?320u:600u),"Invalid full-suite range");
    GeneralFixture c=ReadGeneralFixture(data/suite/RationalFixtureName(suite,id));
    Require(c.rationalDomain&&c.id==id,"Current tests require rational-domain V2 inputs");
    std::string group=FixtureGroupName(c.group);double certificateMs=0;
    std::cout<<"BEGIN_INPUT_CERT"<<std::endl;auto certificateStart=std::chrono::steady_clock::now();
    auto cert=CertifyRationalDomain(c);certificateMs=Ms(certificateStart);
    Require(cert.passed,"V2 input certificate failed: "+cert.message);
    GeneralRationalBoundaryOptions options;
    options.maxTotalMs=maxCaseMs;
    bool supplied=forceFixed||synthetic||(extra&&(id<20||id>=40));
    if(forceFixed||synthetic||(extra&&id<20))options.chartDomain=c.chartDomain;
    if(extra&&c.group==GeneralFixtureGroup::InvalidOptions)
    {if(id%2)options.general.fitting.validation.distanceTolerance=-1;else options.general.fitting.sampling.maxSamples=0;}
    if(extra&&c.group==GeneralFixtureGroup::ValidationBudget)options.maxBoundaryNodes=1;
    if(extra&&c.group==GeneralFixtureGroup::SpherePriority)options.general.similarScoreAbsoluteRatio=1e6;
    const std::string header="suite,id,run,mode,group,source_type,rational,degree_u,degree_v,controls,source_certified,source_bound,certificate_ms,entry,status,accepted,audit_status,audit_sample_max,classify_ms,validate_ms,sampling_ms,recovery_ms,prepare_ms,boundary_ms,total_ms,attempts,chart_attempts,nodes,refinements,fast_checks,fast_accepted,tangent_fallbacks,position_only,fused_errors,boundary_bound,tangent_deg,mother_bound,output_type,raw_candidates,full_scores,message,audit_message,envelope_checks,envelope_accepted,envelope_fallbacks,envelope_position_witnesses,envelope_tangent_witnesses,tangent_residual_checks,tangent_residual_accepted,tangent_residual_interval_fallbacks,max_case_ms,time_budget_exceeded,timeout_stage,within_time_limit,micro_merged_segments,tangent_exemption_count,recovery_3d_refinements\n";
    std::vector<std::string> rows;bool anyFailure=false;
    auto save=[&]{std::ofstream out(csv);Require(bool(out),"Cannot write full-suite CSV");out<<header;for(const auto &row:rows)out<<row;};
    auto attemptPath=csv.parent_path()/(csv.stem().string()+"-attempts.csv");
    std::ofstream attemptsFile(attemptPath);Require(bool(attemptsFile),"Cannot write attempt diagnostics");
    attemptsFile<<"suite,id,run,attempt,type,status,accepted,charts,nodes,refinements,sampling_ms,recovery_ms,prepare_ms,boundary_ms,total_ms,message\n"<<std::setprecision(17);
    std::ofstream exemptionFile(csv.parent_path()/(csv.stem().string()+"-exemptions.csv"));
    exemptionFile<<"suite,id,run,edge,source_first,source_last,target_first,target_last,source_diameter,target_diameter,distance_upper\n"<<std::setprecision(17);
    struct Attempt{sggk::SurfaceType type;GeneralRationalBoundaryReport report;};
    for(unsigned run=0;run<repeats;++run)
    {
        std::cout<<"BEGIN_ALGORITHM "<<run<<std::endl;
        auto start=std::chrono::steady_clock::now();GeneralRationalBoundaryReport r;
        GeneralClassificationDiagnostics classification;std::vector<Attempt> attempts;
        double classify=0,validate=0,sampling=0,recovery=0,prepare=0,boundary=0;
        unsigned charts=0,nodes=0,refinements=0,fast=0,fastAccepted=0,fallback=0,position=0,fused=0;
        unsigned envelopeChecks=0,envelopeAccepted=0,envelopeFallbacks=0,envelopeWitnesses=0,envelopeTangentWitnesses=0;
        unsigned residualChecks=0,residualAccepted=0,residualFallbacks=0;
        unsigned microMerged=0,refined3D=0;
        try
        {
            ExecutionDeadline deadline(maxCaseMs,start);
            SurfaceSampleContext samples(*c.input,options.general.fitting.sampling.maxSamples);
            samples.deadline=&deadline;
            auto attempt=[&](const AnalyticParameters &a)
            {
                r=rb::Validate(samples,a,options);validate+=r.totalMs;sampling+=r.samplingMs;recovery+=r.recoveryMs;
                prepare+=r.preparationMs;boundary+=r.boundaryMs;charts+=r.chartAttempts;nodes+=r.boundaryNodes;
                refinements+=r.mappingRefinements;fast+=r.tangentFastChecks;fastAccepted+=r.tangentFastAccepted;
                fallback+=r.tangentFallbacks;position+=r.positionOnlyNodes;fused+=r.fusedErrorBuilds;
                envelopeChecks+=r.envelopeChecks;envelopeAccepted+=r.envelopeAccepted;envelopeFallbacks+=r.envelopeFallbacks;envelopeWitnesses+=r.envelopePositionWitnesses;envelopeTangentWitnesses+=r.envelopeTangentWitnesses;
                residualChecks+=r.tangentResidualChecks;residualAccepted+=r.tangentResidualAccepted;residualFallbacks+=r.tangentResidualIntervalFallbacks;
                microMerged+=r.microMergedSegments;refined3D+=r.recovery3DRefinements;
                attempts.push_back({a.type,r});
            };
            if(supplied)attempt(c.truth);
            else
            {
                GeneralCandidatePool pool;
                {StageTimer timer(classify);pool=InitializeGeneralCandidates(samples,options.general,classification);}
                for(;;)
                {
                    std::vector<GeneralSeed> batch;
                    {StageTimer timer(classify);batch=NextGeneralCandidateBatch(pool,samples,options.general,classification);}
                    if(batch.empty())break;
                    for(const auto &candidate:batch){attempt(candidate.parameters);if(r.Passed()||r.timeBudgetExceeded)break;}
                    if(r.Passed()||r.timeBudgetExceeded)break;
                }
                if(attempts.empty())r.message="No classified candidate available";
            }
        }
        catch(const DeadlineExceeded &e)
        {r.patch.reset();r.status=HybridValidationStatus::BudgetExceeded;r.timeBudgetExceeded=true;r.timeoutStage=e.stage;r.message=e.what();}
        catch(const std::exception &e){r.patch.reset();r.status=HybridValidationStatus::NotCertified;r.message=std::string("Algorithm exception: ")+e.what();}
        double total=Ms(start);
        if(maxCaseMs>0&&total>maxCaseMs)
        {r.patch.reset();r.status=HybridValidationStatus::BudgetExceeded;r.timeBudgetExceeded=true;r.timeoutStage="case_completion_overrun";r.message="Case wall time exceeded hard limit";}
        bool passed=r.Passed();anyFailure|=!passed;
        std::cout<<"END_ALGORITHM "<<run<<std::endl;
        std::string auditStatus=passed?"pending":"not_applicable",auditMessage;double audit=std::numeric_limits<double>::quiet_NaN();
        auto line=[&]
        {
            std::ostringstream out;out<<std::setprecision(17)<<suite<<','<<id<<','<<run<<",rational,"<<group<<','<<(c.type<5?TypeName(c.truth.type):"Freeform")<<','<<c.input->IsRational()<<','
                <<c.input->DegreeU()<<','<<c.input->DegreeV()<<','<<c.input->ControlPoints().size()*c.input->ControlPoints().front().size()<<','<<int(cert.passed)<<','<<cert.upper<<','<<certificateMs<<','
                <<(options.chartDomain?"fixed_chart":supplied?"supplied_mother":"classified")<<','<<HybridValidationStatusName(r.status)<<','<<passed<<','<<auditStatus<<','<<audit<<','
                <<classify<<','<<validate<<','<<sampling<<','<<recovery<<','<<prepare<<','<<boundary<<','<<total<<','<<attempts.size()<<','<<charts<<','<<nodes<<','<<refinements<<','
                <<fast<<','<<fastAccepted<<','<<fallback<<','<<position<<','<<fused<<','<<r.boundaryDistanceUpperBound<<','<<r.boundaryTangentUpperBoundDeg<<','<<r.outputMotherDistanceUpperBound<<','
                <<(passed?TypeName(r.patch->mother.type):"Unknown")<<','<<classification.rawCandidateCount<<','<<classification.fullScoreCount<<','<<Quote(r.message)<<','<<Quote(auditMessage)<<','<<envelopeChecks<<','<<envelopeAccepted<<','<<envelopeFallbacks<<','<<envelopeWitnesses<<','<<envelopeTangentWitnesses<<','<<residualChecks<<','<<residualAccepted<<','<<residualFallbacks
                <<','<<maxCaseMs<<','<<r.timeBudgetExceeded<<','<<Quote(r.timeoutStage)<<','<<(maxCaseMs<=0||total<=maxCaseMs)<<','<<microMerged<<','<<(passed?r.patch->tangentExemptions.size():0)<<','<<refined3D<<'\n';
            return out.str();
        };
        // 独立 SDK 审计前保存已测算法耗时，避免审计看门狗超时覆盖已经完成的算法结果。
        rows.push_back(line());save();
        for(size_t i=0;i<attempts.size();++i)
        {
            const auto &a=attempts[i];const auto &v=a.report;
            attemptsFile<<suite<<','<<id<<','<<run<<','<<i<<','<<TypeName(a.type)<<','<<HybridValidationStatusName(v.status)<<','<<v.Passed()<<','<<v.chartAttempts<<','<<v.boundaryNodes<<','<<v.mappingRefinements<<','
                <<v.samplingMs<<','<<v.recoveryMs<<','<<v.preparationMs<<','<<v.boundaryMs<<','<<v.totalMs<<','<<Quote(v.message)<<'\n';
        }
        attemptsFile.flush();
        if(passed)for(const auto &piece:r.patch->tangentExemptions)
            exemptionFile<<suite<<','<<id<<','<<run<<','<<piece.edge<<','<<piece.sourceFirst<<','<<piece.sourceLast<<','
                <<piece.targetFirst<<','<<piece.targetLast<<','<<piece.sourceDiameterUpper<<','<<piece.targetDiameterUpper<<','<<piece.distanceUpper<<'\n';
        exemptionFile.flush();
        if(passed)
        {
            std::cout<<"BEGIN_AUDIT "<<run<<std::endl;
            try{audit=Audit(*c.input,r,1e-6);auditStatus="passed";}
            catch(const std::exception &e){auditStatus="unresolved";auditMessage=e.what();}
            rows.back()=line();save();std::cout<<"END_AUDIT "<<run<<std::endl;
        }
    }
    // 真值诊断仅在全部计时运行后执行，不反馈到候选池。
    if(anyFailure&&((suite=="main"&&id<500)||extreme))
    {
        std::cout<<"BEGIN_KNOWN_MOTHER_DIAGNOSTIC"<<std::endl;
        GeneralRationalBoundaryOptions knownOptions;knownOptions.chartDomain=c.chartDomain;knownOptions.maxTotalMs=maxCaseMs;
        auto known=ValidateGeneralRationalBoundaries(*c.input,c.truth,knownOptions);
        std::ofstream out(csv.parent_path()/(csv.stem().string()+"-known.csv"));
        out<<"suite,id,status,accepted,total_ms,charts,nodes,message\n"<<std::setprecision(17)<<suite<<','<<id<<','<<HybridValidationStatusName(known.status)<<','<<known.Passed()<<','
            <<known.totalMs<<','<<known.chartAttempts<<','<<known.boundaryNodes<<','<<Quote(known.message)<<'\n';
        std::cout<<"END_KNOWN_MOTHER_DIAGNOSTIC"<<std::endl;
    }
    if(extra&&id==49)
    {
        // 重新测量重建后的低次内部探针；新的 Greville 点未必保留历史采样盲区，不能照搬旧结论。
        SurfaceSampleContext samples(*c.input);double sampled=0,peak=0;
        for(const auto &p:samples.Greville(false))sampled=std::max(sampled,MotherDistance(c.truth,p.point));
        for(unsigned i=0;i<=10000;++i)peak=std::max(peak,MotherDistance(c.truth,Get(c.input->CalcPoint(double(i)/10000,.5))));
        std::ofstream out(csv.parent_path()/(csv.stem().string()+"-interior.csv"));out<<"greville_max,dense_witness_max\n"<<std::setprecision(17)<<sampled<<','<<peak<<'\n';
    }
}
void Case(const std::filesystem::path &data,const std::filesystem::path &csv,unsigned id,unsigned repeats,bool synthetic)
{
    const std::string suite=synthetic?"synthetic":"main";
    auto fixture=ReadGeneralFixture(data/suite/RationalFixtureName(suite,id));
    Require(fixture.rationalDomain,"Retired dataset: V2 required");
    auto cert=CertifyRationalDomain(fixture);Require(cert.passed,"V2 source certificate failed: "+cert.message);
    auto domain=fixture.chartDomain;
    std::ofstream out(csv);out<<std::setprecision(17);
    out<<"id,synthetic,mode,run,status,accepted,attempts,chart_attempts,nodes,refinements,max_degree,sampling_ms,recovery_ms,prepare_ms,boundary_ms,validate_ms,classify_ms,pipeline_ms,boundary_bound,tangent_deg,mother_bound,audit_sample_max,message,fused_errors,position_only,fast_checks,fast_accepted,tangent_fallbacks,source_degree_u,source_degree_v,source_controls,tangent_residual_checks,tangent_residual_accepted,tangent_residual_interval_fallbacks\n";
    for(unsigned run=0;run<repeats;++run)
    {
        auto start=std::chrono::steady_clock::now();GeneralRationalBoundaryOptions options;
        if(synthetic)options.chartDomain=domain;
        SurfaceSampleContext samples(*fixture.input,options.general.fitting.sampling.maxSamples);
        GeneralRationalBoundaryReport r;unsigned attempts=0,charts=0,nodes=0,refinements=0;
        unsigned fused=0,positionOnly=0,fastChecks=0,fastAccepted=0,fallbacks=0;
        unsigned residualChecks=0,residualAccepted=0,residualFallbacks=0;
        double classify=0,validate=0,sampling=0,recovery=0,prepare=0,boundary=0;bool passed=false;
        auto attempt=[&](const AnalyticParameters &a)
        {
            ++attempts;
            r=rb::Validate(samples,a,options);passed=r.Passed();validate+=r.totalMs;sampling+=r.samplingMs;
            recovery+=r.recoveryMs;prepare+=r.preparationMs;boundary+=r.boundaryMs;
            charts+=r.chartAttempts;nodes+=r.boundaryNodes;refinements+=r.mappingRefinements;
            fused+=r.fusedErrorBuilds;positionOnly+=r.positionOnlyNodes;fastChecks+=r.tangentFastChecks;
            fastAccepted+=r.tangentFastAccepted;fallbacks+=r.tangentFallbacks;
            residualChecks+=r.tangentResidualChecks;residualAccepted+=r.tangentResidualAccepted;residualFallbacks+=r.tangentResidualIntervalFallbacks;
        };
        if(synthetic)attempt(fixture.truth); // 这是显式给定参数图的验证器基准，不是分类耗时。
        else
        {
            GeneralClassificationDiagnostics d;auto t=std::chrono::steady_clock::now();auto pool=InitializeGeneralCandidates(samples,options.general,d);classify+=Ms(t);
            for(;;)
            {
                t=std::chrono::steady_clock::now();auto batch=NextGeneralCandidateBatch(pool,samples,options.general,d);classify+=Ms(t);
                if(batch.empty())break;
                for(const auto &c:batch){attempt(c.parameters);if(passed)break;}
                if(passed)break;
            }
        }
        double pipeline=Ms(start),audit=0;
        if(passed)audit=Audit(*fixture.input,r);
        out<<id<<','<<synthetic<<",rational,"<<run<<','<<HybridValidationStatusName(r.status)<<','<<passed<<','<<attempts<<','<<charts<<','<<nodes<<','<<refinements<<','<<r.maxTargetDegree<<','
           <<sampling<<','<<recovery<<','<<prepare<<','<<boundary<<','<<validate<<','<<classify<<','<<pipeline<<','
           <<r.boundaryDistanceUpperBound<<','<<r.boundaryTangentUpperBoundDeg<<','<<r.outputMotherDistanceUpperBound<<','<<audit<<','<<Quote(r.message)<<','
           <<fused<<','<<positionOnly<<','<<fastChecks<<','<<fastAccepted<<','<<fallbacks<<','<<fixture.input->DegreeU()<<','<<fixture.input->DegreeV()<<','
           <<fixture.input->ControlPoints().size()*fixture.input->ControlPoints().front().size()<<','<<residualChecks<<','<<residualAccepted<<','<<residualFallbacks<<'\n';out.flush();
    }
}
}
int main(int argc,char **argv)
{
    sggk::init();try
    {
        Require(argc>=3,"Usage: --rational-unit data | --rational-case data csv id repeats synthetic --rational-v2-only");
        if(std::string(argv[1])=="--envelope-export"){Require(argc==3,"Invalid export arguments");EnvelopeExport(argv[2]);}
        else if(std::string(argv[1])=="--rational-unit"){Require(argc==3,"Invalid unit arguments");Unit(argv[2]);}
        else if(std::string(argv[1])=="--rational-full-case"||std::string(argv[1])=="--rational-fixed-case")
        {
            Require(argc==9&&std::string(argv[7])=="--rational-v2-only","Invalid full-suite arguments");
            FullCase(argv[2],argv[3],argv[4],unsigned(std::stoul(argv[5])),unsigned(std::stoul(argv[6])),std::string(argv[1])=="--rational-fixed-case",std::stod(argv[8]));
        }
        else
        {
            // 独立命令名与参数数量使废弃前的二进制在进入旧新对照循环之前即失败。
            Require(std::string(argv[1])=="--rational-case"&&argc==8&&std::string(argv[7])=="--rational-v2-only",
                "Invalid rational-only benchmark arguments");
            Case(argv[2],argv[3],unsigned(std::stoul(argv[4])),unsigned(std::stoul(argv[5])),std::stoi(argv[6])!=0);
        }
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
