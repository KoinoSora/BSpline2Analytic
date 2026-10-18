/**
 * @file GeneralRationalRecoveryDiagnostics.cpp
 * @brief 当前有理边界恢复失败的定向诊断。
 * 区分数值搜索上界见证与严格排除证书；只对固定母面/参数图得出的结论不能扩大到所有可能解。
 */
#include "GeneralRationalBoundaryInternal.h"
#include "GeneralDataset.h"
#include "GeneralRationalDataset.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace rb=bspline2analytic::detail::rational_boundary;
namespace
{
void Require(bool x,const char *s){if(!x)throw std::runtime_error(s);}
std::string Quote(std::string s){std::string o="\"";for(char c:s){if(c=='"')o+='"';o+=c=='\n'?' ':c;}return o+'"';}
I Length(IV a){return Sqrt(Sqr(a.x)+Sqr(a.y)+Sqr(a.z));}
I Eval(const BPoly &p,I t)
{
    Require(p.v==0,"Expected one-dimensional polynomial");auto c=p.c;
    for(unsigned n=p.u;n>0;--n)for(unsigned i=0;i<n;++i)c[i]=c[i]*(I(1)-t)+c[i+1]*t;
    return c[0];
}
IV PointAt(const BezierCell &c,I t)
{I w=Eval(c.w,t);Require(w.lo>0,"Point weight unresolved");return {Eval(c.p.x,t)/w,Eval(c.p.y,t)/w,Eval(c.p.z,t)/w};}
IV SourcePoint(const std::vector<BezierCell> &edge,double t)
{
    for(const auto &c:edge)if(c.u0<=t&&t<=c.u1)return PointAt(c,(I(t)-I(c.u0))/(I(c.u1)-I(c.u0)));
    throw std::runtime_error("Source witness outside extracted edge");
}
double HullLower(const BezierCell &c,IV point)
{
    Require(c.p.x.u==c.w.u&&c.p.y.u==c.w.u&&c.p.z.u==c.w.u,"Unequal hull degrees");
    IV box;bool first=true;
    auto add=[&](I &a,I x){if(first)a=x;else{a.lo=std::min(a.lo,x.lo);a.hi=std::max(a.hi,x.hi);}};
    for(unsigned i=0;i<=c.w.u;++i)
    {
        Require(c.w.c[i].lo>0,"Nonpositive hull weight");
        add(box.x,c.p.x.c[i]/c.w.c[i]);add(box.y,c.p.y.c[i]/c.w.c[i]);add(box.z,c.p.z.c[i]/c.w.c[i]);first=false;
    }
    return Length(point-box).lo;
}
struct Exclusion {bool proved=false;double lower=0;unsigned nodes=0;};
Exclusion Exclude(const std::vector<BezierCell> &edge,IV point,double tolerance)
{
    struct Node{BezierCell cell;unsigned depth;};std::vector<Node> stack;
    for(const auto &c:edge)stack.push_back({c,0});Exclusion out;double lower=INFINITY;
    while(!stack.empty())
    {
        if(++out.nodes>4096)return out;
        auto n=std::move(stack.back());stack.pop_back();double d=HullLower(n.cell,point);
        if(d>tolerance){lower=std::min(lower,d);continue;}
        if(n.depth>=40)return out;
        auto halves=Split(n.cell,true);
        stack.push_back({std::move(halves.first),n.depth+1});stack.push_back({std::move(halves.second),n.depth+1});
    }
    out.proved=true;out.lower=lower;return out;
}
// 数值搜索只提供距离上界见证，不构成排除证书。
std::pair<double,double> Nearby(const std::vector<BezierCell> &edge,IV point)
{
    double best=INFINITY,parameter=0;
    for(const auto &c:edge)
    {
        auto distance=[&](double t){return Length(point-PointAt(c,I(t))).hi;};
        double values[65]{};
        for(unsigned i=0;i<=64;++i)
        {values[i]=distance(double(i)/64);if(values[i]<best){best=values[i];parameter=c.u0+(c.u1-c.u0)*i/64;}}
        for(unsigned i=1;i<64;++i)if(values[i]<=values[i-1]&&values[i]<=values[i+1])
        {
            double lo=double(i-1)/64,hi=double(i+1)/64;
            for(unsigned k=0;k<48;++k)
            {double x=lo+(hi-lo)/3,y=hi-(hi-lo)/3;if(distance(x)<distance(y))hi=y;else lo=x;}
            double t=(lo+hi)*.5,d=distance(t);if(d<best){best=d;parameter=c.u0+(c.u1-c.u0)*t;}
        }
    }
    return {best,parameter};
}
void KernelUnit()
{
    BezierCell c{{BPoly(1,0),BPoly(1,0),BPoly(1,0)},0,1,0,1,BPoly(1,0,I(1))};c.p.x.c={I(0),I(1)};
    auto far=Exclude({c},Exact({.5,.01,0}),1e-5);Require(far.proved&&far.lower>1e-5,"Far line witness not excluded");
    Require(!Exclude({c},Exact({.5,1e-6,0}),1e-5).proved,"Near line falsely excluded");
    auto near=Nearby({c},Exact({.5,1e-6,0}));Require(near.first>=1e-6&&near.first<1.001e-6,"Nearest upper witness incorrect");
    c.w.c={I(1),I(2)};c.p.x.c={I(0),I(2)};
    Require(Exclude({c},Exact({.5,.01,0}),1e-5).proved,"Rational line witness not excluded");
    double w=std::sqrt(.5);BezierCell arc{{BPoly(2,0),BPoly(2,0),BPoly(2,0)},0,1,0,1,BPoly(2,0)};
    arc.w.c={I(1),I(w),I(1)};arc.p.x.c={I(1),I(w),I(0)};arc.p.y.c={I(0),I(w),I(1)};
    auto inside=Exclude({arc},Exact({.5,.5,0}),1e-5);
    Require(inside.proved&&inside.nodes>1,"Curved hull subdivision not exercised");
    Require(!Exclude({arc},Exact({1,0,0}),1e-5).proved,"Curve endpoint falsely excluded");
    std::cout<<"Diagnostic exclusion kernel checks passed\n";
}
struct Diagnostic
{
    GeneralFixture fixture;GeneralRationalBoundaryOptions options;SurfaceSampleContext samples;
    std::array<std::vector<BezierCell>,4> source;V sourceOrigin;
    std::ofstream attempts,charts,witnesses;
    Diagnostic(const std::filesystem::path &data,const std::filesystem::path &prefix,unsigned id):
        fixture(ReadGeneralFixture(data/"main"/RationalFixtureName("main",id))),samples(*fixture.input),sourceOrigin(Get(fixture.truth.origin)),
        attempts(prefix.string()+"-candidates.csv"),charts(prefix.string()+"-charts.csv"),witnesses(prefix.string()+"-witnesses.csv")
    {
        Require(fixture.id==id,"Fixture identity mismatch");Require(attempts&&charts&&witnesses,"Cannot create diagnostics");
        auto certificate=CertifyRationalDomain(fixture);Require(certificate.passed,"Source certificate failed");
        source=ExtractHybridBoundaries(*fixture.input,sourceOrigin);
        attempts<<"id,mode,candidate,type,status,accepted,charts,nodes,interior_max,finite_max,boundary_bound,tangent_deg,center_delta,radius_delta,axis_abs_dot,origin_x,origin_y,origin_z,axis_x,axis_y,axis_z,xdir_x,xdir_y,xdir_z,radius,major_radius,half_angle,message\n"<<std::setprecision(17);
        charts<<"id,mode,candidate,type,chart,recovered,reason,edge,sample,source_fraction,value,limit,u0,u1,v0,v1,axis_x,axis_y,axis_z,xdir_x,xdir_y,xdir_z,exception\n"<<std::setprecision(17);
        witnesses<<"id,mode,candidate,type,chart,edge,source_fraction,uv_projection_distance,nearby_upper,nearby_fraction,outside_certified,global_lower,proof_tolerance,proof_nodes,source_box_width\n"<<std::setprecision(17);
    }
    void Probe(const std::string &mode,unsigned candidate,unsigned chart,const rb::Trim &trim,const rb::RecoveryDiagnostics &d)
    {
        struct Choice{unsigned edge;double fraction,distance;};std::vector<Choice> choices;
        for(unsigned e=0;e<4;++e)
        {
            auto a=trim.polygon.corners[rb::Ends[e][0]],b=trim.polygon.corners[rb::Ends[e][1]];
            double dx=b.X()-a.X(),dy=b.Y()-a.Y(),length=dx*dx+dy*dy;if(!(length>0))continue;
            Choice choice{e,0,-1};
            for(const auto &p:d.walks[e])
            {
                double f=std::clamp(((p.second[0]-a.X())*dx+(p.second[1]-a.Y())*dy)/length,0.,1.);
                double distance=Norm(rb::ChartPoint(trim,p.second)-rb::ChartPoint(trim,{a.X()+f*dx,a.Y()+f*dy}));
                if(distance>choice.distance)choice={e,p.first,distance};
            }
            choices.push_back(choice);
        }
        std::sort(choices.begin(),choices.end(),[](const Choice &a,const Choice &b){return a.distance>b.distance;});
        paired::Products products;
        for(const auto &choice:choices)
        {
            unsigned e=choice.edge;const auto &s=*fixture.input;
            double lo=e<2?s.MinParamV():s.MinParamU(),hi=e<2?s.MaxParamV():s.MaxParamU();
            IV p=SourcePoint(source[e],lo+(hi-lo)*choice.fraction)+(Exact(sourceOrigin)-Exact(Get(trim.mother.origin)));
            auto edge=rb::ComposeEdge(trim,e,products);
            for(auto &c:edge)
            {
                unsigned degree=std::max({c.p.x.u,c.p.y.u,c.p.z.u,c.w.u});
                c.p={products.Elevate(c.p.x,degree,0),products.Elevate(c.p.y,degree,0),products.Elevate(c.p.z,degree,0)};
                c.w=products.Elevate(c.w,degree,0);
            }
            auto near=Nearby(edge,p);Exclusion proof;
            double tolerance=options.general.fitting.validation.distanceTolerance;
            // 计入 BuildTarget 使用的完整“输出到理想曲线”误差预算，使排除结论对预算内的输出扰动仍成立。
            double proofTolerance=(I(tolerance)+I(tolerance*.1)).hi;
            if(near.first>=proofTolerance)proof=Exclude(edge,p,proofTolerance);
            double width=std::max({p.x.hi-p.x.lo,p.y.hi-p.y.lo,p.z.hi-p.z.lo});
            witnesses<<fixture.id<<','<<mode<<','<<candidate<<','<<TypeName(trim.mother.type)<<','<<chart<<','<<e<<','<<choice.fraction<<','
                <<choice.distance<<','<<near.first<<','<<near.second<<','<<proof.proved<<','<<proof.lower<<','<<proofTolerance<<','<<proof.nodes<<','<<width<<'\n';
            // 只需严格排除一个源点，即可否定固定参数图的理想直线边在任意单调对应下满足约束。
            if(proof.proved)break;
        }
    }
    GeneralRationalBoundaryReport Analyze(const std::string &mode,unsigned candidate,const AnalyticParameters &a)
    {
        auto r=rb::Validate(samples,a,options);auto origin=Get(a.origin),axis=Get(a.axis),x=Get(a.xDirection);
        if(r.Passed())
        {
            Require(r.boundaryDistanceUpperBound<1e-5&&r.boundaryTangentUpperBoundDeg<.1&&r.outputMotherDistanceUpperBound<1e-6,"Invalid accepted bounds");
            const auto &s=*fixture.input;
            for(unsigned e=0;e<4;++e)for(unsigned j=0;j<=16;++j)
            {
                double f=double(j)/16,u=e<2?(e==0?s.MinParamU():s.MaxParamU()):s.MinParamU()+(s.MaxParamU()-s.MinParamU())*f;
                double v=e<2?s.MinParamV()+(s.MaxParamV()-s.MinParamV())*f:(e==2?s.MinParamV():s.MaxParamV());
                auto sourcePoint=samples.At(u,v);const auto &curve=*r.patch->boundaries[e];double parameter=f;
                V nearest=Get(curve.CalcNearestPoint(Point(sourcePoint.point),parameter));
                Require(Norm(sourcePoint.point-nearest)<1.001e-5,"Counterfactual output distance audit failed");
                if(j&&j<16&&!r.patch->tangentNotApplicable[e])
                {
                    V st=e<2?sourcePoint.dv:sourcePoint.du,tt=Get(curve.CalcDeriv1(parameter));
                    Require(Dot(st,tt)>0&&Norm(Cross(st,tt))/Dot(st,tt)<std::tan(.1*Pi/180),"Counterfactual tangent audit failed");
                }
            }
        }
        attempts<<fixture.id<<','<<mode<<','<<candidate<<','<<TypeName(a.type)<<','<<HybridValidationStatusName(r.status)<<','<<r.Passed()<<','<<r.chartAttempts<<','<<r.boundaryNodes<<','
            <<r.interiorSampleMaxDistance<<','<<r.interiorSampleMaxFiniteDistance<<','<<r.boundaryDistanceUpperBound<<','<<r.boundaryTangentUpperBoundDeg<<','<<Norm(origin-Get(fixture.truth.origin))<<','<<std::abs(a.radius-fixture.truth.radius)<<','
            <<std::abs(Dot(axis,Get(fixture.truth.axis)))<<','<<origin.x<<','<<origin.y<<','<<origin.z<<','<<axis.x<<','<<axis.y<<','<<axis.z<<','<<x.x<<','<<x.y<<','<<x.z<<','
            <<a.radius<<','<<a.majorRadius<<','<<a.halfAngle<<','<<Quote(r.message)<<'\n';
        if(!r.chartAttempts)return r;
        auto proposals=rb::RecoveryCharts(samples,a,options);Require(r.chartAttempts<=proposals.size(),"Unexpected chart count");
        for(unsigned i=0;i<r.chartAttempts;++i)
        {
            rb::Trim trim;rb::RecoveryDiagnostics d;bool recovered=false;std::string exception;
            try{recovered=rb::Recover(samples,proposals[i],options,trim,&d);}catch(const std::exception &e){exception=e.what();}
            rb::Trim plain;bool plainRecovered=false;std::string plainException;
            try{plainRecovered=rb::Recover(samples,proposals[i],options,plain);}catch(const std::exception &e){plainException=e.what();}
            Require(recovered==plainRecovered&&exception==plainException,"Diagnostic sink changed recovery decision");
            auto z=Get(proposals[i].axis),xx=Get(proposals[i].xDirection);
            charts<<fixture.id<<','<<mode<<','<<candidate<<','<<TypeName(a.type)<<','<<i<<','<<recovered<<','<<d.reason<<','<<d.edge<<','<<d.sample<<','<<d.parameter<<','<<d.value<<','<<d.limit;
            for(double v:trim.domain)charts<<','<<v;
            charts<<','<<z.x<<','<<z.y<<','<<z.z<<','<<xx.x<<','<<xx.y<<','<<xx.z<<','<<Quote(exception)<<'\n';
            if(d.reason=="RationalLineWitnessExceeded")Probe(mode,candidate,i,trim,d);
        }
        attempts.flush();charts.flush();witnesses.flush();return r;
    }
    void Run()
    {
        GeneralClassificationDiagnostics classification;auto pool=InitializeGeneralCandidates(samples,options.general,classification);
        unsigned candidate=0;std::optional<AnalyticParameters> firstMatching;
        bool passed=false;
        for(;;)
        {
            auto batch=NextGeneralCandidateBatch(pool,samples,options.general,classification);if(batch.empty())break;
            for(const auto &seed:batch)
            {
                auto r=Analyze("classified",candidate++,seed.parameters);
                if(!firstMatching&&seed.parameters.type==fixture.truth.type&&r.chartAttempts)firstMatching=seed.parameters;
                if(r.Passed()){passed=true;break;}
            }
            if(passed)break;
        }
        options.chartDomain=fixture.chartDomain;
        auto known=Analyze("known_mother",0,fixture.truth); // V2 数据提供的参考参数图只用于计时之外的诊断。
        options.chartDomain.reset();
        if(firstMatching&&fixture.truth.type==sggk::SurfaceType::Sphere)
        {
            auto a=*firstMatching;a.axis=fixture.truth.axis;a.xDirection=fixture.truth.xDirection;
            Analyze("candidate_with_truth_frame",0,a);
            a=fixture.truth;a.axis=firstMatching->axis;a.xDirection=firstMatching->xDirection;
            Analyze("truth_with_candidate_frame",0,a);
        }
        std::cout<<"COMPLETE "<<fixture.id<<" classified="<<passed<<" known="<<known.Passed()<<" attempts="<<candidate<<'\n';
    }
};
}
int main(int argc,char **argv)
{
    sggk::init();try
    {
        Require(IntervalEnvironmentSupported(),"Unsupported interval environment");
        if(argc==2&&std::string(argv[1])=="--diagnostic-unit"){KernelUnit();return 0;}
        Require(argc==6&&std::string(argv[1])=="--diagnose-main"&&std::string(argv[5])=="--rational-v2-only","Invalid rational diagnostic arguments");
        Diagnostic run(argv[2],argv[3],unsigned(std::stoul(argv[4])));run.Run();return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
