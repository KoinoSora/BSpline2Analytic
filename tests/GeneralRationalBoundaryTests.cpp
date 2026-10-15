#include "GeneralRationalBoundaryInternal.h"
#include "GeneralDataset.h"
#include "GeneralExtremeDataset.h"
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
            if(!patch.tangentNotApplicable[e]&&i>0&&i<32)
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
    // Independent de Casteljau evaluation on exact dyadic input coefficients;
    // exercise unequal degrees and both constant/nonconstant weight paths.
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
    IV cross=Range(products.Cross(st,tt));
    Require(products.Dot(st,tt).range().lo>0&&Sqrt(Sqr(cross.x)+Sqr(cross.y)+Sqr(cross.z)).hi<1e-12,
            "Full tangent fallback should accept parallel speed mismatch");
    Require(!std::isfinite(ratio(s,line(-1,0,0))),"Opposite direction accepted by residual criterion");
    Require(!std::isfinite(ratio(line(0,0,0),s)),"Zero source speed accepted");
    // Same straight line in rational quadratic form, with W=1+t.
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
}
void Unit(const std::filesystem::path &data)
{
    ErrorKernelUnit();
    const std::array<double,4> domain{.2,.9,.6,1.2};
    for(unsigned type=0;type<5;++type)
    {
        auto a=Mother(type);auto s=Fixture(a,domain,true);GeneralRationalBoundaryOptions o;o.chartDomain=domain;
        auto r=ValidateGeneralRationalBoundaries(*s,a,o);
        Require(r.Passed(),"Rational skew fixture rejected type "+std::to_string(type)+": "+r.message);Audit(*s,r);
        Require(r.maxTargetDegree<=4&&!r.certifiesInteriorContinuum,"Wrong certificate contract");
        s->UReverse();auto reverse=ValidateGeneralRationalBoundaries(*s,a,o);
        Require(reverse.Passed(),"Reversed rational boundary rejected");Audit(*s,reverse);
        o.maxBoundaryNodes=1;auto budget=ValidateGeneralRationalBoundaries(*s,a,o);
        Require(budget.status==HybridValidationStatus::BudgetExceeded&&!budget.patch,"Budget leaked output");
    }
    auto plane=Mother(0);auto flat=Fixture(plane,domain,false);
    auto variableSpeed=Fixture(plane,domain,false);variableSpeed->UDegreeElevation(3-variableSpeed->DegreeU());
    for(unsigned j=0;j<variableSpeed->ControlPoints().front().size();++j)
        variableSpeed->SetControlPoint(Point(Get(variableSpeed->ControlPoints()[1][j])+Get(plane.xDirection)*1e-6),1,j);
    GeneralRationalBoundaryOptions speedOptions;speedOptions.chartDomain=domain;
    speedOptions.general.fitting.validation.boundaryTangentToleranceDeg=1e-6;
    auto speedReport=ValidateGeneralRationalBoundaries(*variableSpeed,plane,speedOptions);
    Require(speedReport.Passed()&&speedReport.tangentFallbacks>0,"Parallel variable-speed boundary must pass via full tangent fallback");
    GeneralRationalBoundaryOptions bad;bad.general.fitting.validation.distanceTolerance=-1;
    Require(ValidateGeneralRationalBoundaries(*flat,plane,bad).status==HybridValidationStatus::InvalidInput,"Invalid tolerance accepted");
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
    auto oldOnly=ReadGeneralFixture(data/GeneralFixtureFilename(254));GeneralRationalBoundaryOptions fixed;
    std::array<double,4> box{INFINITY,-INFINITY,INFINITY,-INFINITY};
    for(const auto &p:oldOnly.corners){box[0]=std::min(box[0],p.X());box[1]=std::max(box[1],p.X());box[2]=std::min(box[2],p.Y());box[3]=std::max(box[3],p.Y());}
    fixed.chartDomain=box;
    Require(!ValidateGeneralRationalBoundaries(*oldOnly.input,oldOnly.truth,fixed).Passed(),"Old analytic-UV contract substituted for rational-domain contract");
    for(unsigned id:{0u,218u,318u})
    {
        auto f=ReadGeneralFixture(data/GeneralFixtureFilename(id));auto r=ValidateGeneralRationalBoundaries(*f.input,f.truth);
        Require(r.Passed(),"Known rectangular/collapsed fixture failed "+std::to_string(id)+": "+r.message);Audit(*f.input,r);
        if(id)Require(r.patch->tangentNotApplicable[0]||r.patch->tangentNotApplicable[1]||
            r.patch->tangentNotApplicable[2]||r.patch->tangentNotApplicable[3],"Collapsed source tangent not marked");
    }
    // Piecewise chart composition across internal knots, checked against a test-
    // owned inverse circular parameter formula, including a periodic seam.
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
              const std::string &suite,unsigned id,unsigned repeats)
{
    bool synthetic=suite=="synthetic",extra=suite=="extra",extreme=suite=="extreme";
    Require(synthetic||extra||extreme||suite=="main","Unknown rational suite");
    Require(repeats>0&&id<(synthetic?5u:extra?60u:extreme?320u:600u),"Invalid full-suite range");
    GeneralFixture c;std::string group;FixtureCertificate cert;double certificateMs=0;
    const std::array<double,4> domain{.2,.9,.6,1.2};
    if(synthetic)
    {c.truth=Mother(id);c.type=id;c.input=Fixture(c.truth,domain,true,true);c.id=id;group="nurbs_chart_pressure";cert.upper=std::numeric_limits<double>::quiet_NaN();}
    else
    {
        c=ReadGeneralFixture(data/(extreme?GeneralExtremeFilename(id):GeneralFixtureFilename(id,extra)));
        Require(c.id==id,"Fixture identity mismatch");group=FixtureGroupName(c.group);
        std::cout<<"BEGIN_INPUT_CERT"<<std::endl;auto t=std::chrono::steady_clock::now();
        if(extreme)
        {
            auto spec=GeneralExtremeSpecs()[id];group=spec.group;
            cert=CertifyGeneralExtremePair(c,spec.sag>0?std::min(1e-6,spec.sag*.1):1e-6);
        }
        else cert=c.group==GeneralFixtureGroup::Infeasible?CertifyNoQuartic(*c.input,c.exclusionRadius):CertifyAnalyticPair(c,c.pairTolerance);
        certificateMs=Ms(t);Require(cert.passed,"Input certificate failed: "+cert.message);
    }
    GeneralRationalBoundaryOptions options;
    bool supplied=synthetic||(extra&&(id<20||id>=40));
    if(synthetic)options.chartDomain=domain;
    if(extra&&c.group==GeneralFixtureGroup::InvalidOptions)
    {if(id%2)options.general.fitting.validation.distanceTolerance=-1;else options.general.fitting.sampling.maxSamples=0;}
    if(extra&&c.group==GeneralFixtureGroup::ValidationBudget)options.maxBoundaryNodes=1;
    if(extra&&c.group==GeneralFixtureGroup::SpherePriority)options.general.similarScoreAbsoluteRatio=1e6;
    const std::string header="suite,id,run,mode,group,source_type,rational,degree_u,degree_v,controls,source_certified,source_bound,certificate_ms,entry,status,accepted,audit_status,audit_sample_max,classify_ms,validate_ms,sampling_ms,recovery_ms,prepare_ms,boundary_ms,total_ms,attempts,chart_attempts,nodes,refinements,fast_checks,fast_accepted,tangent_fallbacks,position_only,fused_errors,boundary_bound,tangent_deg,mother_bound,output_type,raw_candidates,full_scores,message,audit_message\n";
    std::vector<std::string> rows;bool anyFailure=false;
    auto save=[&]{std::ofstream out(csv);Require(bool(out),"Cannot write full-suite CSV");out<<header;for(const auto &row:rows)out<<row;};
    auto attemptPath=csv.parent_path()/(csv.stem().string()+"-attempts.csv");
    std::ofstream attemptsFile(attemptPath);Require(bool(attemptsFile),"Cannot write attempt diagnostics");
    attemptsFile<<"suite,id,run,attempt,type,status,accepted,charts,nodes,refinements,sampling_ms,recovery_ms,prepare_ms,boundary_ms,total_ms,message\n"<<std::setprecision(17);
    struct Attempt{sggk::SurfaceType type;GeneralRationalBoundaryReport report;};
    for(unsigned run=0;run<repeats;++run)
    {
        std::cout<<"BEGIN_ALGORITHM "<<run<<std::endl;
        auto start=std::chrono::steady_clock::now();GeneralRationalBoundaryReport r;
        GeneralClassificationDiagnostics classification;std::vector<Attempt> attempts;
        double classify=0,validate=0,sampling=0,recovery=0,prepare=0,boundary=0;
        unsigned charts=0,nodes=0,refinements=0,fast=0,fastAccepted=0,fallback=0,position=0,fused=0;
        try
        {
            SurfaceSampleContext samples(*c.input,options.general.fitting.sampling.maxSamples);
            auto attempt=[&](const AnalyticParameters &a)
            {
                r=rb::Validate(samples,a,options);validate+=r.totalMs;sampling+=r.samplingMs;recovery+=r.recoveryMs;
                prepare+=r.preparationMs;boundary+=r.boundaryMs;charts+=r.chartAttempts;nodes+=r.boundaryNodes;
                refinements+=r.mappingRefinements;fast+=r.tangentFastChecks;fastAccepted+=r.tangentFastAccepted;
                fallback+=r.tangentFallbacks;position+=r.positionOnlyNodes;fused+=r.fusedErrorBuilds;
                attempts.push_back({a.type,r});
            };
            if(supplied)attempt(c.truth);
            else
            {
                auto t=std::chrono::steady_clock::now();auto pool=InitializeGeneralCandidates(samples,options.general,classification);classify+=Ms(t);
                for(;;)
                {
                    t=std::chrono::steady_clock::now();auto batch=NextGeneralCandidateBatch(pool,samples,options.general,classification);classify+=Ms(t);
                    if(batch.empty())break;
                    for(const auto &candidate:batch){attempt(candidate.parameters);if(r.Passed())break;}
                    if(r.Passed())break;
                }
                if(attempts.empty())r.message="No classified candidate available";
            }
        }
        catch(const std::exception &e){r.patch.reset();r.status=HybridValidationStatus::NotCertified;r.message=std::string("Algorithm exception: ")+e.what();}
        double total=Ms(start);bool passed=r.Passed();anyFailure|=!passed;
        std::cout<<"END_ALGORITHM "<<run<<std::endl;
        std::string auditStatus=passed?"pending":"not_applicable",auditMessage;double audit=std::numeric_limits<double>::quiet_NaN();
        auto line=[&]
        {
            std::ostringstream out;out<<std::setprecision(17)<<suite<<','<<id<<','<<run<<",rational,"<<group<<','<<(c.type<5?TypeName(c.truth.type):"Freeform")<<','<<c.input->IsRational()<<','
                <<c.input->DegreeU()<<','<<c.input->DegreeV()<<','<<c.input->ControlPoints().size()*c.input->ControlPoints().front().size()<<','<<(synthetic?-1:int(cert.passed))<<','<<cert.upper<<','<<certificateMs<<','
                <<(synthetic?"fixed_chart":supplied?"supplied_mother":"classified")<<','<<HybridValidationStatusName(r.status)<<','<<passed<<','<<auditStatus<<','<<audit<<','
                <<classify<<','<<validate<<','<<sampling<<','<<recovery<<','<<prepare<<','<<boundary<<','<<total<<','<<attempts.size()<<','<<charts<<','<<nodes<<','<<refinements<<','
                <<fast<<','<<fastAccepted<<','<<fallback<<','<<position<<','<<fused<<','<<r.boundaryDistanceUpperBound<<','<<r.boundaryTangentUpperBoundDeg<<','<<r.outputMotherDistanceUpperBound<<','
                <<(passed?TypeName(r.patch->mother.type):"Unknown")<<','<<classification.rawCandidateCount<<','<<classification.fullScoreCount<<','<<Quote(r.message)<<','<<Quote(auditMessage)<<'\n';
            return out.str();
        };
        // Save measured algorithm time before the independent SDK audit, so
        // an audit watchdog cannot erase an already completed algorithm result.
        rows.push_back(line());save();
        for(size_t i=0;i<attempts.size();++i)
        {
            const auto &a=attempts[i];const auto &v=a.report;
            attemptsFile<<suite<<','<<id<<','<<run<<','<<i<<','<<TypeName(a.type)<<','<<HybridValidationStatusName(v.status)<<','<<v.Passed()<<','<<v.chartAttempts<<','<<v.boundaryNodes<<','<<v.mappingRefinements<<','
                <<v.samplingMs<<','<<v.recoveryMs<<','<<v.preparationMs<<','<<v.boundaryMs<<','<<v.totalMs<<','<<Quote(v.message)<<'\n';
        }
        attemptsFile.flush();
        if(passed)
        {
            std::cout<<"BEGIN_AUDIT "<<run<<std::endl;
            try{audit=Audit(*c.input,r,1e-6);auditStatus="passed";}
            catch(const std::exception &e){auditStatus="unresolved";auditMessage=e.what();}
            rows.back()=line();save();std::cout<<"END_AUDIT "<<run<<std::endl;
        }
    }
    // Diagnostic oracle only after all measured runs; never feeds the candidate pool.
    if(anyFailure&&((suite=="main"&&id<500)||extreme))
    {
        std::cout<<"BEGIN_KNOWN_MOTHER_DIAGNOSTIC"<<std::endl;
        auto known=ValidateGeneralRationalBoundaries(*c.input,c.truth);
        std::ofstream out(csv.parent_path()/(csv.stem().string()+"-known.csv"));
        out<<"suite,id,status,accepted,total_ms,charts,nodes,message\n"<<std::setprecision(17)<<suite<<','<<id<<','<<HybridValidationStatusName(known.status)<<','<<known.Passed()<<','
            <<known.totalMs<<','<<known.chartAttempts<<','<<known.boundaryNodes<<','<<Quote(known.message)<<'\n';
        std::cout<<"END_KNOWN_MOTHER_DIAGNOSTIC"<<std::endl;
    }
    if(extra&&id==49)
    {
        // This historical probe demonstrates the sampled-interior contract gap.
        SurfaceSampleContext samples(*c.input);double sampled=0,peak=0;
        for(const auto &p:samples.Greville(false))sampled=std::max(sampled,MotherDistance(c.truth,p.point));
        for(unsigned i=0;i<=10000;++i)peak=std::max(peak,MotherDistance(c.truth,Get(c.input->CalcPoint(double(i)/10000,.5))));
        std::ofstream out(csv.parent_path()/(csv.stem().string()+"-interior.csv"));out<<"greville_max,dense_witness_max\n"<<std::setprecision(17)<<sampled<<','<<peak<<'\n';
    }
}
void Case(const std::filesystem::path &data,const std::filesystem::path &csv,unsigned id,unsigned repeats,bool synthetic)
{
    auto fixture=synthetic?GeneralFixture{}:ReadGeneralFixture(data/GeneralFixtureFilename(id));
    std::array<double,4> domain{.2,.9,.6,1.2};
    if(synthetic){fixture.truth=Mother(id);fixture.input=Fixture(fixture.truth,domain,true,true);fixture.id=id;}
    else
    {auto cert=fixture.group==GeneralFixtureGroup::Infeasible?CertifyNoQuartic(*fixture.input,fixture.exclusionRadius):CertifyAnalyticPair(fixture,fixture.pairTolerance);Require(cert.passed,"Legacy fixture certificate failed");}
    std::ofstream out(csv);out<<std::setprecision(17);
    out<<"id,synthetic,mode,run,status,accepted,attempts,chart_attempts,nodes,refinements,max_degree,sampling_ms,recovery_ms,prepare_ms,boundary_ms,validate_ms,classify_ms,pipeline_ms,boundary_bound,tangent_deg,mother_bound,audit_sample_max,message,fused_errors,position_only,fast_checks,fast_accepted,tangent_fallbacks,source_degree_u,source_degree_v,source_controls\n";
    for(unsigned run=0;run<repeats;++run)
    {
        auto start=std::chrono::steady_clock::now();GeneralRationalBoundaryOptions options;
        if(synthetic)options.chartDomain=domain;
        SurfaceSampleContext samples(*fixture.input,options.general.fitting.sampling.maxSamples);
        GeneralRationalBoundaryReport r;unsigned attempts=0,charts=0,nodes=0,refinements=0;
        unsigned fused=0,positionOnly=0,fastChecks=0,fastAccepted=0,fallbacks=0;
        double classify=0,validate=0,sampling=0,recovery=0,prepare=0,boundary=0;bool passed=false;
        auto attempt=[&](const AnalyticParameters &a)
        {
            ++attempts;
            r=rb::Validate(samples,a,options);passed=r.Passed();validate+=r.totalMs;sampling+=r.samplingMs;
            recovery+=r.recoveryMs;prepare+=r.preparationMs;boundary+=r.boundaryMs;
            charts+=r.chartAttempts;nodes+=r.boundaryNodes;refinements+=r.mappingRefinements;
            fused+=r.fusedErrorBuilds;positionOnly+=r.positionOnlyNodes;fastChecks+=r.tangentFastChecks;
            fastAccepted+=r.tangentFastAccepted;fallbacks+=r.tangentFallbacks;
        };
        if(synthetic)attempt(fixture.truth); // Explicit-chart validator benchmark, not classifier timing.
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
           <<fixture.input->ControlPoints().size()*fixture.input->ControlPoints().front().size()<<'\n';out.flush();
    }
}
}
int main(int argc,char **argv)
{
    sggk::init();try
    {
        Require(argc>=3,"Usage: --rational-unit data | --rational-case data csv id repeats synthetic --rational-only");
        if(std::string(argv[1])=="--rational-unit"){Require(argc==3,"Invalid unit arguments");Unit(argv[2]);}
        else if(std::string(argv[1])=="--rational-full-case")
        {
            Require(argc==8&&std::string(argv[7])=="--rational-only","Invalid full-suite arguments");
            FullCase(argv[2],argv[3],argv[4],unsigned(std::stoul(argv[5])),unsigned(std::stoul(argv[6])));
        }
        else
        {
            // The distinct command and arity make pre-retirement binaries fail
            // before entering their old/new comparison loop.
            Require(std::string(argv[1])=="--rational-case"&&argc==8&&std::string(argv[7])=="--rational-only",
                "Invalid rational-only benchmark arguments");
            Case(argv[2],argv[3],unsigned(std::stoul(argv[4])),unsigned(std::stoul(argv[5])),std::stoi(argv[6])!=0);
        }
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
