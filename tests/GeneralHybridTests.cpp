/**
 * @file GeneralHybridTests.cpp
 * @brief 已废弃 ValidateGeneralHybrid 的历史单测。
 * 保留采样盲区等合同反例，不运行或用作新验证器的性能对照。
 */
#include "GeneralHybridInternal.h"
#include "GeneralDataset.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>

using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace
{
unsigned assertions=0;
void Require(bool value,const std::string &message){++assertions;if(!value)throw std::runtime_error(message);}
double Ms(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();}
std::string Quote(std::string s){std::string out="\"";for(char c:s){if(c=='"')out+='"';out+=c=='\n'||c=='\r'?' ':c;}return out+'"';}
void CheckReport(const GeneralHybridReport &r)
{
    Require(r.totalMs>=0&&std::isfinite(r.totalMs),"Invalid timing");
    Require(r.samplingMs+r.preparationMs+r.recoveryMs+r.boundaryMs<=r.totalMs+1e-6,"Stage timing double count");
    if(r.Passed())Require(r.interiorSampleCount&&std::isfinite(r.boundaryDistanceUpperBound)&&
        std::isfinite(r.boundaryTangentUpperBoundDeg),"Incomplete hybrid certificate");
    else Require(!r.patch,"Failure leaked a patch");
}
void CheckExtraction(const sggk::BSplineSurface &s,V origin)
{
    auto edges=ExtractHybridBoundaries(s,origin);
    for(unsigned e=0;e<4;++e)for(const auto &c:edges[e])for(double t:{0.,.27,.61,1.})
    {
        IV p=Range(Evaluate(c.p,t,0));I w=Evaluate(c.w,t,0).range();
        double f=c.u0+(c.u1-c.u0)*t;
        V expected=Get(s.CalcPoint(e<2?(e==0?s.MinParamU():s.MaxParamU()):f,
                                  e<2?f:(e==2?s.MinParamV():s.MaxParamV())))-origin;
        I x=p.x/w,y=p.y/w,z=p.z/w;
        // SDK 点值计算会独立舍入，比对时只给少量 ulp 余量。
        double slack=1e-10*std::max(1.,Norm(expected));
        Require(x.lo-slack<=expected.x&&expected.x<=x.hi+slack&&y.lo-slack<=expected.y&&expected.y<=y.hi+slack&&
                z.lo-slack<=expected.z&&expected.z<=z.hi+slack,"1D homogeneous extraction disagrees with SDK");
    }
}
sggk::BSplineSurfacePtr PolynomialPatch(double boundaryBump,double interiorBump)
{
    // 构造 z=boundaryBump*v*(1-v)+interiorBump*u*(1-u)*v*(1-v)，再升阶到双三次。
    sggk::Point3DMatrix cp(4,std::vector<sggk::Point3D>(4));
    const double b[]={0,1./3,1./3,0};
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)
        cp[i][j]=sggk::Point3D(i/3.,j/3.,boundaryBump*b[j]+interiorBump*b[i]*b[j]);
    return std::make_shared<sggk::BSplineSurface>(3,3,cp,sggk::RealArray{0,1},sggk::RealArray{0,1},
                                                sggk::UIntArray{4,4},sggk::UIntArray{4,4});
}
AnalyticParameters Plane()
{AnalyticParameters a;a.type=sggk::SurfaceType::Plane;a.origin={0,0,0};a.axis=sggk::Dir3D(0,0,1);a.xDirection=sggk::Dir3D(1,0,0);return a;}
void Unit(const std::filesystem::path &dir)
{
    for(unsigned id:{0u,1u,36u,39u,59u,77u,100u,101u,200u,201u,300u,301u,400u,401u})
    {
        auto c=ReadGeneralFixture(dir/GeneralFixtureFilename(id));
        CheckExtraction(*c.input,Get(c.truth.origin));
        auto r=ValidateGeneralHybrid(*c.input,c.truth);CheckReport(r);
        Require(r.Passed(),"Known feasible fixture failed: "+std::to_string(id)+" "+r.message);
        c.input->UReverse();CheckExtraction(*c.input,Get(c.truth.origin));
        auto reversed=ValidateGeneralHybrid(*c.input,c.truth);CheckReport(reversed);
        Require(reversed.Passed(),"Reversed fixture failed");
    }
    auto flat=PolynomialPatch(0,0);auto a=Plane();
    Require(ValidateGeneralHybrid(*flat,a).Passed(),"Exact plane rejected");
    auto raised=a;raised.origin={0,0,.001};
    Require(!ValidateGeneralHybrid(*flat,raised).Passed(),"Wrong mother accepted");
    GeneralHybridOptions options;options.maxBoundaryPieces=1;
    auto budget=ValidateGeneralHybrid(*flat,a,options);CheckReport(budget);
    Require(budget.status==HybridValidationStatus::BudgetExceeded,"Boundary budget ignored");
    options={};options.general.fitting.validation.distanceTolerance=-1;
    Require(ValidateGeneralHybrid(*flat,a,options).status==HybridValidationStatus::InvalidInput,"Invalid tolerance accepted");
    // 增强 Greville 采样包含中心，可捕获该经典纯 Greville 漏检反例；但这仍不构成连续域认证。
    Require(!ValidateGeneralHybrid(*PolynomialPatch(0,19e-5),a).Passed(),"Interior peak missed by augmented grid");
    // 位置误差很小，但切向误差远大于显式设置的严格角容差。
    options={};options.general.fitting.validation.boundaryTangentToleranceDeg=.00001;
    Require(!ValidateGeneralHybrid(*PolynomialPatch(1e-6,0),a,options).Passed(),"Tangent violation accepted");
    auto warped=PolynomialPatch(0,0);
    const double x[]={0,1./15,.4,1}; // 构造严格递增的 x(u)=0.2u+0.8u²。
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)warped->SetControlPoint({x[i],j/3.,0},i,j);
    auto mapped=ValidateGeneralHybrid(*warped,a);CheckReport(mapped);
    Require(mapped.Passed()&&mapped.mappingRefinements>0,"Nonlinear boundary parameter correspondence failed");
    // 刻意接受的采样盲区用于说明历史 API 合同：全部增强 Greville 点通过，而实际内部峰值超差。
    auto blind=PolynomialPatch(0,0);const double h[]={0,1./15,.4,0},b[]={0,1./3,1./3,0};
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)
        blind->SetControlPoint({i/3.,j/3.,20.75e-5*h[i]*b[j]},i,j);
    auto sampled=ValidateGeneralHybrid(*blind,a);CheckReport(sampled);
    double peak=(1.6+std::sqrt(4.96))/6;
    Require(sampled.Passed()&&!sampled.certifiesInteriorContinuum&&
        std::abs(blind->CalcPoint(peak,.5).Z())>1e-5,"Sampled evidence silently became a continuum claim");
    auto outside=PolynomialPatch(0,0);
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)
        outside->SetControlPoint({i/3.+10*b[i]*b[j],j/3.,0},i,j);
    Require(!ValidateGeneralHybrid(*outside,a).Passed(),"Mother membership substituted for finite patch membership");
    std::cout<<"Hybrid unit assertions="<<assertions<<'\n';
}
int Case(const std::filesystem::path &dir,const std::filesystem::path &csv,unsigned id,bool extra,unsigned repeats)
{
    auto fixture=ReadGeneralFixture(dir/GeneralFixtureFilename(id,extra));
    Require(fixture.id==id&&fixture.extra==extra,"Fixture identity mismatch");
    std::cout<<"BEGIN_INPUT_CERT "<<id<<std::endl;
    auto input=fixture.group==GeneralFixtureGroup::Infeasible?CertifyNoQuartic(*fixture.input,fixture.exclusionRadius):
        CertifyAnalyticPair(fixture,fixture.pairTolerance);
    Require(input.passed,"Input fixture certificate failed: "+input.message);
    std::ofstream out(csv);Require(bool(out),"Cannot write result CSV");out<<std::setprecision(17);
    out<<"id,extra,group,type,rational,controls,run,fixture_certified,status,accepted,output_type,output_audit,output_bound,attempts,samples,boundary_nodes,boundary_segments,mapping_refinements,collapsed_edges,chart_attempts,chart_iterations,classify_ms,sampling_ms,recovery_ms,prepare_ms,boundary_ms,validate_ms,pipeline_ms,interior_max,finite_sample_max,boundary_bound,tangent_deg,message\n";
    for(unsigned run=0;run<repeats;++run)
    {
        GeneralHybridOptions options;
        if(fixture.group==GeneralFixtureGroup::InvalidOptions)options.general.fitting.validation.distanceTolerance=-1;
        if(fixture.group==GeneralFixtureGroup::ValidationBudget)options.maxBoundaryPieces=1;
        std::cout<<"BEGIN_PIPELINE "<<run<<std::endl;
        auto start=std::chrono::steady_clock::now();SurfaceSampleContext samples(*fixture.input,options.general.fitting.sampling.maxSamples);
        GeneralClassificationDiagnostics diagnostic;
        double classify=0,validate=0,sampling=0,recovery=0,prepare=0,boundary=0;
        GeneralHybridReport result;unsigned attempts=0,nodes=0,segments=0,refinements=0,charts=0,chartIterations=0;bool any=false;
        auto classifyStart=std::chrono::steady_clock::now();
        auto pool=InitializeGeneralCandidates(samples,options.general,diagnostic);classify+=Ms(classifyStart);
        for(;;)
        {
            classifyStart=std::chrono::steady_clock::now();
            auto candidates=NextGeneralCandidateBatch(pool,samples,options.general,diagnostic);classify+=Ms(classifyStart);
            if(candidates.empty())break;
            for(const auto &candidate:candidates)
            {
                ++attempts;result=ValidateGeneralHybridWithSamples(*fixture.input,candidate.parameters,options,samples);
                CheckReport(result);validate+=result.totalMs;sampling+=result.samplingMs;recovery+=result.recoveryMs;
                prepare+=result.preparationMs;boundary+=result.boundaryMs;
                nodes+=result.boundaryNodes;segments+=result.boundarySegments;refinements+=result.mappingRefinements;
                charts+=result.chartAttempts;chartIterations+=result.chartIterations;
                if(result.Passed()){any=true;break;}
                if(result.status==HybridValidationStatus::InvalidInput)break;
            }
            if(any||result.status==HybridValidationStatus::InvalidInput)break;
        }
        double pipeline=Ms(start);FixtureCertificate audit;
        std::cout<<"END_PIPELINE "<<run<<std::endl;
        if(any)
        {
            std::cout<<"BEGIN_OUTPUT_AUDIT "<<run<<std::endl;
            auto check=fixture;check.truth=result.patch->Parameters();check.corners=result.patch->Corners();check.bend=0;
            audit=CertifyAnalyticPair(check,1e-5);
        }
        const auto &s=*fixture.input;
        out<<id<<','<<extra<<','<<FixtureGroupName(fixture.group)<<','<<TypeName(fixture.truth.type)<<','<<fixture.rational<<','
            <<s.ControlPoints().size()*s.ControlPoints().front().size()<<','<<run<<",1,"<<HybridValidationStatusName(result.status)<<','
            <<any<<','<<(any?TypeName(result.patch->Parameters().type):"Unknown")<<','<<audit.passed<<','<<audit.upper<<','<<attempts<<','<<result.interiorSampleCount<<','<<nodes<<','
            <<segments<<','<<refinements<<','<<result.collapsedEdges<<','<<charts<<','<<chartIterations<<','
            <<classify<<','<<sampling<<','<<recovery<<','<<prepare<<','<<boundary<<','<<validate<<','<<pipeline<<','
            <<result.interiorSampleMaxDistance<<','<<result.interiorSampleMaxFiniteDistance<<','<<result.boundaryDistanceUpperBound<<','
            <<result.boundaryTangentUpperBoundDeg<<','<<Quote(result.message)<<'\n';out.flush();
        Require(!any||audit.passed,"Accepted output failed independent audit");
    }
    return 0;
}
}
int main(int argc,char **argv)
{
    sggk::init();int code=2;
    try
    {
        Require(argc>=3,"Usage: --unit data | --case data csv id extra repetitions");
        if(std::string(argv[1])=="--unit"){Unit(argv[2]);code=0;}
        else {Require(argc==7,"Invalid case arguments");code=Case(argv[2],argv[3],unsigned(std::stoul(argv[4])),std::stoi(argv[5])!=0,unsigned(std::stoul(argv[6])));}
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';}
    return code;
}
