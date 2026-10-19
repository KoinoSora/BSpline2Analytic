/** @file AnalyticExclusionTests.cpp
 * @brief 排除证书正例保护、容差保护、局部变形与一次调度回归。
 */
#include "AnalyticExclusionInternal.h"
#include "SequentialConversion.h"
#include "ReferenceCases.h"
#include <Foundation/init.h>
#include <cfenv>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
namespace
{
unsigned checks=0;
void Require(bool b,const char *message){++checks;if(!b)throw std::runtime_error(message);}
sggk::BSplineSurfacePtr Copy(const sggk::BSplineSurface &s,sggk::Point3DMatrix cp)
{
    return std::make_shared<sggk::BSplineSurface>(s.DegreeU(),s.DegreeV(),cp,s.Weights(),s.KnotsU(),s.KnotsV(),s.MultsU(),s.MultsV());
}
AnalyticExclusionReport Check(const sggk::BSplineSurface &s,double tolerance=1e-5,double local=1000)
{
    ExecutionDeadline deadline(0);AnalyticExclusionReport r;
    CertifyAnalyticExclusion(s,tolerance,deadline,r,local);return r;
}
sggk::BSplineSurfacePtr Bump(double scale,double translation)
{
    sggk::Point3DMatrix cp(7,std::vector<sggk::Point3D>(7));
    sggk::RealMatrix weights(7,std::vector<double>(7,1));
    unsigned state=12345;
    for(unsigned i=0;i<7;++i)for(unsigned j=0;j<7;++j)
    {
        state=1664525u*state+1013904223u;
        double z=(i&&i<6&&j&&j<6)?double(state%1001)/1000.-.5:0.;
        // 固定正交旋转、均匀缩放和平移；边界仍共面，内部有局部变形。
        double x=double(i)/3.-1,y=double(j)/3.-1;
        cp[i][j]=sggk::Point3D(translation+scale*(.6*x-.8*y),-translation+scale*(.8*x+.6*y),translation+scale*z);
    }
    return std::make_shared<sggk::BSplineSurface>(3,3,cp,weights,sggk::RealArray{0,.5,1},sggk::RealArray{0,.5,1},sggk::UIntArray{4,3,4},sggk::UIntArray{4,3,4});
}
void Run()
{
    auto cases=tests::MakeCases();
    for(const auto &c:cases)
    {
        auto r=Check(*c.spline);
        Require(!r.Excluded(),"Exact analytic surface falsely excluded");
        Require(r.status==AnalyticExclusionStatus::Inconclusive,"Exact test did not complete proof attempt");
        auto points=c.spline->ControlPoints();unsigned i=0;
        for(auto &row:points)for(auto &p:row){p=sggk::Point3D(p.X(),p.Y(),p.Z()+(i++%2?1:-1)*2.5e-6);}
        auto near=Copy(*c.spline,points);
        Require(!Check(*near).Excluded(),"Within-tolerance perturbation falsely excluded");
    }
    auto accepted=ConvertIsoThenFallback(*cases.front().spline);
    Require(accepted.Passed()&&!accepted.fallbackAttempted&&accepted.exclusion.attempts==0,"Iso success did not skip exclusion");
    for(auto values:std::array<std::array<double,2>,3>{{{1.,0.},{1000.,1e6},{.1,-1e4}}})
    {
        auto source=Bump(values[0],values[1]);const double tolerance=1e-5*values[0];auto proof=Check(*source,tolerance);
        std::cout<<"Deformation scale="<<values[0]<<" lower="<<proof.distanceLowerBound<<" rho="<<proof.residualNormUpper<<" reason="<<proof.message<<'\n';
        Require(proof.Excluded()&&proof.distanceLowerBound>tolerance,"Interior deformation not excluded");
        Require(!Check(*source,proof.distanceLowerBound*2).Excluded(),"Configured tolerance ignored");
    }
    auto source=Bump(1,0);
    Require(Check(*source,1e-5,0).status==AnalyticExclusionStatus::LocalBudgetExceeded,"Local budget not conservative");
    auto rankDeficient=Copy(*cases.front().spline,cases.front().spline->ControlPoints());
    Require(!Check(*rankDeficient).Excluded(),"Rank deficiency rejected");
    std::fesetround(FE_DOWNWARD);auto unsupported=Check(*source);std::fesetround(FE_TONEAREST);
    Require(!unsupported.Excluded()&&unsupported.message.find("environment")!=std::string::npos,"Unsupported rounding accepted");
    bool deadlineRaised=false;ExecutionDeadline expired(1e-9);AnalyticExclusionReport aborted;
    try{CertifyAnalyticExclusion(*source,1e-5,expired,aborted);}catch(const DeadlineExceeded &){deadlineRaised=true;}
    Require(deadlineRaised&&!aborted.Excluded(),"Global deadline swallowed");
    auto r=ConvertIsoThenFallback(*source);
    Require(!r.Passed()&&!r.timeBudgetExceeded&&r.exclusion.Excluded()&&r.exclusion.attempts==1&&r.optimizationIterations==0,"Certified rejection did not precede optimization");
    Require(!r.isoPatch&&!r.exact.surface&&!r.fallback.patch,"Rejection published geometry");
    SequentialConversionOptions disabled;disabled.fallback.general.fitting.refinement.maxIterations=0;
    Require(ConvertIsoThenFallback(*source,disabled).exclusion.attempts==0,"No-optimization path ran exclusion");
    // 椭圆柱不属于五类，但属于排除器覆盖的二次多项式超集，必须无结论后继续优化。
    auto cp=cases[1].spline->ControlPoints();for(auto &row:cp)for(auto &p:row)p=sggk::Point3D(p.X()*1.5,p.Y(),p.Z());
    auto ellipse=Copy(*cases[1].spline,cp);auto fallback=ConvertIsoThenFallback(*ellipse);
    Require(!fallback.exclusion.Excluded()&&fallback.exclusion.attempts==1&&fallback.optimizationIterations>0,"Inconclusive gate did not resume original optimizer once");
    std::cout<<"Analytic exclusion checks passed: "<<checks<<'\n';
}
}
int main()
{
    try{sggk::init();Run();sggk::fini();return 0;}
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
