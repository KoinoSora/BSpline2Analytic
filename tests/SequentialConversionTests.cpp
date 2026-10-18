/**
 * @file SequentialConversionTests.cpp
 * @brief 统一串联入口的分支、限时与逐例结果驱动。
 * 先认证输入，计时仅包围完整 API，随后独立审计输出；预算无效、等参失败和真实返回超时分别检查。
 */
#include "SequentialConversion.h"
#include "Dataset500.h"
#include "GeneralRationalDataset.h"
#include "MotherDistanceInternal.h"
#include "GeneralRationalBoundaryInternal.h"
#include <Foundation/init.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace rb=bspline2analytic::detail::rational_boundary;
namespace
{
using Row=std::map<std::string,std::string>;
void Require(bool b,const std::string &m) { if(!b) throw std::runtime_error(m); }
template<class T> void Set(Row &r,const std::string &k,T v)
{ std::ostringstream s;s<<std::setprecision(17)<<v;r[k]=s.str(); }
#include "SequentialOutputAudit.h"
double Ms(std::chrono::steady_clock::time_point t)
{return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count();}
std::string Quote(std::string s)
{std::string q="\"";for(char c:s){if(c=='"')q+='"';q+=(c=='\n'||c=='\r')?' ':c;}return q+'"';}
int Run(const std::filesystem::path &root,const std::string &suite,unsigned id,
        unsigned runs,const std::filesystem::path &csv,bool unit)
{
    sggk::BSplineSurfacePtr input;
    sggk::SurfaceType truth;
    bool certified=false; double sourceBound=0;
    if(suite=="iso")
    {
        auto c=ReadDatasetCase(root/"dataset500"/CaseFilename(id));
        CheckDataset500Contract(c,id);
        auto cert=BoundPair(*c.source,*c.input,c.refitted?1e-6:1e-9);
        certified=cert.passed;sourceBound=cert.upper;input=c.input;truth=c.truth.type;
    }
    else
    {
        Require(suite=="main"||suite=="extra","Non-extreme suites only");
        auto c=ReadGeneralFixture(root/"general-rational-v2"/suite/RationalFixtureName(suite,id));
        Require(c.rationalDomain&&c.id==id,"V2 input required");
        auto cert=CertifyRationalDomain(c);
        certified=cert.passed;sourceBound=cert.upper;input=c.input;truth=c.truth.type;
    }
    Require(certified,"Independent source certificate failed");
    if(unit)
    {
        auto invalid=SequentialConversionOptions{};invalid.maxTotalMs=51;
        auto a=ConvertIsoThenFallback(*input,invalid);
        Require(!a.Passed()&&!a.isoAttempted&&!a.fallbackAttempted,"Invalid budget accepted");
        invalid.maxTotalMs=0;
        Require(!ConvertIsoThenFallback(*input,invalid).Passed(),"Zero budget accepted");
        invalid.maxTotalMs=1e-9;
        auto b=ConvertIsoThenFallback(*input,invalid);
        Require(b.timeBudgetExceeded&&!b.Passed()&&!b.isoPatch&&!b.fallback.patch&&!b.exact.surface,"Deadline published success/output");
        SequentialConversionOptions normal;
        if(suite=="main"&&id==100)
        { normal.exact.distanceTolerance=-1; normal.iso.validation.distanceTolerance=-1; }
        auto c=ConvertIsoThenFallback(*input,normal);
        Require(c.isoAttempted,"Iso branch skipped");
        if(suite=="iso"&&id==0)Require(c.Passed()&&!c.fallbackAttempted,"Successful iso did not short circuit");
        if(suite=="main"&&id==100)Require(c.fallbackAttempted,"Failed iso did not enter fallback");
        std::cout<<"Sequential unit checks passed: "<<suite<<':'<<id<<'\n';return 0;
    }
    std::ofstream out(csv);Require(bool(out),"Cannot open CSV");
    out<<"suite,id,run,source_certified,source_bound,source_type,output_type,rational,iso_attempted,fallback_attempted,path,accepted,time_budget_exceeded,timeout_stage,iso_ms,fallback_ms,api_ms,total_ms,under_10ms,within_50ms,iso_candidates,fallback_candidates,optimization_iterations,audit,iso_message,message,iso_budget_exceeded\n";
    out<<std::setprecision(17);
    for(unsigned run=0;run<runs;++run)
    {
        // 该区间严禁日志、文件读写、真值检查或输出审计，只测完整算法调用。
        const auto start=std::chrono::steady_clock::now();
        auto r=ConvertIsoThenFallback(*input);
        const double elapsed=Ms(start);
        if(elapsed>50){r.accepted=false;r.timeBudgetExceeded=true;r.timeoutStage="api_return_overrun";}
        Require(r.isoAttempted,"Every input must try iso first");
        Require(r.Passed()||r.fallbackAttempted||r.timeBudgetExceeded,"Failure skipped fallback");
        std::string audit="not_applicable",auditMessage,outputType;
        if(r.Passed())
        {
            try
            {
                Row values;
                if(r.path=="general_rational")
                {
                    Require(r.fallback.Passed(),"Missing fallback certificate");
                    outputType=TypeName(r.fallback.patch->mother.type);
                    AuditGeneral(*input,r.fallback,1e-6);
                    // 补充采样审计，不声称内部连续域已认证。
                    AuditPosition(*input,r.fallback.patch->mother,1e-5,values);
                }
                else if(r.path=="iso_mother")
                {
                    Require(r.isoPatch&&r.mother.Passed(),"Missing iso certificate");
                    outputType=TypeName(r.isoPatch->parameters.type);
                    AuditPosition(*input,r.isoPatch->parameters,r.mother.motherDistance.upperBound,values);
                    AuditEdges(*input,*r.isoPatch,r.mother,values);
                }
                else
                {
                    Require(r.exact.surface&&r.exact.domain&&r.exact.validation.level==VerificationLevel::BernsteinBound,"Missing exact certificate");
                    outputType=TypeName(r.exact.parameters.type);
                    AuditPosition(*input,r.exact.parameters,r.exact.validation.surfaceDistanceBound,values);
                }
                if(suite=="iso"||(suite=="main"&&id<500))
                    Require(outputType==TypeName(truth),"Unexpected output type");
                audit="passed";
            }
            catch(const std::exception &e){audit="failed";auditMessage=e.what();}
        }
        out<<suite<<','<<id<<','<<run<<','<<certified<<','<<sourceBound<<','<<TypeName(truth)<<','<<outputType<<','<<input->IsRational()<<','
           <<r.isoAttempted<<','<<r.fallbackAttempted<<','<<r.path<<','<<r.Passed()<<','<<r.timeBudgetExceeded<<','<<Quote(r.timeoutStage)<<','
           <<r.isoMs<<','<<r.fallbackMs<<','<<r.totalMs<<','<<elapsed<<','<<(elapsed<10)<<','<<(elapsed<=50)<<','
           <<r.isoCandidates<<','<<r.fallbackCandidates<<','<<r.optimizationIterations<<','<<audit<<','<<Quote(r.isoMessage)<<','<<Quote(r.message+"; "+auditMessage)<<','<<r.isoTimeBudgetExceeded<<'\n';
        out.flush();
    }
    return 0; // Python 调度器负责检查全部预期，包括负例对照。
}
}
int main(int argc,char **argv)
{
    if(argc!=7){std::cerr<<"data suite id runs csv run|unit\n";return 2;}
    sggk::init();int result=2;
    try{result=Run(argv[1],argv[2],unsigned(std::stoul(argv[3])),unsigned(std::stoul(argv[4])),argv[5],std::string(argv[6])=="unit");}
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';}
    sggk::fini();return result;
}
