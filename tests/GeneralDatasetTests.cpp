#include "GeneralDataset.h"
#include "GeneralInternal.h"
#include <Foundation/init.h>
#include <Geometry/2D/BSplineCurve2D.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace bspline2analytic;
using namespace bspline2analytic::tests;
namespace
{
void Require(bool ok,const std::string &message){if(!ok)throw std::runtime_error(message);}
std::string Quote(std::string text){std::string out="\"";for(char c:text){if(c=='"')out+='"';out+=c=='\n'||c=='\r'?' ':c;}return out+'"';}
void CheckMembership(const std::filesystem::path &dir,bool extra)
{
    unsigned bins[5][5][2]{},collapsed[5]{},negativeKinds[6][2]{},groups[8]{};
    for(unsigned i=0;i<(extra?60u:600u);++i)
    {
        auto c=ReadGeneralFixture(dir/GeneralFixtureFilename(i,extra));
        Require(c.id==i && c.extra==extra,"Fixture identity mismatch");++groups[unsigned(c.group)];
        if(extra)continue;
        if(i<500)
        {
            Require(c.group==GeneralFixtureGroup::Feasible&&c.type<5,"Positive label mismatch");
            const auto &cp=c.input->ControlPoints();size_t count=cp.size()*cp[0].size();
            const unsigned lo[]={50,80,110,140,170},hi[]={79,109,139,169,200};
            Require(count>=lo[c.bucket]&&count<=hi[c.bucket],"Control point quota mismatch");
            ++bins[c.type][c.bucket][c.rational?0:1];if(c.collapsedMask)++collapsed[c.type];
            Require(c.bend==0&&c.pairTolerance<=1e-6,"Main positive oracle weakened");
        }
        else
        {
            Require(c.group==GeneralFixtureGroup::Infeasible&&c.exclusionRadius>=2e-5,"Negative label/radius mismatch");
            ++negativeKinds[c.type][c.rational?0:1];
        }
    }
    if(extra)
    {
        Require(groups[2]==20&&groups[3]==10&&groups[4]==10&&groups[5]==10&&groups[6]==5&&groups[7]==5,"Supplementary group quota mismatch");return;
    }
    for(unsigned t=0;t<5;++t)
    {
        Require(collapsed[t]==((t==2||t==3)?20u:0u),"Point-degenerate quota mismatch");
        for(unsigned b=0;b<5;++b)for(unsigned r=0;r<2;++r)Require(bins[t][b][r]==10,"Positive type/budget/representation quota mismatch");
    }
    for(unsigned t=0;t<6;++t)for(unsigned r=0;r<2;++r)Require(negativeKinds[t][r]==(t==5?20u:6u),"Negative quota mismatch");
}
FixtureCertificate InputCertificate(const GeneralFixture &c)
{
    return c.group==GeneralFixtureGroup::Infeasible?CertifyNoQuartic(*c.input,c.exclusionRadius):CertifyAnalyticPair(c,c.pairTolerance);
}
FixtureCertificate AuditOutput(const GeneralFixture &c,const GeneralConversionResult &r)
{
    FixtureCertificate out;
    if(!r.patch || !r.validation.Passed()){out.message="Success lacks finite patch or continuous production certificate";return out;}
    const auto &p=*r.patch;
    GeneralFixture recovered=c;recovered.truth=p.Parameters();recovered.corners=p.Corners();recovered.bend=0;
    out=CertifyAnalyticPair(recovered,1e-5);
    if(!out.passed)return out;
    const unsigned ends[4][2]={{0,2},{1,3},{0,1},{2,3}};
    for(unsigned e=0;e<4;++e)
    {
        const auto &curve=p.PCurves()[e];auto spline=std::dynamic_pointer_cast<sggk::BSplineCurve2D>(curve.sdkCurve);
        if(!spline || spline->Degree()!=1){out.passed=false;out.message="SDK PCurve is not degree one";return out;}
        auto a=p.Corners()[ends[e][0]],b=p.Corners()[ends[e][1]];
        if(curve.start.X()!=a.X()||curve.start.Y()!=a.Y()||curve.end.X()!=b.X()||curve.end.Y()!=b.Y())
        {out.passed=false;out.message="Canonical PCurves do not share exact corners";return out;}
        for(unsigned j=0;j<=16;++j)
        {
            double t=double(j)/16;auto uv=spline->CalcPoint(spline->MinParam()+(spline->MaxParam()-spline->MinParam())*t);
            auto world=p.MotherSurface()->CalcPoint(uv.X(),uv.Y());
            double error=detail::Norm(detail::Get(world)-detail::Get(p.CalcBoundaryPoint(BoundaryId(e),t)));
            if(!(error<1e-7)){out.passed=false;out.message="SDK/canonical PCurve evaluation mismatch";return out;}
        }
    }
    return out;
}
void CacheChecks(const GeneralFixture &c)
{
    detail::SurfaceSampleContext context(*c.input);
    Require(context.Evaluations()==0,"Sampling context eagerly evaluates source");
    auto p=context.At(.37,.61,sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
    auto count=context.Evaluations();
    auto q=context.At(.37,.61).point;
    Require(count==context.Evaluations()&&detail::Norm(p-q)==0,"Derivative upgrade duplicates cached position");
    context.Greville(false);count=context.Evaluations();context.Greville(false);
    Require(count==context.Evaluations(),"Repeated Greville view evaluates again");
    const auto coarse=context.CoarseScores(16);count=context.Evaluations();
    Require(coarse==context.CoarseScores(16)&&count==context.Evaluations(),"Coarse view evaluates source again");
    Require(coarse.size()==std::min(size_t(16),context.Scores().size()),"Wrong coarse subset size");
    for(size_t i:coarse)
        Require(std::count(coarse.begin(),coarse.end(),i)==1&&
            std::find(context.Scores().begin(),context.Scores().end(),i)!=context.Scores().end(),"Coarse points are not distinct held-out scores");
    Require(context.CoarseScores(20000)==context.Scores(),"Small score set is not retained in full");
    context.Boundaries();count=context.Evaluations();context.Boundaries();
    Require(count==context.Evaluations(),"Repeated boundary view evaluates again");
    Require(!CertifyNoQuartic(*c.input,2e-5).passed,"Negative oracle excluded a certified positive");
}
// Exhaust both pools without converting the wrong mothers. This exercises the
// backup scheduler on positive inputs even when conversion succeeds immediately.
void ScoringChecks(const GeneralFixture &c)
{
    using namespace detail;
    SurfaceSampleContext stagedContext(*c.input),fullContext(*c.input);
    GeneralConversionOptions staged,full;full.coarseScorePointCount=0;
    GeneralClassificationDiagnostics a,b;
    auto pool=InitializeGeneralCandidates(stagedContext,staged,a);
    auto reference=InitializeGeneralCandidates(fullContext,full,b);
    Require(a.seedAttemptCount==40&&b.seedAttemptCount==40,"Seed budget was refilled after failure/deduplication");
    unsigned seen=0,batches=0;
    for(;;)
    {
        auto x=NextGeneralCandidateBatch(pool,stagedContext,staged,a);
        auto y=NextGeneralCandidateBatch(reference,fullContext,full,b);
        Require(x.size()==y.size(),"Staged/eager batch sizes differ: fixture="+std::to_string(c.id)+
            ", batch="+std::to_string(batches)+", sizes="+std::to_string(x.size())+"/"+std::to_string(y.size())+
            ", unique="+std::to_string(a.scoreEvaluationCount)+"/"+std::to_string(b.scoreEvaluationCount)+
            ", completed="+std::to_string(a.completedScoreCount)+"/"+std::to_string(b.completedScoreCount));
        if(x.empty())break;
        ++batches;seen+=unsigned(x.size());
        for(size_t i=0;i<x.size();++i)
        {
            const auto &p=x[i].parameters,&q=y[i].parameters;
            Require(p.type==q.type&&Norm(Get(p.origin)-Get(q.origin))==0&&Norm(Get(p.axis)-Get(q.axis))==0&&
                p.radius==q.radius&&p.majorRadius==q.majorRadius&&p.halfAngle==q.halfAngle,"Staged ranking differs from full scoring");
            const auto &s=x[i].score,&t=y[i].score;
            Require(s.interiorMax==t.interiorMax&&s.boundaryMax==t.boundaryMax&&s.normalMax==t.normalMax&&
                std::abs(s.interiorRms-t.interiorRms)<=1e-13*std::max(1.0,t.interiorRms),"Completed score differs from full scoring");
        }
    }
    Require(seen==a.completedScoreCount&&batches==a.batchCount,"Backup pool skipped/repeated a scored candidate");
    Require(a.seedAttemptCount==40,"Backup pool generated new seed groups");
    auto coarse=stagedContext.CoarseScores(16).size(),scores=stagedContext.Scores().size();
    Require(a.interiorPointEvaluations<=a.coarseScoreCount*coarse+a.fullScoreCount*(scores-coarse),"Coarse points were scored twice");
    for(const auto &family:pool.families)for(const auto &candidate:family)
        Require(candidate.stage==GeneralScoreStage::Consumed||candidate.stage==GeneralScoreStage::Failed,"Unexhausted backup candidate");
    std::cout<<"Scoring check "<<c.id<<": "<<seen<<" candidates, "<<batches<<" batches, backups="<<a.backupCandidateCount<<std::endl;
}
void BoundaryAxisChecks(const GeneralFixture &c)
{
    using namespace detail;
    const auto &s=*c.input;
    auto cells=ExtractCells(s);
    auto check=[&](const BezierCell &b,unsigned edge,double t) {
        bool along=edge>=2;
        double u=along?b.u0+(b.u1-b.u0)*t:b.u0;
        double v=along?b.v0:b.v0+(b.v1-b.v0)*t;
        auto p=Range(Evaluate(b.p,along?t:0,along?0:t));
        I w=Evaluate(b.w,along?t:0,along?0:t).range();
        auto center=[](I a) {return (a.lo+a.hi)*.5;};
        V actual{center(p.x/w),center(p.y/w),center(p.z/w)};
        sggk::Point3D point;
        auto d=s.CalcPartialD1(u,v,u==s.MaxParamU()?sggk::DerivType::LEFT:sggk::DerivType::RIGHT,
            v==s.MaxParamV()?sggk::DerivType::LEFT:sggk::DerivType::RIGHT,point);
        Require(Norm(actual-Get(point))<1e-9,"Boundary endpoint/subdivision parameter mismatch");
        auto dp=Range(Evaluate(Derivative(b.p,along),along?t:0,along?0:t));
        I dw=Evaluate(Derivative(b.w,along),along?t:0,along?0:t).range();
        V derivative{center((dp.x*w-p.x*dw)/Sqr(w)),center((dp.y*w-p.y*dw)/Sqr(w)),center((dp.z*w-p.z*dw)/Sqr(w))};
        V expected=Get(d[along?0:1])*(along?b.u1-b.u0:b.v1-b.v0);
        Require(Norm(derivative-expected)<1e-8,"Boundary derivative uses wrong parameter axis");
    };
    for(const auto &cell:cells) for(unsigned edge=0;edge<4;++edge)
    {
        auto b=BoundaryCell(cell,BoundaryId(edge));
        check(b,edge,0);check(b,edge,1);check(b,edge,.37);
        auto children=Split(b,edge>=2);
        check(children.first,edge,.37);check(children.second,edge,.63);
    }
}
int Run(const std::filesystem::path &dir,const std::filesystem::path &csvPath,bool extra,unsigned first,unsigned last,bool certificatesOnly,bool compare=false)
{
    Require(!compare||(!extra&&!certificatesOnly&&last<=500),"Comparison is restricted to main positives");
    Require(first<last&&last<=(extra?60u:600u),"Invalid test range");CheckMembership(dir,extra);
    if(!extra&&first==0)
    {
        CacheChecks(ReadGeneralFixture(dir/GeneralFixtureFilename(0)));
        for(unsigned id:{0u,100u,200u,216u,300u,316u,339u,400u})
            ScoringChecks(ReadGeneralFixture(dir/GeneralFixtureFilename(id)));
        BoundaryAxisChecks(ReadGeneralFixture(dir/GeneralFixtureFilename(0)));
        for(unsigned type=1;type<5;++type)
        {
            auto positive=ReadGeneralFixture(dir/GeneralFixtureFilename(type*100));
            BoundaryAxisChecks(positive);
            Require(!CertifyNoQuartic(*positive.input,2e-5).passed,"Negative oracle excluded a curved positive");
        }
        for(unsigned id:{216u,217u,218u,219u,316u,317u,318u,319u})
            BoundaryAxisChecks(ReadGeneralFixture(dir/GeneralFixtureFilename(id)));
    }
    std::ofstream csv(csvPath);Require(bool(csv),"Cannot open test CSV");csv<<std::setprecision(17);
    if(certificatesOnly)csv<<"id,extra,file,certificate_passed,pair_bound,rank_margin,message\n";
    else csv<<"id,extra,file,group,source_type,rational,controls,collapsed_mask,fixture_certified,fixture_bound,rank_margin,expected,success,passed,output_type,output_audit,output_bound,status,pcurve_failure,candidates,attempts,iterations,validation_pieces,classify_ms,optimize_ms,validate_ms,total_ms,sphere_attempted,torus_attempted,known_candidate_passed,chart_attempts,chart_iterations,raw_candidates,duplicate_candidates,score_evaluations,configuration,seed_attempts,coarse_points,coarse_scores,full_scores,completed_scores,deferred_candidates,backup_candidates,batches,interior_point_evaluations,boundary_point_evaluations,message\n";
    unsigned passed=0,success=0,failures=0;
    for(unsigned i=first;i<last;++i)
    {
        auto c=ReadGeneralFixture(dir/GeneralFixtureFilename(i,extra));auto certificate=InputCertificate(c);
        Require(certificate.passed,"Fixture oracle failed: "+GeneralFixtureFilename(i,extra)+" "+certificate.message);
        if(certificatesOnly)
        {
            ++passed;csv<<i<<','<<extra<<','<<GeneralFixtureFilename(i,extra)<<",1,"<<certificate.upper<<','<<certificate.rankMargin<<','<<Quote(certificate.message)<<'\n';
            std::cout<<"Certified "<<i<<' '<<certificate.message<<std::endl;continue;
        }
        // Rotate all six orders across fixtures; fresh contexts inside each
        // ConvertGeneral call, with input/output certificates outside timing.
        const unsigned orders[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        for(unsigned run=0;run<(compare?3u:1u);++run)
        {
            unsigned configuration=compare?orders[i%6][run]:2;
            const char *configurationName[]={"64-full","8-full","8-staged"};
            GeneralConversionOptions options;
            if(configuration==0)options.seedAttemptsPerType=64;
            if(configuration<2)options.coarseScorePointCount=0;
            if(c.group==GeneralFixtureGroup::InvalidOptions)
            {
                if(i%2)options.fitting.validation.distanceTolerance=-1;else options.fitting.sampling.maxSamples=0;
            }
            if(c.group==GeneralFixtureGroup::ValidationBudget){options.fitting.validation.maxValidationPieces=1;options.fitting.refinement.maxRounds=0;}
            // Policy probes deliberately broaden 'similar' to exercise the retained
            // sphere -> torus path; main 600 always use unchanged default options.
            if(c.group==GeneralFixtureGroup::SpherePriority)options.similarScoreAbsoluteRatio=1e6;
            std::cout<<"BEGIN "<<GeneralFixtureFilename(i,extra)<<' '<<FixtureGroupName(c.group)<<' '<<configurationName[configuration]<<std::endl;
            auto result=ConvertGeneral(*c.input,options); // ONLY conversion entry used by this suite.
            bool succeeded=result.status==ConversionStatus::Success;success+=succeeded;
            FixtureCertificate audit;if(succeeded)audit=AuditOutput(c,result);
            bool sphere=false,torus=false,sphereBeforeTorus=false;
            for(const auto &attempt:result.attempts){if(attempt.type==sggk::SurfaceType::Sphere)sphere=true;if(attempt.type==sggk::SurfaceType::Torus){torus=true;sphereBeforeTorus|=sphere;}}
            bool expectedSuccess=c.group==GeneralFixtureGroup::Feasible||c.group==GeneralFixtureGroup::Periodic||c.group==GeneralFixtureGroup::SpherePriority||c.group==GeneralFixtureGroup::Tolerance;
            bool ok=expectedSuccess ? succeeded&&audit.passed : !succeeded&&!result.patch;
            if(!extra&&i<500&&succeeded)ok&=result.patch->Parameters().type==GeneralDatasetType(c.type);
            int known=-1;
            if(c.group==GeneralFixtureGroup::NonlinearPCurve||c.group==GeneralFixtureGroup::Tolerance)
            {
                detail::SurfaceSampleContext context(*c.input);GeneralConversionOptions fixed;
                if(c.group==GeneralFixtureGroup::Tolerance&&i==49)
                {
                    auto origin=detail::Get(c.truth.origin),normal=detail::Get(c.truth.axis);
                    double sampled=0,peak=0;
                    for(const auto &point:context.Greville(false))sampled=std::max(sampled,std::abs(detail::Dot(point.point-origin,normal)));
                    for(unsigned k=0;k<=10000;++k)peak=std::max(peak,std::abs(detail::Dot(detail::Get(c.input->CalcPoint(double(k)/10000,.5))-origin,normal)));
                    Require(sampled<1e-5&&peak>1e-5,"Blind-spot fixture no longer straddles tolerance");
                }
                auto report=detail::ValidateGeneral(*c.input,c.truth,fixed,context);
                known=report.Passed()?1:0;
                if(c.group==GeneralFixtureGroup::NonlinearPCurve)
                {
                    // This label concerns the specified mother. Another independently
                    // certified analytic solution is allowed; it is not a global negative.
                    ok=!report.Passed()&&report.pcurveFailure==PCurveFailure::NonlinearPCurve&&(!succeeded||audit.passed);
                }
                else ok&=(i<45?report.Passed():!report.Passed());
            }
            if(c.group==GeneralFixtureGroup::SpherePriority)
            {
                if(i<35)ok&=succeeded&&result.patch->Parameters().type==sggk::SurfaceType::Sphere;
                else ok&=sphereBeforeTorus&&succeeded&&result.patch->Parameters().type==sggk::SurfaceType::Torus;
            }
            if(c.group==GeneralFixtureGroup::InvalidOptions)ok&=result.status==ConversionStatus::InvalidInput;
            if(succeeded)ok&=result.validation.distanceUpperBound<1e-5&&result.validation.boundaryHausdorffUpperBound<1e-5;
            if(ok)++passed;else ++failures;
            std::string message=result.message;if(succeeded&&!audit.passed)message+="; independent output audit: "+audit.message;
            const auto &time=result.timing;const auto &statistics=result.statistics;
            csv<<i<<','<<extra<<','<<GeneralFixtureFilename(i,extra)<<','<<FixtureGroupName(c.group)<<','<<(c.type<5?TypeName(GeneralDatasetType(c.type)):"Freeform")<<','<<c.rational<<','<<c.input->ControlPoints().size()*c.input->ControlPoints()[0].size()<<','<<c.collapsedMask<<",1,"<<certificate.upper<<','<<certificate.rankMargin<<','<<(expectedSuccess?"success":c.group==GeneralFixtureGroup::NonlinearPCurve?"reject_known_pcurve":"failure")<<','<<succeeded<<','<<ok<<','<<(succeeded?TypeName(result.patch->Parameters().type):"Unknown")<<','<<audit.passed<<','<<audit.upper<<','<<unsigned(result.status)<<','<<PCurveFailureName(result.validation.pcurveFailure)<<','<<result.candidates.size()<<','<<result.attempts.size()<<','<<statistics.optimizationIterations<<','<<statistics.validationPieces<<','<<time.classification.elapsedMs<<','<<time.optimization.elapsedMs<<','<<time.validation.elapsedMs<<','<<time.totalMs<<','<<sphere<<','<<torus<<','<<known<<','<<result.chartAttempts<<','<<result.chartIterations<<','<<result.classification.rawCandidateCount<<','<<result.classification.duplicateCandidateCount<<','<<result.classification.scoreEvaluationCount<<','<<configurationName[configuration]<<','<<result.classification.seedAttemptCount<<','<<options.coarseScorePointCount<<','<<result.classification.coarseScoreCount<<','<<result.classification.fullScoreCount<<','<<result.classification.completedScoreCount<<','<<result.classification.deferredCandidateCount<<','<<result.classification.backupCandidateCount<<','<<result.classification.batchCount<<','<<result.classification.interiorPointEvaluations<<','<<result.classification.boundaryPointEvaluations<<','<<Quote(message)<<'\n';csv.flush();
            std::cout<<"END "<<i<<" expected="<<expectedSuccess<<" success="<<succeeded<<" passed="<<ok<<" ms="<<time.totalMs<<" "<<message<<std::endl;
        }
    }
    Require(bool(csv),"CSV write failed");std::cout<<"Completed "<<(last-first)*(compare?3:1)<<"; passed="<<passed<<"; successes="<<success<<"; mismatches="<<failures<<std::endl;
    return failures?1:0;
}
void Diagnose(const std::filesystem::path &dir,unsigned id)
{
    auto c=ReadGeneralFixture(dir/GeneralFixtureFilename(id));
    auto certificate=InputCertificate(c);Require(certificate.passed,"Diagnostic fixture is uncertified");
    GeneralConversionOptions options;
    std::cout<<std::setprecision(17)<<"fixture="<<id<<" pair_bound="<<certificate.upper<<'\n';
    if(c.group==GeneralFixtureGroup::Feasible)
    {
        detail::SurfaceSampleContext context(*c.input);
        auto known=detail::ValidateGeneral(*c.input,c.truth,options,context);
        std::cout<<"known_mother passed="<<known.Passed()<<" geometry="<<unsigned(known.geometry)<<" boundary="<<unsigned(known.boundaryHausdorff)<<" pcurve="<<unsigned(known.pcurveLinearity)<<" failure="<<PCurveFailureName(known.pcurveFailure)<<" message="<<known.message<<'\n';
    }
    // Separate black-box call. The known mother above is never supplied to it.
    {
        detail::SurfaceSampleContext context(*c.input);
        auto seeds=detail::ClassifyGeneral(context,options);
        if(!seeds.empty())
        {
            detail::GeneralTrim trim;GeneralValidationReport report;
            bool recovered=detail::RecoverLinearPCurves(context,options,seeds.front().parameters,trim,report);
            std::cout<<"initial_chart="<<recovered<<" message="<<report.message<<'\n';
            if(!recovered) for(unsigned edge=0;edge<4;++edge)
            {
                const auto &line=trim.edges[edge];
                std::cout<<"edge "<<edge<<" collapsed="<<line.collapsed3D<<" start="<<line.start.X()<<','<<line.start.Y()<<" end="<<line.end.X()<<','<<line.end.Y()<<'\n';
                const auto &walk=trim.walks[edge];
                if(!walk.empty()) for(size_t j:{size_t(0),walk.size()/2,walk.size()-1})
                    std::cout<<"  uv "<<walk[j].first<<' '<<walk[j].second[0]<<','<<walk[j].second[1]<<'\n';
            }
        }
    }
    auto result=ConvertGeneral(*c.input,options);
    for(size_t i=0;i<result.candidates.size();++i)
    {
        const auto &cnd=result.candidates[i];std::cout<<"seed "<<i<<' '<<TypeName(cnd.type)<<" interior="<<cnd.interiorMax<<" boundary="<<cnd.boundaryMax<<'\n';
    }
    for(size_t i=0;i<result.attempts.size();++i)
    {
        const auto &attempt=result.attempts[i];std::cout<<"attempt "<<i<<' '<<TypeName(attempt.type)<<" success="<<attempt.succeeded<<" failure="<<PCurveFailureName(attempt.pcurveFailure)<<" message="<<attempt.message<<'\n';
    }
}
}
int main(int argc,char **argv)
{
    sggk::init();int code=2;
    try
    {
        Require(argc>=3,"Usage: GeneralDatasetTests --generate dir [first last] | --run/--extras/--certify/--compare-positive dir csv [first last] | --diagnose dir id");
        std::string mode=argv[1];
        if(mode=="--generate"){GenerateGeneralDataset(argv[2],argc>3?unsigned(std::stoul(argv[3])):0,argc>4?unsigned(std::stoul(argv[4])):660);code=0;}
        else if(mode=="--diagnose"){Require(argc==4,"Expected --diagnose directory id");Diagnose(argv[2],unsigned(std::stoul(argv[3])));code=0;}
        else
        {
            Require(argc>=4&&(mode=="--run"||mode=="--extras"||mode=="--certify"||mode=="--compare-positive"),"Invalid test command");
            bool extra=mode=="--extras";code=Run(argv[2],argv[3],extra,argc>4?unsigned(std::stoul(argv[4])):0,argc>5?unsigned(std::stoul(argv[5])):extra?60:mode=="--compare-positive"?500:600,mode=="--certify",mode=="--compare-positive");
        }
    }
    catch(const std::exception &e){std::cerr<<"GENERAL DATASET ERROR: "<<e.what()<<std::endl;}
    sggk::fini();return code;
}
