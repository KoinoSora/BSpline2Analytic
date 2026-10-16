#include "GeneralExtremeDataset.h"
#include "GeneralInternal.h"
#include <Foundation/init.h>
#include <Geometry/2D/BSplineCurve2D.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace bspline2analytic;
using namespace bspline2analytic::tests;
using namespace bspline2analytic::detail;
namespace
{
void Require(bool ok,const std::string &message){if(!ok)throw std::runtime_error(message);}
std::string Quote(std::string s){std::string out="\"";for(char c:s){if(c=='"')out+='"';out+=(c=='\n'||c=='\r')?' ':c;}return out+'"';}
void SeedDiagnostics(const std::filesystem::path &csv,const GeneralConversionResult &result)
{
    std::ofstream seeds(csv.parent_path()/(csv.stem().string()+"-seeds.csv"));
    std::ofstream attempts(csv.parent_path()/(csv.stem().string()+"-attempts.csv"));
    Require(bool(seeds)&&bool(attempts),"Cannot write seed diagnostics");
    seeds<<std::setprecision(17)<<"group,type,attempt,samples,multi_point,failure,generated,duplicates,scored,fully_scored,selected,relative_rank,axis_separation,construction_residual\n";
    const auto &groups=result.classification.seedGroups;
    for(size_t i=0;i<groups.size();++i)
    {
        const auto &d=groups[i];
        seeds<<i<<','<<TypeName(d.type)<<','<<d.attempt<<','<<d.sampleCount<<','<<d.multiPoint<<','
             <<GeneralSeedFailureName(d.failure)<<','<<d.generated<<','<<d.duplicates<<','<<d.scored<<','
             <<d.fullyScored<<','<<d.selected<<','<<d.relativeRank<<','<<d.axisSeparation<<','<<d.constructionResidual<<'\n';
    }
    attempts<<"type,seed_group,iterations,succeeded,pcurve_failure,message\n";
    for(const auto &a:result.attempts)
        attempts<<TypeName(a.type)<<','<<a.seedGroup<<','<<a.optimizationIterations<<','<<a.succeeded<<','
                <<PCurveFailureName(a.pcurveFailure)<<','<<Quote(a.message)<<'\n';
}
FixtureCertificate Audit(const GeneralFixture &c,const GeneralConversionResult &r)
{
    FixtureCertificate out;
    if(!r.patch||!r.validation.Passed()){out.message="Missing successful finite patch validation";return out;}
    auto recovered=c;recovered.truth=r.patch->Parameters();
    recovered.corners=r.patch->Corners();
    out=CertifyGeneralExtremePair(recovered,1e-5);
    if(!out.passed)return out;
    const unsigned ends[4][2]={{0,2},{1,3},{0,1},{2,3}};
    for(unsigned e=0;e<4;++e)
    {
        const auto &edge=r.patch->PCurves()[e];auto spline=std::dynamic_pointer_cast<sggk::BSplineCurve2D>(edge.sdkCurve);
        if(!spline||spline->Degree()!=1){out.passed=false;out.message="Nonlinear SDK PCurve";return out;}
        auto a=r.patch->Corners()[ends[e][0]],b=r.patch->Corners()[ends[e][1]];
        if(edge.start.X()!=a.X()||edge.start.Y()!=a.Y()||edge.end.X()!=b.X()||edge.end.Y()!=b.Y())
        {out.passed=false;out.message="PCurve endpoint mismatch";return out;}
        for(unsigned k=0;k<=16;++k)
        {
            double t=double(k)/16;auto uv=spline->CalcPoint(spline->MinParam()+(spline->MaxParam()-spline->MinParam())*t);
            double error=Norm(Get(r.patch->MotherSurface()->CalcPoint(uv.X(),uv.Y()))-Get(r.patch->CalcBoundaryPoint(BoundaryId(e),t)));
            if(!(error<1e-5)){out.passed=false;out.message="SDK PCurve world evaluation exceeds tolerance";return out;}
        }
    }
    return out;
}
int Run(const std::filesystem::path &dir,const std::filesystem::path &csv,unsigned id)
{
    Require(id<320,"Invalid ID");auto spec=GeneralExtremeSpecs()[id];
    auto c=ReadGeneralFixture(dir/GeneralExtremeFilename(id));Require(c.id==id,"ID mismatch");
    std::cout<<"BEGIN_INPUT_CERT "<<id<<std::endl;
    double tolerance=spec.sag>0?std::min(1e-6,spec.sag*.1):1e-6;
    auto input=CertifyGeneralExtremePair(c,tolerance);
    GeneralConversionResult result;FixtureCertificate output;std::string status="input_unqualified",message=input.message;
    bool converted=false;
    auto save=[&] {
        std::ofstream file(csv);Require(bool(file),"Cannot write result");file<<std::setprecision(17);
        file<<"id,group,variant,type,rational,qualified,source_bound,status,output_type,output_bound,conversion_status,pcurve_failure,classify_ms,optimize_ms,validate_ms,total_ms,geometry_bound,boundary_bound,attempts,iterations,validation_pieces,chart_attempts,chart_iterations,seed_groups,raw_candidates,duplicate_candidates,full_scores,torus_groups,torus_rank_deficient,torus_no_real_axis,torus_candidates,torus_dual_groups,message\n";
        auto &t=result.timing;auto &d=result.classification;auto &tor=d.torus;
        file<<id<<','<<spec.group<<','<<spec.variant<<','<<TypeName(c.truth.type)<<','<<c.rational<<','<<input.passed<<','<<input.upper<<','<<status<<','<<(result.patch?TypeName(result.patch->Parameters().type):"Unknown")<<','<<output.upper<<','<<unsigned(result.status)<<','<<PCurveFailureName(result.validation.pcurveFailure)<<',';
        if(converted)file<<t.classification.elapsedMs<<','<<t.optimization.elapsedMs<<','<<t.validation.elapsedMs<<','<<t.totalMs;
        else file<<"nan,nan,nan,nan";
        file<<','<<result.validation.distanceUpperBound<<','<<result.validation.boundaryHausdorffUpperBound<<','<<result.attempts.size()<<','<<result.statistics.optimizationIterations<<','<<result.statistics.validationPieces<<','<<result.chartAttempts<<','<<result.chartIterations<<','<<d.seedAttemptCount<<','<<d.rawCandidateCount<<','<<d.duplicateCandidateCount<<','<<d.fullScoreCount<<','<<tor.sampleGroups<<','<<tor.rankDeficientGroups<<','<<tor.noRealAxisGroups<<','<<tor.candidateCount<<','<<tor.dualCandidateGroups<<','<<Quote(message)<<'\n';
    };
    if(!input.passed){save();return 1;}
    std::cout<<"BEGIN_CONVERT "<<id<<std::endl;
    result=ConvertGeneral(*c.input);converted=true;message=result.message;
    SeedDiagnostics(csv,result);
    std::cout<<"END_CONVERT "<<id<<std::endl;
    bool success=result.status==ConversionStatus::Success;
    status=success?"audit_pending":"conversion_failed";save();
    if(!success){Require(!result.patch,"Failed conversion published a patch");return 1;}
    std::cout<<"BEGIN_OUTPUT_CERT "<<id<<std::endl;
    output=Audit(c,result);
    status=output.passed?(c.truth.type==result.patch->Parameters().type?"same_type_success":"alternate_type_success"):"output_audit_unresolved";
    if(!output.passed)message+="; "+output.message;
    save();std::cout<<"END_OUTPUT_CERT "<<id<<' '<<status<<std::endl;
    return output.passed?0:1;
}
void CertificateChecks(const std::filesystem::path &dir)
{
    for(unsigned id:{128u,144u,272u,278u})
    {
        auto c=ReadGeneralFixture(dir/GeneralExtremeFilename(id));
        Require(CertifyGeneralExtremePair(c,1e-6).passed,"New fixture certificate failed");
        if(c.truth.type==sggk::SurfaceType::Plane)
        {
            c.truth.origin=Point(Get(c.truth.origin)+Get(c.truth.axis)*.001);
            Require(!CertifyGeneralExtremePair(c,1e-5).passed,"Certificate accepted an offset plane");
        }
        else
        {
            c.corners.fill(c.corners[0]);
            Require(!CertifyGeneralExtremePair(c,1e-5).passed,"Certificate accepted a point chart");
        }
    }
    std::cout<<"Certificate checks passed on new triangle, pole and translated-plane fixtures"<<std::endl;
}
void SphereTriangleChecks(const std::filesystem::path &dir)
{
    // Only inputs from the sphere_triangle group participate in these checks.
    for(unsigned id:{128u,144u,160u})
    {
        auto c=ReadGeneralFixture(dir/GeneralExtremeFilename(id));
        auto cells=ExtractCells(*c.input);
        GeneralTrim trim;trim.parameters=c.truth;trim.corners=c.corners;
        GeneralValidationReport report;
        Require(TryPairedTrimBound(*c.input,cells,trim,GeneralConversionOptions{}.fitting.validation,report),
            "Paired certificate rejected a valid sphere triangle/pole chart");
        Require(report.Passed(),"Paired sphere certificate did not pass all checks");
        using UV=sggk::Point2D;
        const double tiny=std::ldexp(1.0,-40);
        const std::array<UV,4> invalid[]={
            {UV{0,0},UV{0,0},UV{0,0},UV{0,0}}, // Entire domain is a point.
            {UV{0,0},UV{1,1},UV{2,2},UV{2,2}}, // Collapsed side, collinear vertices.
            {UV{0,0},UV{1,0},UV{1,1},UV{0,1}}, // Folded quadrilateral.
            {UV{0,0},UV{1,0},UV{0,1},UV{0,0}}, // Nonadjacent duplicate corners.
            {UV{0,0},UV{1,0},UV{tiny,1},UV{0,1}} // Near collapse, with a fold.
        };
        for(const auto &corners:invalid)
        {
            trim.corners=corners;report={};
            Require(!TryPairedTrimBound(*c.input,cells,trim,GeneralConversionOptions{}.fitting.validation,report),
                "Paired certificate accepted an invalid UV domain");
            Require(report.piecesVisited==0,"Invalid UV domain entered geometric subdivision");
        }
    }
    std::cout<<"Sphere triangle checks passed: 3 valid domains and 15 invalid-domain checks"<<std::endl;
}
void RadiusArcChecks(const std::filesystem::path &dir)
{
    unsigned accepted=0,rejected=0;
    for(unsigned id:{8u,24u,38u,44u})
    {
        auto c=ReadGeneralFixture(dir/GeneralExtremeFilename(id));
        auto frame=FrameOf(c.truth);
        std::vector<SurfaceSample> points;
        for(double u:{.1,.5,.9}) for(double v:{.1,.5,.9})
        {
            double a=0,b=0;double w[]={(1-u)*(1-v),u*(1-v),(1-u)*v,u*v};
            for(unsigned k=0;k<4;++k) {a+=w[k]*c.corners[k].X();b+=w[k]*c.corners[k].Y();}
            SurfaceSample p;p.point=frame.global(CanonicalPoint(c.truth,a,b));
            p.normal=frame.vector(CanonicalNormal(c.truth,a,b));p.normalValid=true;
            points.push_back(p);
        }
        double scale=0;for(const auto &p:points) scale=std::max(scale,Norm(p.point-points.front().point));
        for(unsigned flip=0;flip<2;++flip)
        {
            if(flip) for(auto &p:points) p.normal=p.normal*-1;
            GeneralSeedGroupDiagnostic d;
            auto seed=ConditionedGeneralSeed(c.truth.type,points,scale,d);
            // Held-out positions, including ends: never used in constructing this seed.
            for(double u:{0.,.25,.75,1.}) for(double v:{0.,.25,.75,1.})
            {
                double a=0,b=0;double w[]={(1-u)*(1-v),u*(1-v),(1-u)*v,u*v};
                for(unsigned k=0;k<4;++k){a+=w[k]*c.corners[k].X();b+=w[k]*c.corners[k].Y();}
                V point=frame.global(CanonicalPoint(c.truth,a,b));
                Require(Norm(point-ProjectToMother(seed,point))<1e-5,"Conditioned short-arc seed misses held-out truth positions");
            }
            ++accepted;
        }
        for(unsigned mode=0;mode<2;++mode)
        {
            auto invalid=points;
            for(auto &p:invalid)
                if(mode==0) p.normal=points.front().normal;
                else p=points.front();
            GeneralSeedGroupDiagnostic d;bool failed=false;
            try {ConditionedGeneralSeed(c.truth.type,invalid,scale,d);}
            catch(const Failure &) {failed=true;}
            Require(failed && d.failure!=GeneralSeedFailure::None,"Degenerate seed group was not diagnosed");
            ++rejected;
        }
        SurfaceSampleContext context(*c.input);GeneralClassificationDiagnostics d;
        auto pool=InitializeGeneralCandidates(context,GeneralConversionOptions{},d);
        Require(d.seedAttemptCount==40 && d.seedGroups.size()==40,"Changed eight-group construction budget");
        unsigned multi=0;
        for(const auto &g:d.seedGroups) if(g.type==c.truth.type)
        {
            multi+=g.multiPoint;
            if(g.multiPoint) Require(g.generated==1 && g.sampleCount<=12,"Missing bounded multi-point seed");
        }
        Require(multi==1,"Repeated multi-point construction");
    }
    std::cout<<"Radius arc checks passed: "<<accepted<<" oriented seeds, "<<rejected
             <<" degenerate groups, four fixed-budget sample pools"<<std::endl;
}
}
int main(int argc,char **argv)
{
    std::cerr << "Retired analytic-UV fallback dataset entry. Use GeneralRationalDatasetTests and rational-domain-v2.\n";
    return 1;

    sggk::init();int code=2;
    try
    {
        Require(argc>=3,"Usage: --generate dir [first last] | --run dir csv id | --checks dir | --sphere-triangle-checks dir | --radius-arc-checks dir");
        if(std::string(argv[1])=="--generate")
        {GenerateGeneralExtreme(argv[2],argc>3?unsigned(std::stoul(argv[3])):0,argc>4?unsigned(std::stoul(argv[4])):320);code=0;}
        else if(std::string(argv[1])=="--checks"){CertificateChecks(argv[2]);code=0;}
        else if(std::string(argv[1])=="--sphere-triangle-checks"){SphereTriangleChecks(argv[2]);code=0;}
        else if(std::string(argv[1])=="--radius-arc-checks"){RadiusArcChecks(argv[2]);code=0;}
        else {Require(argc==5&&std::string(argv[1])=="--run","Invalid command");code=Run(argv[2],argv[3],unsigned(std::stoul(argv[4])));}
    }
    catch(const std::exception &e){std::cerr<<"GENERAL EXTREME ERROR: "<<e.what()<<std::endl;}
    sggk::fini();return code;
}
