// Documentation evidence: three selected cases, always through the sequential API.
// Source certification, independent output audits and plotting exports are untimed.
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
namespace {
using Row=std::map<std::string,std::string>;
void Require(bool b,const std::string &m){if(!b)throw std::runtime_error(m);}
template<class T> void Set(Row &r,const std::string &k,T v)
{std::ostringstream s;s<<std::setprecision(17)<<v;r[k]=s.str();}
#include "SequentialOutputAudit.h"
void Export(const sggk::BSplineSurface &source,const GeneralBoundaryPatch &p,const std::filesystem::path &file)
{
    std::ofstream out(file);Require(bool(out),"Cannot write geometry");
    out<<std::setprecision(17)<<"kind,edge,i,j,x,y,z,mother_distance\n";
    auto frame=FrameOf(p.mother);
    auto write=[&](const char *kind,int edge,unsigned i,unsigned j,V point,double d){
        auto local=frame.local(point);out<<kind<<','<<edge<<','<<i<<','<<j<<','<<local.x<<','<<local.y<<','<<local.z<<','<<d<<'\n';};
    SurfaceSampleContext samples(source);
    rb::Trim trim;trim.mother=p.mother;trim.domain=p.chartDomain;trim.polygon.corners=p.rationalCorners;
    AnalyticPatch ap;ap.parameters=p.mother;ap.domain=p.chartDomain;
    trim.chart=paired::BuildChart(ap,false);paired::PlaceChart(trim.chart,p.mother,Get(p.mother.origin));
    constexpr unsigned n=50;
    for(unsigned i=0;i<=n;++i)for(unsigned j=0;j<=n;++j)
    {
        double u=double(i)/n,v=double(j)/n;
        auto point=samples.At(source.MinParamU()+(source.MaxParamU()-source.MinParamU())*u,
            source.MinParamV()+(source.MaxParamV()-source.MinParamV())*v,sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
        write("source",-1,i,j,point,MotherDistance(p.mother,point));
        const auto &c=p.rationalCorners;
        double x=(1-u)*(1-v)*c[0].X()+u*(1-v)*c[1].X()+(1-u)*v*c[2].X()+u*v*c[3].X();
        double y=(1-u)*(1-v)*c[0].Y()+u*(1-v)*c[1].Y()+(1-u)*v*c[2].Y()+u*v*c[3].Y();
        write("target",-1,i,j,rb::ChartPoint(trim,{x,y}),0);
    }
    for(unsigned e=0;e<4;++e)for(unsigned i=0;i<=200;++i)
    {
        double f=double(i)/200;
        double u=e<2?(e?source.MaxParamU():source.MinParamU()):source.MinParamU()+(source.MaxParamU()-source.MinParamU())*f;
        double v=e<2?source.MinParamV()+(source.MaxParamV()-source.MinParamV())*f:(e==2?source.MinParamV():source.MaxParamV());
        write("source_edge",int(e),i,0,samples.At(u,v,sggk::DerivType::NONE,sggk::DerivType::NONE,false).point,0);
        write("output_edge",int(e),i,0,Get(p.boundaries[e]->CalcPoint(f)),0);
    }
    for(unsigned e=0;e<4;++e)if(p.tangentNotApplicable[e])
        write("collapsed",int(e),0,0,Get(p.vertices[rb::Ends[e][0]]),0);
}
int Run(const std::filesystem::path &data,const std::filesystem::path &output)
{
    std::filesystem::create_directories(output);
    std::ofstream runs(output/"runs.csv");Require(bool(runs),"Cannot write runs");
    runs<<std::setprecision(17)<<"case,run,type,degree_u,degree_v,controls,source_bound,iso_attempted,fallback_attempted,path,accepted,audit,total_ms,iso_ms,fallback_ms,under_10ms,within_50ms,boundary_bound,tangent_bound_deg,output_mother_bound,audit_edge_max,audit_position_max,collapsed_edges,micro_exemptions,optimization_iterations\n";
    bool all=true;
    for(const std::string name:{"sphere-317","cone-217","torus-triangle"})
    {
        sggk::BSplineSurfacePtr input;sggk::SurfaceType truth;double sourceBound=0;
        if(name!="torus-triangle")
        {
            unsigned id=name=="sphere-317"?317:217;
            auto c=ReadDatasetCase(data/"dataset500"/CaseFilename(id));CheckDataset500Contract(c,id);
            auto cert=BoundPair(*c.source,*c.input,c.refitted?1e-6:1e-9);
            Require(cert.passed,"Source certificate failed");sourceBound=cert.upper;input=c.input;truth=c.truth.type;
        }
        else
        {
            auto c=ReadGeneralFixture(data/"general-rational-v2/main/main-400.b2a");
            c.id=1000;c.rational=true;c.group=GeneralFixtureGroup::Feasible;c.extra=false;c.bend=0;c.collapsedMask=8;
            // Source order: 00,10,01,11. Vmax collapses to a single vertex.
            // These angular corners specify the independent fixture generator's chart.
            c.corners={sggk::Point2D{.2,.2},sggk::Point2D{1.4,.2},sggk::Point2D{.2,1.1},sggk::Point2D{.2,1.1}};
            c=RebuildRationalFixture(c);
            Require(c.collapsedMask==8,"Torus edge did not collapse");
            auto cert=CertifyRationalDomain(c);Require(cert.passed,"Torus source certificate failed");
            WriteGeneralFixture(c,output/"torus-triangle.b2a");
            sourceBound=cert.upper;input=c.input;truth=c.truth.type;
        }
        Require(input->DegreeU()<=5&&input->DegreeV()<=5&&input->ControlPoints().size()*input->ControlPoints()[0].size()<=200,"Representation limit");
        for(unsigned run=0;run<3;++run)
        {
            const auto start=std::chrono::steady_clock::now();
            auto result=ConvertIsoThenFallback(*input);
            double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
            if(ms>50){result.accepted=false;result.timeBudgetExceeded=true;}
            Require(result.isoAttempted,"Iso was skipped");
            bool audited=false;Row values;double edgeMax=INFINITY;unsigned collapsed=0;
            if(result.Passed()&&result.path=="general_rational")
            {
                try {
                    Require(result.fallback.patch->mother.type==truth,"Wrong mother type");
                    edgeMax=AuditGeneral(*input,result.fallback,1e-6);
                    AuditPosition(*input,result.fallback.patch->mother,1e-5,values);audited=true;
                    for(bool b:result.fallback.patch->tangentNotApplicable)collapsed+=b;
                }catch(const std::exception &e){std::cerr<<name<<" audit: "<<e.what()<<'\n';}
                if(run==0&&audited)Export(*input,*result.fallback.patch,output/(name+"-geometry.csv"));
            }
            const auto &r=result.fallback;
            runs<<name<<','<<run<<','<<TypeName(truth)<<','<<input->DegreeU()<<','<<input->DegreeV()<<','<<input->ControlPoints().size()*input->ControlPoints()[0].size()<<','<<sourceBound<<','
                <<result.isoAttempted<<','<<result.fallbackAttempted<<','<<result.path<<','<<result.Passed()<<','<<audited<<','<<ms<<','<<result.isoMs<<','<<result.fallbackMs<<','<<(ms<10)<<','<<(ms<=50)<<','
                <<r.boundaryDistanceUpperBound<<','<<r.boundaryTangentUpperBoundDeg<<','<<r.outputMotherDistanceUpperBound<<','<<edgeMax<<','<<values["audit_position_max"]<<','<<collapsed<<','
                <<(r.patch?r.patch->tangentExemptions.size():0)<<','<<result.optimizationIterations<<'\n';runs.flush();
            std::cout<<name<<" run="<<run<<" accepted="<<result.Passed()<<" audit="<<audited<<" ms="<<ms<<" path="<<result.path<<" message="<<result.message<<'\n';
            all&=result.Passed()&&audited;
        }
    }
    return all?0:1;
}
}
int main(int argc,char **argv)
{
    if(argc!=3){std::cerr<<"data output\n";return 2;}
    sggk::init();int status=2;
    try{status=Run(argv[1],argv[2]);}catch(const std::exception &e){std::cerr<<e.what()<<'\n';}
    sggk::fini();return status;
}
