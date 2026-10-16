#include "GeneralRationalDataset.h"
#include "GeneralExtremeDataset.h"
#include "ApproximateInternal.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
void Require(bool ok,const std::string &s){if(!ok)throw std::runtime_error(s);}
void RequireRepresentation(const GeneralFixture &c)
{
    const auto &s=*c.input;
    Require(s.DegreeU()>=3&&s.DegreeU()<=5&&s.DegreeV()>=3&&s.DegreeV()<=5,
            "V2 inputs must have degree 3 through 5 in both directions");
    Require(s.ControlPoints().size()*s.ControlPoints().front().size()<=200,"V2 input exceeds 200 controls");
    Require((c.extra&&c.group==GeneralFixtureGroup::Tolerance)||c.pairTolerance<=1e-5,
            "Positive construction tolerance exceeds 1e-5");
}
double Mid(I a){Require(a.finite(),"Nonfinite fixture coefficient");return a.lo+(a.hi-a.lo)*.5;}
BPoly Elevate(const BPoly &p,unsigned u,unsigned v){Require(p.u<=u&&p.v<=v,"Fixture degree too small");return p*BPoly(u-p.u,v-p.v,I(1));}
double B(unsigned n,unsigned i,double x)
{double c=1;for(unsigned k=1;k<=i;++k)c*=double(n+1-k)/k;return c*std::pow(x,int(i))*std::pow(1-x,int(n-i));}
std::vector<std::vector<double>> Matrix(unsigned n)
{
    std::vector<std::vector<double>> a(n+1,std::vector<double>(n+1));
    for(unsigned i=0;i<=n;++i)for(unsigned j=0;j<=n;++j)a[i][j]=B(n,j,(1-std::cos(Pi*i/n))*.5);return a;
}
std::vector<V> Solve(const std::vector<std::vector<double>> &a,const std::vector<V> &values)
{
    std::vector<V> out(a.front().size());for(unsigned d=0;d<3;++d)
    {
        std::vector<double> b,x;for(auto p:values)b.push_back(d==0?p.x:d==1?p.y:p.z);
        Require(SolveLeastSquares(a,b,x),"Rational-domain fixture interpolation rank deficiency");
        for(size_t i=0;i<x.size();++i)(d==0?out[i].x:d==1?out[i].y:out[i].z)=x[i];
    }return out;
}
void Snap(GeneralFixture &c,unsigned mask,const RationalReference &ref)
{
    unsigned nu=unsigned(c.input->ControlPoints().size()),nv=unsigned(c.input->ControlPoints()[0].size());
    for(unsigned e=0;e<4;++e)if(mask&(1u<<e))
    {
        auto p=Point(ref.Value(e<2?double(e):0,e<2?0:double(e-2)));
        for(unsigned j=0;j<(e<2?nv:nu);++j)c.input->SetControlPoint(p,e<2?(e?nu-1:0):j,e<2?j:(e==2?0:nv-1));
    }
    c.collapsedMask=FixtureCollapsedMask(*c.input);
}
struct FitAxis
{
    unsigned degree;std::vector<double> knots,flat,points;sggk::UIntArray mults;
    FitAxis(unsigned p,const std::vector<double> &breaks,const std::vector<double> &seams):degree(p),knots(breaks)
    {
        for(size_t k=0;k<knots.size();++k)
        {
            unsigned m=(k==0||k+1==knots.size())?p+1:
                std::find(seams.begin(),seams.end(),knots[k])!=seams.end()?p-1:1;
            mults.push_back(m);flat.insert(flat.end(),m,knots[k]);
        }
        size_t n=flat.size()-p-1;points.resize(n);
        for(size_t i=0;i<n;++i){double sum=0;for(unsigned j=1;j<=p;++j)sum+=flat[i+j];points[i]=sum/p;}
        points.front()=0;points.back()=1;
    }
    std::vector<double> BasisAt(double t) const
    {
        std::vector<double> a(flat.size()-1);
        if(t==1){a.resize(points.size());a.back()=1;return a;}
        for(size_t i=0;i<a.size();++i)a[i]=flat[i]<=t&&t<flat[i+1]?1.:0.;
        for(unsigned d=1;d<=degree;++d)
        {
            std::vector<double> b(a.size()-1);
            for(size_t i=0;i<b.size();++i)
            {
                if(flat[i+d]>flat[i])b[i]+=(t-flat[i])/(flat[i+d]-flat[i])*a[i];
                if(flat[i+d+1]>flat[i+1])b[i]+=(flat[i+d+1]-t)/(flat[i+d+1]-flat[i+1])*a[i+1];
            }
            a=std::move(b);
        }
        return a;
    }
    std::vector<double> Samples() const
    {
        std::vector<double> out=points;
        for(size_t k=0;k+1<knots.size();++k)for(unsigned i=0;i<=12;++i)
            out.push_back(knots[k]+(knots[k+1]-knots[k])*(1-std::cos(Pi*i/12))*.5);
        std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());return out;
    }
    std::vector<V> FitSamples(const std::vector<double> &samples,const std::vector<V> &values) const
    {
        // Fix endpoints exactly, then solve the overdetermined interior fit.
        // Sampling selects coefficients only; publication still requires the
        // independent continuous source certificate.
        std::vector<std::vector<double>> a;std::vector<V> b;
        for(size_t i=1;i+1<samples.size();++i)
        {
            auto row=BasisAt(samples[i]);b.push_back(values[i]-values.front()*row.front()-values.back()*row.back());
            a.emplace_back(row.begin()+1,row.end()-1);
        }
        auto cp=Solve(a,b);cp.insert(cp.begin(),values.front());cp.push_back(values.back());return cp;
    }
};
sggk::BSplineSurfacePtr FitPolynomial(const GeneralFixture &f,unsigned p,unsigned q,
    const std::array<std::vector<double>,2> &breaks,const sggk::BSplineSurface *probe=nullptr)
{
    auto local=f;local.truth.origin={0,0,0};RationalReference ref(local);auto seams=ref.TensorBreaks();
    FitAxis u(p,breaks[0],seams[0]),v(q,breaks[1],seams[1]);size_t nu=u.points.size(),nv=v.points.size();
    Require(nu*nv<=200,"Degree-five spline control budget exhausted");
    auto value=[&](double s,double t){return probe?Get(probe->CalcPoint(s,t))-Get(f.truth.origin):ref.Value(s,t);};
    auto us=u.Samples(),vs=v.Samples();V anchor=value(.5,.5);
    std::vector<std::vector<V>> values(nu,std::vector<V>(vs.size()));
    for(size_t j=0;j<vs.size();++j)
    {
        std::vector<V> col;for(double s:us)col.push_back(value(s,vs[j])-anchor);
        col=u.FitSamples(us,col);for(size_t i=0;i<nu;++i)values[i][j]=col[i];
    }
    sggk::Point3DMatrix cp(nu,std::vector<sggk::Point3D>(nv));
    for(size_t i=0;i<nu;++i){auto row=v.FitSamples(vs,values[i]);for(size_t j=0;j<nv;++j)cp[i][j]=Point(Get(f.truth.origin)+anchor+row[j]);}
    return std::make_shared<sggk::BSplineSurface>(p,q,cp,u.knots,v.knots,u.mults,v.mults);
}
BPoly DegreeBounded(const BPoly &poly,unsigned p,unsigned q)
{
    if(poly.u<=p&&poly.v<=q)return Elevate(poly,p,q);
    // Local homogeneous degree reduction; the final stored surface must still
    // pass the independent continuous certificate before it can be published.
    auto au=Matrix(p),av=Matrix(q);std::vector<std::vector<double>> values(p+1,std::vector<double>(q+1));
    for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
        values[i][j]=Mid(Evaluate(poly,(1-std::cos(Pi*i/p))*.5,(1-std::cos(Pi*j/q))*.5).range());
    for(unsigned j=0;j<=q;++j){std::vector<double> b,x;for(unsigned i=0;i<=p;++i)b.push_back(values[i][j]);Require(SolveLeastSquares(au,b,x),"Homogeneous U fit failed");for(unsigned i=0;i<=p;++i)values[i][j]=x[i];}
    BPoly out(p,q);
    for(unsigned i=0;i<=p;++i){std::vector<double> x;Require(SolveLeastSquares(av,values[i],x),"Homogeneous V fit failed");for(unsigned j=0;j<=q;++j)out.at(i,j)=I(x[j]);}
    return out;
}
sggk::BSplineSurfacePtr Fit(const GeneralFixture &f,unsigned p,unsigned q,
                           std::array<std::vector<double>,2> breaks)
{
    if(!f.rational)return FitPolynomial(f,p,q,breaks);
    auto local=f;local.truth.origin={0,0,0};RationalReference ref(local);
    V anchor=ref.Value(.5,.5);unsigned us=unsigned(breaks[0].size()-1),vs=unsigned(breaks[1].size()-1);
    std::vector<BezierCell> targets;
    if(f.rational)
    {
        for(unsigned i=0;i<us;++i)for(unsigned j=0;j<vs;++j)
        {
            auto candidates=ref.On(breaks[0][i],breaks[0][i+1],breaks[1][j],breaks[1][j+1],anchor);
            V actual=ref.Value((breaks[0][i]+breaks[0][i+1])*.5,(breaks[1][j]+breaks[1][j+1])*.5)-anchor;
            size_t best=0;double error=INFINITY;
            for(size_t k=0;k<candidates.size();++k)
            {
                auto &c=candidates[k];auto n=Range(Evaluate(c.p,.5,.5));I w=Evaluate(c.w,.5,.5).range();
                double d=Norm(V{Mid(n.x/w),Mid(n.y/w),Mid(n.z/w)}-actual);if(d<error){error=d;best=k;}
            }
            auto c=candidates[best];
            if(f.type==0)
            {
                BPoly factor(1,1);for(unsigned a=0;a<2;++a)for(unsigned b=0;b<2;++b)
                    factor.at(a,b)=I(1+.2*breaks[0][i+a])*I(1-.15*breaks[1][j+b]);
                c.p=c.p*factor;c.w=c.w*factor;
            }
            p=std::min(5u,std::max({p,c.p.x.u,c.p.y.u,c.p.z.u,c.w.u}));q=std::min(5u,std::max({q,c.p.x.v,c.p.y.v,c.p.z.v,c.w.v}));targets.push_back(std::move(c));
        }
    }
    unsigned nu=us*p+1,nv=vs*q+1;Require(nu*nv<=200,"Reconstructed input exceeds 200 controls");
    sggk::Point3DMatrix cp(nu,std::vector<sggk::Point3D>(nv));sggk::RealMatrix weights(nu,std::vector<double>(nv,1));
    auto au=Matrix(p),av=Matrix(q);
    for(unsigned su=0;su<us;++su)for(unsigned sv=0;sv<vs;++sv)
    {
        if(f.rational)
        {
            auto &c=targets[su*vs+sv];auto x=DegreeBounded(c.p.x,p,q),y=DegreeBounded(c.p.y,p,q),z=DegreeBounded(c.p.z,p,q),w=DegreeBounded(c.w,p,q);
            for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
            {
                double d=Mid(w.at(i,j));Require(d>0,"Reconstructed coefficient weight nonpositive");
                weights[su*p+i][sv*q+j]=d;
                cp[su*p+i][sv*q+j]=Point(Get(f.truth.origin)+anchor+V{Mid(x.at(i,j))/d,Mid(y.at(i,j))/d,Mid(z.at(i,j))/d});
            }
        }
        else
        {
            std::vector<std::vector<V>> values(p+1,std::vector<V>(q+1));
            for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
            {
                double u=breaks[0][su]+(breaks[0][su+1]-breaks[0][su])*(1-std::cos(Pi*i/p))*.5;
                double v=breaks[1][sv]+(breaks[1][sv+1]-breaks[1][sv])*(1-std::cos(Pi*j/q))*.5;
                values[i][j]=ref.Value(u,v)-anchor;
            }
            for(unsigned j=0;j<=q;++j){std::vector<V> col;for(unsigned i=0;i<=p;++i)col.push_back(values[i][j]);col=Solve(au,col);for(unsigned i=0;i<=p;++i)values[i][j]=col[i];}
            for(unsigned i=0;i<=p;++i){auto row=Solve(av,values[i]);for(unsigned j=0;j<=q;++j)cp[su*p+i][sv*q+j]=Point(Get(f.truth.origin)+anchor+row[j]);}
        }
    }
    auto knots=[](const std::vector<double> &b,unsigned degree){sggk::UIntArray m(b.size(),degree);m.front()=m.back()=degree+1;return m;};
    auto mu=knots(breaks[0],p),mv=knots(breaks[1],q);
    if(f.rational)return std::make_shared<sggk::BSplineSurface>(p,q,cp,weights,breaks[0],breaks[1],mu,mv);
    return std::make_shared<sggk::BSplineSurface>(p,q,cp,breaks[0],breaks[1],mu,mv);
}
void Pad(GeneralFixture &c,unsigned nu,unsigned nv)
{
    // Preserve the original complexity floor where the new fit allows it.
    for(unsigned d=0;d<2;++d)
    {
        unsigned target=d?nv:nu,k=0;
        auto count=[&](){return d?c.input->ControlPoints().front().size():c.input->ControlPoints().size();};
        while(count()<target)
        {
            size_t a=c.input->ControlPoints().size()+(d==0),b=c.input->ControlPoints().front().size()+(d==1);if(a*b>200)break;
            double t=std::fmod(++k*.6180339887498949,1.);if(d)c.input->InsertVKnots(t,1);else c.input->InsertUKnots(t,1);
        }
    }
}
bool RefineFit(const GeneralFixture &c,const RationalReference &ref,unsigned p,unsigned q,
               std::array<std::vector<double>,2> &breaks)
{
    auto seams=ref.TensorBreaks();double best=-1;unsigned chosen=0;double knot=0;
    for(unsigned axis=0;axis<2;++axis)for(size_t k=0;k+1<breaks[axis].size();++k)
    {
        double a=breaks[axis][k],b=breaks[axis][k+1],mid=(a+b)*.5;
        auto next=breaks;next[axis].insert(next[axis].begin()+k+1,mid);
        size_t count=c.rational?(next[0].size()-1)*c.input->DegreeU()+1:FitAxis(p,next[0],seams[0]).points.size();
        count*=c.rational?(next[1].size()-1)*c.input->DegreeV()+1:FitAxis(q,next[1],seams[1]).points.size();
        if(count>200)continue;
        auto error=[&](double t,double other){double u=axis?other:t,v=axis?t:other;return Get(c.input->CalcPoint(u,v))-ref.Value(u,v);};
        double score=0;
        for(double other:{0.,.25,.5,.75,1.})
        {
            V ea=error(a,other),eb=error(b,other);
            for(double t:{.25,.5,.75})score=std::max(score,Norm(error(a+(b-a)*t,other)-(ea*(1-t)+eb*t)));
        }
        if(score>best){best=score;chosen=axis;knot=mid;}
    }
    if(best<0)return false;
    std::cout<<"KNOT "<<c.id<<" axis="<<chosen<<" t="<<knot<<" score="<<best<<" spans="<<breaks[0].size()-1<<','<<breaks[1].size()-1<<std::endl;
    auto &v=breaks[chosen];v.insert(std::lower_bound(v.begin(),v.end(),knot),knot);return true;
}
std::string Expectation(const GeneralFixture &c)
{
    if(c.group==GeneralFixtureGroup::Infeasible)return "reject_infeasible";
    if(c.group==GeneralFixtureGroup::NonlinearPCurve)return "boundary_probe";
    if(c.group==GeneralFixtureGroup::InvalidOptions)return "invalid_options";
    if(c.group==GeneralFixtureGroup::ValidationBudget)return "reject_budget";
    if(c.group==GeneralFixtureGroup::Tolerance)return c.id==49?"sampled_interior_probe":c.id>=45?"reject_known_mother":"accept";
    return "accept";
}
void Manifest(const std::filesystem::path &dir,const std::string &suite,unsigned size)
{
    for(unsigned i=0;i<size;++i)if(!std::filesystem::exists(dir/RationalFixtureName(suite,i)))return;
    std::ofstream out(dir/"manifest.csv");out<<"id,file,contract,group,expectation,type,rational,degree_u,degree_v,controls,collapsed_mask,source_bound,source_tangent_deg\n"<<std::setprecision(17);
    for(unsigned i=0;i<size;++i)
    {
        auto c=ReadGeneralFixture(dir/RationalFixtureName(suite,i));Require(c.rationalDomain,"V1 fixture in V2 manifest");auto &s=*c.input;
        out<<i<<','<<RationalFixtureName(suite,i)<<",rational-domain-v2,"<<FixtureGroupName(c.group)<<','<<Expectation(c)<<','<<c.type<<','<<c.rational<<','<<s.DegreeU()<<','<<s.DegreeV()<<','<<s.ControlPoints().size()*s.ControlPoints().front().size()<<','<<c.collapsedMask<<','<<c.storedUpper<<','<<c.storedTangentDeg<<'\n';
    }
}
}
GeneralFixture RebuildRationalFixture(const GeneralFixture &legacy,bool pressure)
{
    GeneralFixture c=legacy;c.rationalDomain=true;c.storedUpper=-1;c.storedTangentDeg=-1;
    c.chartDomain={INFINITY,-INFINITY,INFINITY,-INFINITY};
    for(auto p:c.corners){c.chartDomain[0]=std::min(c.chartDomain[0],p.X());c.chartDomain[1]=std::max(c.chartDomain[1],p.X());c.chartDomain[2]=std::min(c.chartDomain[2],p.Y());c.chartDomain[3]=std::max(c.chartDomain[3],p.Y());}
    for(auto &p:c.corners)p={(p.X()-c.chartDomain[0])/(c.chartDomain[1]-c.chartDomain[0]),(p.Y()-c.chartDomain[2])/(c.chartDomain[3]-c.chartDomain[2])};
    if(c.group==GeneralFixtureGroup::NonlinearPCurve)
    {
        // Leave room for the prescribed bend INSIDE the rational chart.
        for(auto &p:c.corners)p={.05+.9*p.X(),.05+.85*p.Y()};c.bend=.12;
    }
    unsigned originalType=c.type;if(c.type==5)c.type=0;
    RationalReference ref(c);auto breaks=ref.TensorBreaks();unsigned mask=legacy.collapsedMask;
    unsigned p=3,q=3;
    if(pressure)p=q=5;
    else if(!c.rational&&c.type!=0)
    {
        // Only angularly varying source directions need quintic interpolation.
        for(unsigned axis=0;axis<2;++axis)
        {
            bool curved=false;
            for(unsigned d=0;d<2;++d)if(d==0||c.type==2||c.type==4)
            {
                auto coordinate=[&](unsigned k){return d?c.corners[k].Y():c.corners[k].X();};
                curved|=coordinate(0)!=coordinate(axis?2:1)||coordinate(axis?1:2)!=coordinate(3)||(d&&c.bend&&axis==0);
            }
            (axis?q:p)=curved?5:3;
        }
    }
    unsigned nu=pressure?10:unsigned(legacy.input->ControlPoints().size()),nv=pressure?20:unsigned(legacy.input->ControlPoints().front().size());
    if(c.group==GeneralFixtureGroup::Infeasible){nu=nv=10;}
    if(c.group==GeneralFixtureGroup::Tolerance&&c.id==49)
    {
        // Keep the prescribed interior probe, approximate its source geometry
        // with a quintic/cubic spline instead of retaining the degree-seven net.
        breaks[0]={0,.5,1};c.input=FitPolynomial(c,5,3,breaks,legacy.input.get());c.collapsedMask=0;return c;
    }
    RationalSourceCertificate certificate;unsigned deformationAttempt=0;
    auto refine=[&]()
    {
        if(RefineFit(c,ref,p,q,breaks))return true;
        // Retain tighter tiny-arc tolerances. Ordinary geometry never exceeds
        // the requested 1e-5; intentional offset probes keep their own labels.
        if(c.pairTolerance>=1e-6&&c.pairTolerance<1e-5){c.pairTolerance=1e-5;return true;}
        return false;
    };
    for(unsigned attempt=0;attempt<100;++attempt)
    {
        c.input=Fit(c,p,q,breaks);Pad(c,nu,nv);Snap(c,mask,ref);
        if(c.group==GeneralFixtureGroup::Tolerance)
        {
            V delta=Get(c.truth.axis)*(c.id<45?8e-6:1.2e-5);
            for(unsigned i=0;i<c.input->ControlPoints().size();++i)for(unsigned j=0;j<c.input->ControlPoints()[0].size();++j)c.input->SetControlPoint(Point(Get(c.input->ControlPoints()[i][j])+delta),i,j);
        }
        if(c.group==GeneralFixtureGroup::Infeasible)
        {
            // Certify the undeformed base first: negative labels must not rely
            // on a low-order boundary approximation missing the new contract.
            auto base=c;base.group=GeneralFixtureGroup::Feasible;
            auto baseCertificate=CertifyRationalDomain(base,10000,true);
            if(!baseCertificate.passed)
            {
                std::cout<<"REFIT_NEGATIVE_BASE "<<c.id<<" p="<<p<<" q="<<q<<" reason="<<baseCertificate.message<<std::endl;
                Require(refine(),"Negative reference cannot meet degree/control/tolerance constraints");
                continue;
            }
            // Keep the reference boundary, deform only the interior. Labels
            // depend on a degree-four exclusion proof, never converter results.
            const auto cp=c.input->ControlPoints();double amplitude=std::ldexp(1.,int(deformationAttempt++));
            for(unsigned i=1;i+1<cp.size();++i)for(unsigned j=1;j+1<cp[0].size();++j)
            {
                unsigned h=c.seed+1664525u*(i*101+j);h^=h>>16;double a=(double(h%10007)/5003.5-1)*amplitude;
                V axis{std::sin(i*2.1+j),std::cos(i+j*1.7),1};c.input->SetControlPoint(Point(Get(cp[i][j])+axis*a),i,j);
            }
            certificate=CertifyRationalDomain(c);if(certificate.passed){c.type=originalType;break;}
        }
        else
        {
            certificate=CertifyRationalDomain(c,10000);if(certificate.passed)break;
            std::cout<<"REFIT "<<c.id<<" p="<<p<<" q="<<q<<" reason="<<certificate.message<<std::endl;
            if(!refine())
            {
                double sampled=0;for(unsigned i=0;i<=64;++i)for(unsigned j=0;j<=64;++j)sampled=std::max(sampled,Norm(Get(c.input->CalcPoint(i/64.,j/64.))-ref.Value(i/64.,j/64.)));
                std::cout<<"FIT_LIMIT "<<c.id<<" sampled="<<sampled<<" controls="<<c.input->ControlPoints().size()<<','<<c.input->ControlPoints().front().size()<<std::endl;break;
            }
        }
    }
    Require(certificate.passed,"Reconstruction certificate failed id="+std::to_string(c.id)+": "+certificate.message);
    c.storedUpper=c.group==GeneralFixtureGroup::Infeasible?-1:certificate.upper;c.storedTangentDeg=certificate.tangentDeg;c.storedRankMargin=certificate.rankMargin;
    return c;
}
std::string RationalFixtureName(const std::string &suite,unsigned id)
{std::ostringstream s;s<<suite<<'-'<<std::setfill('0')<<std::setw(3)<<id<<".b2a";return s.str();}
void GenerateRationalDataset(const std::filesystem::path &main,const std::filesystem::path &extreme,
                            const std::filesystem::path &root,const std::string &suite,unsigned first,unsigned last)
{
    unsigned size=suite=="main"?600:suite=="extra"?60:suite=="extreme"?320:suite=="synthetic"?5:0;
    Require(first<last&&last<=size,"Invalid V2 generation suite/range");auto dir=root/suite;std::filesystem::create_directories(dir);
    auto specs=GeneralExtremeSpecs();
    for(unsigned id=first;id<last;++id)
    {
        auto path=dir/RationalFixtureName(suite,id);Require(!std::filesystem::exists(path),"Refusing to overwrite V2 input "+path.string());
        GeneralFixture old;
        if(suite=="synthetic")old=ReadGeneralFixture(main/GeneralFixtureFilename(id*100+4));
        else old=ReadGeneralFixture((suite=="extreme"?extreme:main)/(suite=="extreme"?GeneralExtremeFilename(id):GeneralFixtureFilename(id,suite=="extra")));
        if(suite=="synthetic"){old.id=id;old.rational=true;old.extra=false;old.group=GeneralFixtureGroup::Feasible;}
        if(suite=="extreme"&&specs[id].sag>0)old.pairTolerance=std::min(1e-6,specs[id].sag*.1);
        auto c=RebuildRationalFixture(old,suite=="synthetic");RequireRepresentation(c);
        for(const auto &row:c.input->ControlPoints())for(auto p:row)for(double x:{p.X(),p.Y(),p.Z()})Require(std::isfinite(x)&&std::abs(x)<=1e6,"Input control hull exceeds world range");
        // Write into a staging name. Only fully re-read, certified files publish.
        auto temp=path;temp+=".pending";WriteGeneralFixture(c,temp);auto loaded=ReadGeneralFixture(temp);auto cert=CertifyRationalDomain(loaded);
        Require(cert.passed,"Disk certificate failed: "+cert.message);loaded.storedUpper=loaded.group==GeneralFixtureGroup::Infeasible?-1:cert.upper;loaded.storedTangentDeg=cert.tangentDeg;
        WriteGeneralFixture(loaded,temp);std::filesystem::rename(temp,path);
        std::cout<<"GENERATED "<<suite<<':'<<id<<" source="<<cert.upper<<" tangent="<<cert.tangentDeg<<" pieces="<<cert.pieces<<std::endl;
    }
    Manifest(dir,suite,size);
}
void AuditRationalDataset(const std::filesystem::path &root,const std::string &suite,unsigned first,unsigned last)
{
    auto target=root/suite/("audit-"+std::to_string(first)+"-"+std::to_string(last)+".csv");auto pending=target;pending+=".pending";
    std::ofstream out(pending);Require(bool(out),"Cannot write source audit");
    out<<"id,source_certified,source_bound,source_tangent_deg,rank_margin,pieces,sampled_pair_error,edge0_length_sampled,edge1_length_sampled,edge2_length_sampled,edge3_length_sampled,opposite_u_mid_distance,opposite_v_mid_distance,max_edge_sag_sampled,control_world_absmax,boundary_source_bound\n"<<std::setprecision(17);
    for(unsigned id=first;id<last;++id)
    {
        auto c=ReadGeneralFixture(root/suite/RationalFixtureName(suite,id));Require(c.rationalDomain&&c.id==id,"V2 identity mismatch");RequireRepresentation(c);auto cert=CertifyRationalDomain(c);
        Require(cert.passed,"V2 source audit failed "+suite+':'+std::to_string(id)+": "+cert.message);
        auto point=[&](double u,double v){return Get(c.input->CalcPoint(u,v));};double world=0,sampled=0,lengths[4]{},sag=0;
        for(auto &row:c.input->ControlPoints())for(auto p:row)world=std::max({world,std::abs(p.X()),std::abs(p.Y()),std::abs(p.Z())});
        if(c.group!=GeneralFixtureGroup::Infeasible)
        {
            RationalReference ref(c);
            for(unsigned i=0;i<=16;++i)for(unsigned j=0;j<=16;++j)sampled=std::max(sampled,Norm(point(i/16.,j/16.)-ref.Value(i/16.,j/16.)));
            Require(sampled<=cert.upper+64*std::numeric_limits<double>::epsilon()*(1+world),"SDK samples contradict continuous source bound");
        }
        for(unsigned e=0;e<4;++e)
        {
            auto edge=[&](double t){return e<2?point(double(e),t):point(t,double(e-2));};
            V a=edge(0),b=edge(1),previous=a;
            for(unsigned i=1;i<=128;++i){double t=i/128.;V p=edge(t);lengths[e]+=Norm(p-previous);previous=p;sag=std::max(sag,Norm(p-(a*(1-t)+b*t)));}
        }
        out<<id<<",1,"<<(c.group==GeneralFixtureGroup::Infeasible?-1:cert.upper)<<','<<cert.tangentDeg<<','<<cert.rankMargin<<','<<cert.pieces<<','<<sampled;
        for(double x:lengths)out<<','<<x;
        out<<','<<Norm(point(1,.5)-point(0,.5))<<','<<Norm(point(.5,1)-point(.5,0))<<','<<sag<<','<<world<<','<<cert.boundaryUpper<<'\n';
        std::cout<<"AUDITED "<<suite<<':'<<id<<std::endl;
    }
    out.close();std::filesystem::rename(pending,target);
    unsigned size=suite=="main"?600:suite=="extra"?60:suite=="extreme"?320:suite=="synthetic"?5:0;
    if(first==0&&last==size)Manifest(root/suite,suite,size);
}
}
