#include "GeneralDataset.h"
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
void Require(bool ok, const std::string &what) { if (!ok) throw std::runtime_error(what); }
double Random(unsigned seed, unsigned channel)
{
    uint32_t x = seed + 0x9e3779b9u * (channel + 1);
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16;
    return double(x) / 4294967296.0;
}
V World(const AnalyticParameters &a, V p)
{
    V z = Unit(Get(a.axis)), x = Get(a.xDirection);
    x = Unit(x - z * Dot(x, z));
    return Get(a.origin) + x * p.x + Cross(z, x) * p.y + z * p.z;
}
V LocalTarget(const AnalyticParameters &a, double u, double v)
{
    switch (a.type)
    {
    case sggk::SurfaceType::Plane: return {u, v, 0};
    case sggk::SurfaceType::Cylinder: return {a.radius * std::cos(u), a.radius * std::sin(u), v};
    case sggk::SurfaceType::Cone: return {v * std::tan(a.halfAngle) * std::cos(u), v * std::tan(a.halfAngle) * std::sin(u), v};
    case sggk::SurfaceType::Sphere: return {a.radius * std::sin(v) * std::cos(u), a.radius * std::sin(v) * std::sin(u), -a.radius * std::cos(v)};
    case sggk::SurfaceType::Torus: return {(a.majorRadius + a.radius * std::cos(v)) * std::cos(u), (a.majorRadius + a.radius * std::cos(v)) * std::sin(u), a.radius * std::sin(v)};
    default: throw std::runtime_error("Unknown fixture mother");
    }
}
std::array<double, 2> Chart(const GeneralFixture &c, double s, double t)
{
    std::array<double, 2> p{};
    const double w[] = {(1-s)*(1-t), s*(1-t), (1-s)*t, s*t};
    for (unsigned i = 0; i < 4; ++i) { p[0] += w[i]*c.corners[i].X(); p[1] += w[i]*c.corners[i].Y(); }
    p[1] += c.bend*s*(1-s);
    return p;
}
double Bernstein(unsigned n, unsigned i, double t)
{
    double b = 1;
    for (unsigned j = 1; j <= i; ++j) b *= double(n + 1 - j) / j;
    return b * std::pow(t, int(i)) * std::pow(1-t, int(n-i));
}
sggk::BSplineSurfacePtr Interpolate(const GeneralFixture &c, unsigned degree)
{
    const unsigned n = degree + 1;
    const V origin = Get(c.truth.origin);
    double wa = c.rational ? .18 + .2*Random(c.seed, 21) : 0;
    double wb = c.rational ? -.12 - .1*Random(c.seed, 22) : 0;
    std::vector<double> nodes(n);
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    for (unsigned i = 0; i < n; ++i)
    {
        nodes[i] = (1-std::cos(Pi*i/degree))*.5;
        for (unsigned j = 0; j < n; ++j) matrix[i][j] = Bernstein(degree,j,nodes[i]);
    }
    std::vector<std::vector<V>> values(n, std::vector<V>(n));
    for (unsigned i = 0; i < n; ++i) for (unsigned j = 0; j < n; ++j)
        values[i][j] = (FixtureTarget(c,nodes[i],nodes[j])-origin)*((1+wa*nodes[i])*(1+wb*nodes[j]));
    auto solve = [&](const std::vector<V> &v) {
        std::vector<V> out(n);
        for (unsigned d = 0; d < 3; ++d)
        {
            std::vector<double> rhs(n), x;
            for (unsigned i = 0; i < n; ++i) rhs[i] = d == 0 ? v[i].x : d == 1 ? v[i].y : v[i].z;
            Require(SolveLeastSquares(matrix,rhs,x),"Fixture interpolation is rank deficient");
            for (unsigned i = 0; i < n; ++i) (d == 0 ? out[i].x : d == 1 ? out[i].y : out[i].z) = x[i];
        }
        return out;
    };
    for (unsigned j = 0; j < n; ++j)
    {
        std::vector<V> v(n); for (unsigned i = 0; i < n; ++i) v[i] = values[i][j];
        v = solve(v); for (unsigned i = 0; i < n; ++i) values[i][j] = v[i];
    }
    sggk::Point3DMatrix cp(n, std::vector<sggk::Point3D>(n));
    sggk::RealMatrix weights(n, std::vector<double>(n));
    for (unsigned i = 0; i < n; ++i)
    {
        auto row = solve(values[i]);
        for (unsigned j = 0; j < n; ++j)
        {
            weights[i][j] = (1+wa*double(i)/degree)*(1+wb*double(j)/degree);
            cp[i][j] = Point(origin+row[j]/weights[i][j]);
        }
    }
    for (unsigned i : {0u,degree}) for (unsigned j : {0u,degree})
        cp[i][j] = Point(FixtureTarget(c,double(i)/degree,double(j)/degree));
    const sggk::RealArray k{0,1}; const sggk::UIntArray m{n,n};
    if (c.rational) return std::make_shared<sggk::BSplineSurface>(degree,degree,cp,weights,k,k,m,m);
    return std::make_shared<sggk::BSplineSurface>(degree,degree,cp,k,k,m,m);
}
void SnapCollapsed(GeneralFixture &c)
{
    auto &s = *c.input;
    unsigned nu = unsigned(s.ControlPoints().size()), nv = unsigned(s.ControlPoints()[0].size());
    const unsigned ends[4][2] = {{0,2},{1,3},{0,1},{2,3}};
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        auto a = c.corners[ends[edge][0]], b = c.corners[ends[edge][1]];
        bool collapsed = (c.truth.type == sggk::SurfaceType::Sphere && a.Y()==b.Y() &&
            (a.Y()==0 || a.Y()==Pi)) || (c.truth.type == sggk::SurfaceType::Cone && a.Y()==0 && b.Y()==0);
        if (!collapsed || c.bend != 0) continue;
        V local = c.truth.type == sggk::SurfaceType::Cone ? V{} : V{0,0,a.Y()==0 ? -c.truth.radius : c.truth.radius};
        auto p = Point(World(c.truth,local));
        for (unsigned j = 0; j < (edge<2 ? nv : nu); ++j)
            s.SetControlPoint(p,edge<2 ? (edge==0 ? 0 : nu-1) : j,edge<2 ? j : (edge==2 ? 0 : nv-1));
    }
    c.collapsedMask = FixtureCollapsedMask(s);
}
void Pad(GeneralFixture &c)
{
    const std::array<std::array<unsigned,2>,5> counts{{{7,9},{9,11},{10,13},{11,14},{10,20}}};
    unsigned nu=counts[c.bucket][0], nv=counts[c.bucket][1];
    if (c.type==0 && c.bucket==0 && c.variant==0) { nu=5; nv=10; }
    if (c.variant%3==1 && c.bucket<4) --nv;
    for (unsigned axis=0;axis<2;++axis)
    {
        auto size=[&] { return axis ? c.input->ControlPoints()[0].size() : c.input->ControlPoints().size(); };
        unsigned target=axis ? nv : nu, i=1;
        Require(size()<=target,"Initial fit exceeds control budget");
        while(size()<target)
        {
            double t=std::fmod(double(i++)*.6180339887498949,1.0);
            if(c.variant&1) t=std::pow(t,1.3);
            if(axis) c.input->InsertVKnots(t,1); else c.input->InsertUKnots(t,1);
        }
    }
    SnapCollapsed(c);
}
GeneralFixture Positive(unsigned id)
{
    GeneralFixture c;
    c.id=id; c.type=id/100; c.bucket=(id%100)/20; c.variant=(id%20)/2;
    c.rational=id%2==0; c.seed=0x600b2a26u+id/2;
    auto r=[&](unsigned k){return Random(c.seed,k);};
    auto &a=c.truth; a.type=GeneralDatasetType(c.type);
    a.origin=Point({8*(r(0)-.5),8*(r(1)-.5),8*(r(2)-.5)});
    V z=Unit(V{r(3)-.5,r(4)-.5,.4+r(5)});
    a.axis=Direction(z); a.xDirection=Direction(Perpendicular(z));
    a.radius=.8+1.8*r(6); a.majorRadius=a.radius*(1.6+1.8*r(7)); a.halfAngle=.3+.5*r(8);
    double u=.2+3*r(9), v=.5+.9*r(10), du=.18+.1*r(11), dv=.18+.1*r(12);
    if(c.type==0) {u=-.4;v=-.5;du=.6+1.8*r(11);dv=.7+1.2*r(12);}
    if(c.type==1 || c.type==3) {v=.6+.6*r(10);dv=.5+.8*r(12);}
    // Type order is plane, cylinder, sphere, cone, torus.
    if(c.variant>=8 && (c.type==1 || c.type==4)) u=2*Pi-du*.5;
    if(c.variant>=8 && c.type==4) v=2*Pi-dv*.5;
    if(c.variant>=8 && c.type==2) v=c.variant==8 ? 0 : Pi-dv;
    if(c.variant>=8 && c.type==3) v=0;
    c.corners={sggk::Point2D{u,v},sggk::Point2D{u+du,v},sggk::Point2D{u,v+dv},sggk::Point2D{u+du,v+dv}};
    if(c.variant>=2 && c.variant<=7)
    {
        c.corners[1]={u+du,v+dv*(.10+.10*r(13))};
        c.corners[2]={u+du*(.12+.13*r(14)),v+dv};
        c.corners[3]={u+du*(1.05+.18*r(15)),v+dv*(1.08+.14*r(16))};
    }
    // Orient the chart before fitting so the oracle includes the same mapping.
    if(c.seed&1) {std::swap(c.corners[0],c.corners[1]);std::swap(c.corners[2],c.corners[3]);}
    if(c.seed&2) {std::swap(c.corners[0],c.corners[2]);std::swap(c.corners[1],c.corners[3]);}
    if(c.seed&4) std::swap(c.corners[1],c.corners[2]);
    c.input=Interpolate(c,c.type==0 ? 3 : 6); Pad(c);
    return c;
}
GeneralFixture Negative(unsigned index, unsigned attempt)
{
    GeneralFixture c; c.id=500+index; c.group=GeneralFixtureGroup::Infeasible;
    c.type=index<60 ? index/12 : 5; c.rational=index%2==0;
    c.seed=0x6bad0026u+index; c.variant=attempt;
    auto r=[&](unsigned k){return Random(c.seed,k);};
    auto &a=c.truth; a.type=GeneralDatasetType(c.type<5?c.type:0);
    a.origin=Point({4*r(0),4*r(1),4*r(2)}); a.axis=Direction({.1+r(3),.1+r(4),1});
    a.xDirection=Direction(Perpendicular(Get(a.axis)));
    a.radius=8+4*r(5); a.majorRadius=2.4*a.radius; a.halfAngle=.45+.25*r(6);
    const unsigned n=10,degree=3;
    double du=c.type==0 || c.type==5 ? 12 : 1.6, dv=c.type==1 || c.type==3 || c.type==0 || c.type==5 ? 10 : 1.4;
    c.corners={sggk::Point2D{.2,.4},sggk::Point2D{.2+du,.4},sggk::Point2D{.2,.4+dv},sggk::Point2D{.2+du,.4+dv}};
    sggk::RealArray knots{0}; sggk::UIntArray mults{degree+1};
    for(unsigned k=1;k<n-degree;++k) {knots.push_back(double(k)/(n-degree));mults.push_back(1);}
    knots.push_back(1);mults.push_back(degree+1);
    std::vector<double> expanded;for(size_t i=0;i<knots.size();++i) expanded.insert(expanded.end(),mults[i],knots[i]);
    std::vector<double> g(n);for(unsigned i=0;i<n;++i)for(unsigned j=1;j<=degree;++j)g[i]+=expanded[i+j]/degree;
    sggk::Point3DMatrix cp(n,std::vector<sggk::Point3D>(n));
    sggk::RealMatrix w(n,std::vector<double>(n,1));
    double amplitude=(c.type==5?2.0:.18)*std::pow(1.5,attempt);
    for(unsigned i=0;i<n;++i)for(unsigned j=0;j<n;++j)
    {
        auto uv=Chart(c,g[i],g[j]); V p=LocalTarget(a,uv[0],uv[1]);
        V normal;
        if(c.type==0 || c.type==5) normal={0,0,1};
        else if(c.type==1) normal={std::cos(uv[0]),std::sin(uv[0]),0};
        else if(c.type==2) normal=Unit(p);
        else if(c.type==3) normal={std::cos(a.halfAngle)*std::cos(uv[0]),std::cos(a.halfAngle)*std::sin(uv[0]),-std::sin(a.halfAngle)};
        else normal={std::cos(uv[1])*std::cos(uv[0]),std::cos(uv[1])*std::sin(uv[0]),std::sin(uv[1])};
        double bump=2*r(100+i*n+j)-1;
        if(index%3==0 && (i==0 || j==0 || i==n-1 || j==n-1)) bump=0;
        p=p+normal*(amplitude*bump);
        if(c.type==5) {p.x+=.3*std::sin(2*p.y);p.y+=.3*std::cos(p.x);}
        cp[i][j]=Point(World(a,p));
        if(c.rational)w[i][j]=.9+.2*r(300+i*n+j);
    }
    if(c.rational)c.input=std::make_shared<sggk::BSplineSurface>(degree,degree,cp,w,knots,knots,mults,mults);
    else c.input=std::make_shared<sggk::BSplineSurface>(degree,degree,cp,knots,knots,mults,mults);
    c.exclusionRadius=2e-5;return c;
}
void MakeSamplingBlindSpot(GeneralFixture &c)
{
    // An interior Bernstein bump peaks between the degree-seven Greville nodes.
    // Its amplitude straddles epsilon: cached samples pass, a separate witness
    // exceeds epsilon. This is a known-plane rejection, not a global negative.
    const unsigned p=7,q=2;
    auto height=[](double u){return Bernstein(7,2,u)-.35*Bernstein(7,4,u);};
    double sampled=std::abs(height(.5)),peak=0;
    for(unsigned i=0;i<=p;++i)sampled=std::max(sampled,std::abs(height(double(i)/p)));
    for(unsigned i=0;i<=10000;++i)peak=std::max(peak,std::abs(height(double(i)/10000)));
    Require(peak>sampled*1.001,"Blind-spot fixture has no sampling gap");
    double amplitude=1e-5/((sampled+peak)*.5);
    c.rational=false;c.collapsedMask=0;
    c.corners={sggk::Point2D{0,0},sggk::Point2D{1,0},sggk::Point2D{0,1},sggk::Point2D{1,1}};
    sggk::Point3DMatrix cp(p+1,std::vector<sggk::Point3D>(q+1));
    for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
    {
        double h=(i==2?1:i==4?-.35:0)*(j==1?2:0)*amplitude;
        cp[i][j]=Point(World(c.truth,{double(i)/p,double(j)/q,h}));
    }
    c.input=std::make_shared<sggk::BSplineSurface>(p,q,cp,sggk::RealArray{0,1},sggk::RealArray{0,1},sggk::UIntArray{p+1,p+1},sggk::UIntArray{q+1,q+1});
}
void WriteV(std::ostream &o,V p){o<<p.x<<' '<<p.y<<' '<<p.z<<'\n';}
V ReadV(std::istream &s){V v;s>>v.x>>v.y>>v.z;return v;}
void Save(const GeneralFixture &c,const std::filesystem::path &path)
{
    std::ofstream o(path);Require(bool(o),"Cannot create fixture");o<<std::setprecision(17);
    o<<"B2A_GENERAL_V1\n"<<c.id<<' '<<unsigned(c.group)<<' '<<c.type<<' '<<c.bucket<<' '<<c.variant<<' '<<c.seed<<' '<<c.rational<<' '<<c.extra<<' '<<c.collapsedMask<<'\n';
    WriteV(o,Get(c.truth.origin));WriteV(o,Get(c.truth.axis));WriteV(o,Get(c.truth.xDirection));
    o<<unsigned(c.truth.type)<<' '<<c.truth.radius<<' '<<c.truth.majorRadius<<' '<<c.truth.halfAngle<<'\n';
    for(const auto &p:c.corners)o<<p.X()<<' '<<p.Y()<<'\n';
    o<<c.bend<<' '<<c.storedUpper<<' '<<c.exclusionRadius<<' '<<c.storedRankMargin<<' '<<c.pairTolerance<<'\n';
    const auto &s=*c.input;const auto &p=s.ControlPoints();const auto &w=s.Weights();
    o<<s.DegreeU()<<' '<<s.DegreeV()<<' '<<p.size()<<' '<<p[0].size()<<' '<<!w.empty()<<'\n';
    auto knots=[&](const auto &k,const auto &m){o<<k.size()<<'\n';for(size_t i=0;i<k.size();++i)o<<k[i]<<' '<<m[i]<<'\n';};
    knots(s.KnotsU(),s.MultsU());knots(s.KnotsV(),s.MultsV());
    for(size_t i=0;i<p.size();++i)for(size_t j=0;j<p[i].size();++j){auto v=Get(p[i][j]);o<<v.x<<' '<<v.y<<' '<<v.z;if(!w.empty())o<<' '<<w[i][j];o<<'\n';}
    Require(bool(o),"Fixture write failed");
}
}
sggk::SurfaceType GeneralDatasetType(unsigned n)
{
    const sggk::SurfaceType t[]={sggk::SurfaceType::Plane,sggk::SurfaceType::Cylinder,sggk::SurfaceType::Sphere,sggk::SurfaceType::Cone,sggk::SurfaceType::Torus};
    if(n>=5)throw std::runtime_error("Invalid fixture type");return t[n];
}
void WriteGeneralFixture(const GeneralFixture &c,const std::filesystem::path &path) { Save(c,path); }
const char *FixtureGroupName(GeneralFixtureGroup g)
{
    const char *names[]={"feasible","infeasible","nonlinear_pcurve","periodic_singular","sphere_priority","tolerance","invalid_options","validation_budget"};
    return unsigned(g)<8?names[unsigned(g)]:"invalid";
}
std::string GeneralFixtureFilename(unsigned id,bool extra)
{
    std::ostringstream s;s<<(extra?"general-extra-":"general-")<<std::setfill('0')<<std::setw(3)<<id<<".b2a";return s.str();
}
V FixtureTarget(const GeneralFixture &c,double u,double v)
{
    auto p=Chart(c,u,v);return World(c.truth,LocalTarget(c.truth,p[0],p[1]));
}
unsigned FixtureCollapsedMask(const sggk::BSplineSurface &s)
{
    unsigned mask=0,nu=unsigned(s.ControlPoints().size()),nv=unsigned(s.ControlPoints()[0].size());
    for(unsigned e=0;e<4;++e)
    {
        V first;bool same=true;
        for(unsigned j=0;j<(e<2?nv:nu);++j)
        {
            auto p=Get(s.ControlPoints()[e<2?(e==0?0:nu-1):j][e<2?j:(e==2?0:nv-1)]);
            if(!j)first=p;else same&=p.x==first.x && p.y==first.y && p.z==first.z;
        }
        if(same)mask|=1u<<e;
    }
    return mask;
}
GeneralFixture ReadGeneralFixture(const std::filesystem::path &path)
{
    std::ifstream in(path);std::string tag;in>>tag;Require(tag=="B2A_GENERAL_V1","Invalid general fixture header: "+path.string());
    GeneralFixture c;unsigned group,type;in>>c.id>>group>>c.type>>c.bucket>>c.variant>>c.seed>>c.rational>>c.extra>>c.collapsedMask;
    Require(group<8 && c.type<=5 && c.bucket<5,"Invalid fixture metadata");c.group=GeneralFixtureGroup(group);
    c.truth.origin=Point(ReadV(in));c.truth.axis=Direction(ReadV(in));c.truth.xDirection=Direction(ReadV(in));
    in>>type>>c.truth.radius>>c.truth.majorRadius>>c.truth.halfAngle;c.truth.type=sggk::SurfaceType(type);
    for(auto &p:c.corners){double u=0,v=0;in>>u>>v;Require(std::isfinite(u)&&std::isfinite(v),"Nonfinite chart corner");p={u,v};}
    in>>c.bend>>c.storedUpper>>c.exclusionRadius>>c.storedRankMargin>>c.pairTolerance;
    unsigned p=0,q=0,nu=0,nv=0;bool weights=false;in>>p>>q>>nu>>nv>>weights;
    Require(!in.fail()&&nu&&nv&&nu<=200&&nv<=200&&nu*nv<=400&&p<nu&&q<nv,"Invalid fixture dimensions");
    sggk::RealArray ku,kv;sggk::UIntArray mu,mv;
    auto knots=[&](auto &k,auto &m,unsigned degree,unsigned count){unsigned n,total=0;in>>n;Require(n>=2&&n<=201,"Invalid knots");k.resize(n);m.resize(n);for(unsigned i=0;i<n;++i){in>>k[i]>>m[i];Require(std::isfinite(k[i])&&m[i]>0&&m[i]<=degree+1&&(!i||k[i]>k[i-1]),"Invalid knot value");total+=m[i];}Require(k.front()==0&&k.back()==1&&total==count+degree+1,"Invalid knot range");};
    knots(ku,mu,p,nu);knots(kv,mv,q,nv);
    sggk::Point3DMatrix cp(nu,std::vector<sggk::Point3D>(nv));sggk::RealMatrix w;
    if(weights)w.assign(nu,std::vector<double>(nv));
    for(unsigned i=0;i<nu;++i)for(unsigned j=0;j<nv;++j){auto v=ReadV(in);Require(Finite(v),"Nonfinite control point");cp[i][j]=Point(v);if(weights){in>>w[i][j];Require(w[i][j]>0&&std::isfinite(w[i][j]),"Invalid positive weight");}}
    Require(!in.fail(),"Truncated fixture");
    if(weights)c.input=std::make_shared<sggk::BSplineSurface>(p,q,cp,w,ku,kv,mu,mv);
    else c.input=std::make_shared<sggk::BSplineSurface>(p,q,cp,ku,kv,mu,mv);
    Require(c.input->IsRational()==c.rational,"Rational representation mismatch");
    Require(FixtureCollapsedMask(*c.input)==c.collapsedMask,"Collapsed-edge metadata mismatch");return c;
}
void GenerateGeneralDataset(const std::filesystem::path &dir,unsigned first,unsigned last)
{
    Require(first<last && last<=660,"Invalid generation range");std::filesystem::create_directories(dir);
    for(unsigned i=first;i<last;++i)
    {
        GeneralFixture c;FixtureCertificate certificate;
        if(i<500)c=Positive(i);
        else if(i<600)
        {
            for(unsigned attempt=0;attempt<10;++attempt)
            {
                c=Negative(i-500,attempt);certificate=CertifyNoQuartic(*c.input,c.exclusionRadius);
                if(certificate.passed)break;
                std::cout<<"negative "<<i<<" construction "<<attempt<<": "<<certificate.message<<std::endl;
            }
            Require(certificate.passed,"No certified negative: "+std::to_string(i));c.storedRankMargin=certificate.rankMargin;
        }
        else
        {
            unsigned j=i-600;
            c=Positive(((j<20?j/4:j%5)*100)+((j%5)*20)+4+(j%2));
            c.extra=true;c.id=j;
            if(j<20){c.group=GeneralFixtureGroup::NonlinearPCurve;c.bend=.12;c.input=Interpolate(c,c.type==0?3:6);Pad(c);}
            else if(j<30){c=Positive((j%5)*100+18+j%2);c.extra=true;c.id=j;c.group=GeneralFixtureGroup::Periodic;}
            else if(j<40){c=Positive((j<35?2:4)*100+40+4+j%2);c.extra=true;c.id=j;c.group=GeneralFixtureGroup::SpherePriority;}
            else if(j<50)
            {
                c=Positive((j%5)*20+4+j%2);c.extra=true;c.id=j;c.group=GeneralFixtureGroup::Tolerance;
                auto cp=c.input->ControlPoints();V offset=Get(c.truth.axis)*(j<45?8e-6:1.2e-5);
                for(unsigned u=0;u<cp.size();++u)for(unsigned v=0;v<cp[u].size();++v)c.input->SetControlPoint(Point(Get(cp[u][v])+offset),u,v);
                c.pairTolerance=2e-5;
                if(j==49)MakeSamplingBlindSpot(c);
            }
            else c.group=j<55?GeneralFixtureGroup::InvalidOptions:GeneralFixtureGroup::ValidationBudget;
        }
        if(c.group!=GeneralFixtureGroup::Infeasible)
        {
            certificate=CertifyAnalyticPair(c,c.pairTolerance);
            Require(certificate.passed,"Positive certificate failed "+std::to_string(i)+": "+certificate.message);
            c.storedUpper=certificate.upper;
        }
        const auto path=dir/GeneralFixtureFilename(c.id,c.extra);Save(c,path);
        auto loaded=ReadGeneralFixture(path);
        auto recheck=c.group==GeneralFixtureGroup::Infeasible?CertifyNoQuartic(*loaded.input,c.exclusionRadius):CertifyAnalyticPair(loaded,c.pairTolerance);
        Require(recheck.passed,"Reloaded certificate failed: "+path.string());
        std::cout<<"Generated "<<i+1<<"/660 "<<FixtureGroupName(c.group)<<" certificate="<<(c.group==GeneralFixtureGroup::Infeasible?recheck.rankMargin:recheck.upper)<<std::endl;
    }
    // Each manifest has its own membership. Historical manifest.csv stays untouched.
    for(bool extra:{false,true})
    {
        bool complete=true;unsigned count=extra?60:600;
        for(unsigned i=0;i<count;++i)complete&=std::filesystem::exists(dir/GeneralFixtureFilename(i,extra));
        if(!complete)continue;
        std::ofstream out(dir/(extra?"manifest-general-extra60.csv":"manifest-general600.csv"));
        out<<"id,file,group,type,rational,controls,degree_u,degree_v,bucket,variant,seed,collapsed_mask,pair_upper,exclusion_radius,rank_margin\n"<<std::setprecision(17);
        for(unsigned i=0;i<count;++i){auto c=ReadGeneralFixture(dir/GeneralFixtureFilename(i,extra));out<<i<<','<<GeneralFixtureFilename(i,extra)<<','<<FixtureGroupName(c.group)<<','<<(c.type<5?TypeName(GeneralDatasetType(c.type)):"Freeform")<<','<<c.rational<<','<<c.input->ControlPoints().size()*c.input->ControlPoints()[0].size()<<','<<c.input->DegreeU()<<','<<c.input->DegreeV()<<','<<c.bucket<<','<<c.variant<<','<<c.seed<<','<<c.collapsedMask<<','<<c.storedUpper<<','<<c.exclusionRadius<<','<<c.storedRankMargin<<'\n';}
        Require(bool(out),"Manifest write failed");
    }
}
}
