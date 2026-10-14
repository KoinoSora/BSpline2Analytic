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
void Require(bool ok,const std::string &message){if(!ok)throw std::runtime_error(message);}
using UV=sggk::Point2D;
GeneralExtremeSpec Base(unsigned type,const std::string &group)
{
    GeneralExtremeSpec s;s.group=group;s.fixture.type=type;s.fixture.truth.type=GeneralDatasetType(type);
    s.fixture.truth.axis=Direction({0,0,1});s.fixture.truth.xDirection=Direction({1,0,0});
    s.fixture.truth.radius=1;s.fixture.truth.halfAngle=.2;
    s.fixture.corners={UV{0,0},UV{.3,0},UV{0,.3},UV{.3,.3}};return s;
}
void Quad(GeneralExtremeSpec &s,double u,double v,double du,double dv,double skew=0,double sideU=0)
{s.fixture.corners={UV{u,v},UV{u+du,v+skew},UV{u+sideU,v+dv},UV{u+du+sideU,v+skew+dv}};}
void Directions(GeneralExtremeSpec &s,unsigned seed)
{
    V z=Unit(V{.23+.017*(seed%7),-.41+.03*(seed%11),.83});
    s.fixture.truth.axis=Direction(z);s.fixture.truth.xDirection=Direction(Perpendicular(z));
    s.fixture.seed=seed;s.orientation=seed%8;
    auto original=s.fixture.corners;
    for(unsigned i=0;i<4;++i)
    {
        unsigned u=i%2,v=i/2;if(s.orientation&1)u=1-u;if(s.orientation&2)v=1-v;
        if(s.orientation&4)std::swap(u,v);s.fixture.corners[i]=original[u+2*v];
    }
}
double AngularSpan(const GeneralFixture &c,unsigned axis)
{
    double span=0;
    for(unsigned k=0;k<2;++k)
    {
        unsigned a=axis?k:2*k,b=a+(axis?2:1);
        if(c.type)span=std::max(span,std::abs(c.corners[a].X()-c.corners[b].X()));
        if(c.type==2||c.type==4)span=std::max(span,std::abs(c.corners[a].Y()-c.corners[b].Y()));
    }
    return span;
}
double Bernstein(unsigned n,unsigned i,double t)
{
    double value=1;for(unsigned j=1;j<=i;++j)value*=double(n+1-j)/j;
    return value*std::pow(t,int(i))*std::pow(1-t,int(n-i));
}
sggk::BSplineSurfacePtr Fit(const GeneralFixture &c,unsigned p,unsigned q,unsigned us,unsigned vs)
{
    auto local=c;local.truth.origin={0,0,0};
    V reference=FixtureTarget(local,.5,.5);
    const unsigned nu=us*p+1,nv=vs*q+1;
    sggk::Point3DMatrix cp(nu,std::vector<sggk::Point3D>(nv));
    sggk::RealMatrix weights(nu,std::vector<double>(nv));
    auto matrix=[](unsigned degree) {
        std::vector<std::vector<double>> a(degree+1,std::vector<double>(degree+1));
        for(unsigned i=0;i<=degree;++i)for(unsigned j=0;j<=degree;++j)
            a[i][j]=Bernstein(degree,j,(1-std::cos(Pi*i/degree))*.5);
        return a;
    };
    auto au=matrix(p),av=matrix(q);
    auto solve=[](const std::vector<std::vector<double>> &a,const std::vector<V> &values) {
        std::vector<V> out(values.size());
        for(unsigned d=0;d<3;++d)
        {
            std::vector<double> rhs,x;for(V v:values)rhs.push_back(d==0?v.x:d==1?v.y:v.z);
            Require(SolveLeastSquares(a,rhs,x),"Interpolation rank deficiency");
            for(size_t i=0;i<x.size();++i)(d==0?out[i].x:d==1?out[i].y:out[i].z)=x[i];
        }
        return out;
    };
    auto weight=[&](double u,double v){return c.rational?(1+.2*u)*(1-.15*v):1.;};
    for(unsigned su=0;su<us;++su)for(unsigned sv=0;sv<vs;++sv)
    {
        std::vector<std::vector<V>> values(p+1,std::vector<V>(q+1));
        for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
        {
            double u=(su+(1-std::cos(Pi*i/p))*.5)/us,v=(sv+(1-std::cos(Pi*j/q))*.5)/vs;
            values[i][j]=(FixtureTarget(local,u,v)-reference)*weight(u,v);
        }
        for(unsigned j=0;j<=q;++j)
        {
            std::vector<V> column;for(unsigned i=0;i<=p;++i)column.push_back(values[i][j]);
            column=solve(au,column);for(unsigned i=0;i<=p;++i)values[i][j]=column[i];
        }
        for(unsigned i=0;i<=p;++i)
        {
            auto row=solve(av,values[i]);
            for(unsigned j=0;j<=q;++j)
            {
                unsigned u=su*p+i,v=sv*q+j;
                weights[u][v]=weight((su+double(i)/p)/us,(sv+double(j)/q)/vs);
                cp[u][v]=Point(reference+row[j]/weights[u][v]);
            }
        }
    }
    for(unsigned i:{0u,nu-1})for(unsigned j:{0u,nv-1})cp[i][j]=Point(FixtureTarget(local,double(i)/(nu-1),double(j)/(nv-1)));
    const unsigned ends[4][2]={{0,2},{1,3},{0,1},{2,3}};
    for(unsigned edge=0;edge<4;++edge)
    {
        auto a=c.corners[ends[edge][0]],b=c.corners[ends[edge][1]];
        bool collapsed=a.X()==b.X()&&a.Y()==b.Y();
        collapsed|=c.type==2&&a.Y()==b.Y()&&(a.Y()==0||a.Y()==Pi);
        collapsed|=c.type==3&&a.Y()==0&&b.Y()==0;
        if(!collapsed)continue;
        unsigned corner=ends[edge][0];auto point=Point(FixtureTarget(local,corner%2,corner/2));
        for(unsigned j=0;j<(edge<2?nv:nu);++j)cp[edge<2?(edge?nu-1:0):j][edge<2?j:(edge==2?0:nv-1)]=point;
    }
    auto knots=[](unsigned pieces,unsigned degree) {
        std::pair<sggk::RealArray,sggk::UIntArray> out;
        for(unsigned i=0;i<=pieces;++i){out.first.push_back(double(i)/pieces);out.second.push_back(i==0||i==pieces?degree+1:degree);}return out;
    };
    auto ku=knots(us,p),kv=knots(vs,q);
    if(c.rational)return std::make_shared<sggk::BSplineSurface>(p,q,cp,weights,ku.first,kv.first,ku.second,kv.second);
    return std::make_shared<sggk::BSplineSurface>(p,q,cp,ku.first,kv.first,ku.second,kv.second);
}
void Pad(GeneralFixture &c)
{
    unsigned i=0;
    while(c.input->ControlPoints().size()*c.input->ControlPoints()[0].size()<50)
    {
        bool alongU=c.input->ControlPoints().size()<=c.input->ControlPoints()[0].size();
        double t=std::fmod(++i*.6180339887498949,1.);
        if(alongU)c.input->InsertUKnots(t,1);else c.input->InsertVKnots(t,1);
    }
    // A nonuniform representation of the same geometry, within the fixed budget.
    size_t nu=c.input->ControlPoints().size(),nv=c.input->ControlPoints()[0].size();
    if(c.seed%3==0&&(nu+1)*nv<=200)c.input->InsertUKnots(1e-6,1);
    if(c.seed%3==1&&nu*(nv+1)<=200)c.input->InsertVKnots(.017,1);
    nu=c.input->ControlPoints().size();nv=c.input->ControlPoints()[0].size();
    const unsigned ends[4][2]={{0,2},{1,3},{0,1},{2,3}};
    for(unsigned edge=0;edge<4;++edge)
    {
        auto a=c.corners[ends[edge][0]],b=c.corners[ends[edge][1]];
        bool collapsed=a.X()==b.X()&&a.Y()==b.Y();
        collapsed|=c.type==2&&a.Y()==b.Y()&&(a.Y()==0||a.Y()==Pi);
        if(!collapsed)continue;
        unsigned corner=ends[edge][0];auto point=Point(FixtureTarget(c,corner%2,corner/2));
        for(unsigned j=0;j<(edge<2?nv:nu);++j)c.input->SetControlPoint(point,unsigned(edge<2?(edge?nu-1:0):j),unsigned(edge<2?j:(edge==2?0:nv-1)));
    }
}
double EdgeLength(const GeneralFixture &c,unsigned edge,unsigned intervals)
{
    auto point=[&](double t){return FixtureTarget(c,edge<2?double(edge):t,edge<2?t:double(edge-2));};
    V previous=point(0);double length=0;
    for(unsigned i=1;i<=intervals;++i){auto next=point(double(i)/intervals);length+=Norm(next-previous);previous=next;}return length;
}
}
std::vector<GeneralExtremeSpec> GeneralExtremeSpecs()
{
    std::vector<GeneralExtremeSpec> geometry;
    auto add=[&](GeneralExtremeSpec s) {Directions(s,unsigned(geometry.size())+31);geometry.push_back(s);};
    for(unsigned type=1;type<5;++type)for(double rho:{100.,1000.})for(unsigned level=0;level<4;++level)
    {
        auto s=Base(type,"radius_arc");s.variant=level?"sag_ratio_"+std::to_string(level):"short_edge_0.010001";
        double shortestRadius=type==2?rho*std::cos(.5/rho):type==3?rho-std::sin(.9):
            type==4?rho*.9+rho*.1*std::cos(5/rho):rho;
        double theta=level?4*std::asin(std::sqrt(std::pow(10.,int(level)-2)*1e-5/(2*rho))):.010001/shortestRadius;
        s.sag=2*rho*std::pow(std::sin(theta/4),2);s.length=rho*theta;s.width=1;
        s.fixture.truth.radius=rho;
        if(type==1)Quad(s,-theta/2,0,theta,1);
        if(type==2){Quad(s,-theta/2,Pi/2-.5/rho,theta,1/rho);}
        if(type==3){s.fixture.truth.halfAngle=.9;double z=rho/std::tan(.9);Quad(s,-theta/2,z-std::cos(.9),theta,std::cos(.9));}
        if(type==4){s.fixture.truth.radius=rho*.1;s.fixture.truth.majorRadius=rho*.9;Quad(s,-theta/2,-.5/s.fixture.truth.radius,theta,1/s.fixture.truth.radius);}
        add(s);
    }
    for(unsigned i=0;i<16;++i)
    {
        auto s=Base(1,"cylinder_helix");s.variant=i%2?"left_seam":"right";
        double radius=std::pow(10.,int(i/4)-1),theta=std::array<double,4>{.3,1.2,3.,6.}[i%4];
        double h=radius*theta*(.25+.25*(i%4)),length=std::hypot(radius*theta,h),width=std::max(.01,.02*length);
        s.fixture.truth.radius=radius;s.length=length;s.width=width;
        Quad(s,i%2?2*Pi+.1:.2,.1,i%2?-theta:theta,width*length/(radius*theta),h);add(s);
    }
    for(unsigned i=0;i<16;++i)
    {
        auto s=Base(1,"slender_cylinder_helix");s.variant="aspect_"+std::to_string(i%4);
        double length=std::array<double,4>{10,100,500,1000}[i/4],theta=std::array<double,4>{.6,1.2,3.,5.8}[i%4];
        double radius=length/(2*theta),h=length*std::sqrt(.75);
        double width=std::max(.01,length/std::pow(10.,2+int(i%4)));
        s.fixture.truth.radius=radius;s.length=length;s.width=width;
        Quad(s,2*Pi-.2,.03,theta,2*width,h);add(s);
    }
    for(unsigned kind=0;kind<3;++kind)for(unsigned i=0;i<8;++i)
    {
        auto s=Base(2,"sphere_triangle");s.variant=kind==0?"regular_uv":kind==1?"pole":"slender";
        double radius=std::array<double,4>{.1,1,100,1000}[i/2];if(kind==2&&radius<1)radius=1;
        s.fixture.truth.radius=radius;
        double u=i%2?2*Pi-.05:.2,v=.8,du=kind==2?.02/(radius*std::sin(v)):.5,dv=radius==1000?.4:.6;
        if(kind==1)
        {
            double pole=i%2?Pi:0;Quad(s,u,pole,.7,i%2?-dv:dv);s.width=radius*std::sin(dv)*.7;
        }
        else {s.fixture.corners={UV{u,v},UV{u+du,v},UV{u+.2*du,v+dv},UV{u+.2*du,v+dv}};s.width=radius*std::sin(v)*du;}
        s.length=radius*dv;add(s);
    }
    for(unsigned i=0;i<16;++i)
    {
        auto s=Base(3,"cone_helix");s.variant=i%2?"left_seam":"right";
        double z=std::array<double,4>{1,5,20,50}[i/4],a=std::array<double,4>{.03,.08,.25,.6}[i%4];
        double theta=std::array<double,4>{.3,1.5,3.,5.8}[i%4],h=.6*z,r=z*std::tan(a),speed=std::hypot(r*theta,h/std::cos(a));
        double width=std::max(.01,.025*speed),delta=width*speed*std::cos(a)/(r*theta);
        s.fixture.truth.halfAngle=a;s.length=speed;s.width=width;
        Quad(s,i%2?2*Pi+.1:.1,z,i%2?-theta:theta,delta,h);add(s);
    }
    for(unsigned i=0;i<16;++i)
    {
        auto s=Base(3,"slender_cone");s.variant=i%4<2?"generator_band":"spiral_band";
        double z=std::array<double,4>{1,5,30,100}[i/4],h=std::array<double,4>{10,50,300,800}[i/4];
        double a=std::array<double,4>{.02,.05,.15,.4}[i%4],width=.0101,turn=i%4<2?0:.3;
        s.fixture.truth.halfAngle=a;s.length=h/std::cos(a);s.width=width;
        Quad(s,2*Pi-.1,z,turn,0,h,width/(z*std::tan(a)));add(s);
    }
    for(unsigned i=0;i<16;++i)
    {
        auto s=Base(4,"torus_extreme");double major=1,minor=.1,du=.15,dv=.2,u=.1,v=.3;
        if(i<4){major=std::array<double,4>{1,10,100,1000}[i];minor=i?.01:.1;dv=std::max(.2,.01/minor);s.variant="radius_ratio";}
        else if(i<8){minor=std::pow(10.,int(i)-5);major=minor+.01;u=-.6;v=Pi-.3;du=1.2;dv=.6;s.variant="inner_neck";}
        else if(i<12){major=std::pow(10.,int(i)-9);minor=.2*major;u=2*Pi-.2;v=2*Pi-.3;du=.4;dv=.6;s.variant="double_seam";}
        else {major=std::pow(10.,int(i)-12);minor=.2*major;du=.01;dv=std::max(.01,.01/minor);s.variant="clustered_normals";}
        s.fixture.truth.majorRadius=major;s.fixture.truth.radius=minor;
        s.length=std::max(.01,(major+minor*std::cos(v+dv/2))*du);s.width=minor*dv;
        Quad(s,u,v,du,dv);add(s);
    }
    Require(geometry.size()==136,"Base matrix size mismatch");
    for(unsigned kind=0;kind<6;++kind)
    {
        GeneralExtremeSpec base;
        if(kind==0){base=Base(0,"world_coordinates");Quad(base,0,0,1000,.01);base.length=1000;base.width=.01;Directions(base,300);}
        else {unsigned from=std::array<unsigned,5>{40,55,66,96,123}[kind-1];base=geometry[from];}
        base.group="world_coordinates";base.variant="geometry_"+std::to_string(kind);
        for(double d:{0.,1e3,1e6,1e9})
        {
            auto s=base;s.translation=d;s.fixture.truth.origin=Point(V{d,-.73*d,.41*d});geometry.push_back(s);
        }
    }
    Require(geometry.size()==160,"Geometry matrix size mismatch");
    std::vector<GeneralExtremeSpec> cases;
    for(auto s:geometry)for(bool rational:{true,false})
    {
        s.fixture.id=unsigned(cases.size());s.fixture.rational=rational;s.fixture.bucket=s.fixture.id%5;
        s.fixture.pairTolerance=1e-6;cases.push_back(s);
    }
    return cases;
}
std::string GeneralExtremeFilename(unsigned id)
{std::ostringstream s;s<<"extreme-"<<std::setfill('0')<<std::setw(3)<<id<<".b2a";return s.str();}
void GenerateGeneralExtreme(const std::filesystem::path &dir,unsigned first,unsigned last)
{
    Require(first<last&&last<=320,"Invalid generation range");auto specs=GeneralExtremeSpecs();
    std::filesystem::create_directories(dir);
    for(unsigned i=first;i<last;++i)
    {
        auto spec=specs[i];auto &c=spec.fixture;auto path=dir/GeneralExtremeFilename(i);
        Require(!std::filesystem::exists(path),"Refusing to overwrite "+path.string());
        auto local=c;local.truth.origin={0,0,0};
        double au=AngularSpan(c,0),av=AngularSpan(c,1);
        unsigned us=std::max(1u,unsigned(std::ceil(au/1.5))),vs=std::max(1u,unsigned(std::ceil(av/1.5)));
        FixtureCertificate certificate;
        const double tolerance=spec.sag>0?std::min(1e-6,spec.sag*.1):1e-6;
        for(unsigned degree:{5u,7u,9u,11u,13u})
        {
            unsigned p=au==0?(c.rational?2u:1u):degree,q=av==0?(c.rational?2u:1u):degree;
            if((us*p+1)*(vs*q+1)>200)continue;
            local.input=Fit(local,p,q,us,vs);Pad(local);
            certificate=CertifyGeneralExtremePair(local,tolerance,3000);
            std::cout<<"FIT "<<i<<" degree="<<p<<','<<q<<" passed="<<certificate.passed<<' '<<certificate.message<<std::endl;
            if(certificate.passed)break;
        }
        Require(bool(local.input),"No representation fits the control budget");
        c.input=local.input;
        V offset=Get(c.truth.origin);
        for(unsigned u=0;u<c.input->ControlPoints().size();++u)for(unsigned v=0;v<c.input->ControlPoints()[0].size();++v)
            c.input->SetControlPoint(Point(Get(c.input->ControlPoints()[u][v])+offset),u,v);
        c.collapsedMask=FixtureCollapsedMask(*c.input);c.storedUpper=certificate.passed?certificate.upper:-1;
        WriteGeneralFixture(c,path);
        auto loaded=ReadGeneralFixture(path);auto disk=CertifyGeneralExtremePair(loaded,spec.sag>0?tolerance:1e-6);
        c.storedUpper=disk.passed?disk.upper:-1;WriteGeneralFixture(c,path);
        std::cout<<"GENERATED "<<i<<" qualified="<<disk.passed<<" upper="<<disk.upper<<' '<<disk.message<<std::endl;
    }
    bool complete=true;for(unsigned i=0;i<320;++i)complete&=std::filesystem::exists(dir/GeneralExtremeFilename(i));
    if(!complete)return;
    std::ofstream manifest(dir/"manifest.csv");manifest<<std::setprecision(17);
    manifest<<"id,file,group,variant,type,rational,controls,degree_u,degree_v,collapsed_mask,length,width,aspect,sag,translation,coordinate_ulp,source_bound,qualified,edge_min,edge_max,edge_measure_change,orientation,seed\n";
    for(const auto &spec:specs)
    {
        auto c=ReadGeneralFixture(dir/GeneralExtremeFilename(spec.fixture.id));const auto &s=*c.input;
        double ulp=0,lo=INFINITY,hi=0,change=0;
        for(const auto &row:s.ControlPoints())for(auto point:row)for(double x:{point.X(),point.Y(),point.Z()})ulp=std::max(ulp,std::nextafter(std::abs(x),INFINITY)-std::abs(x));
        for(unsigned e=0;e<4;++e)
        {
            if(c.collapsedMask&(1u<<e))continue;
            auto local=c;local.truth.origin={0,0,0};double length=EdgeLength(local,e,2048);
            lo=std::min(lo,length);hi=std::max(hi,length);change=std::max(change,std::abs(length-EdgeLength(local,e,1024)));
        }
        manifest<<c.id<<','<<GeneralExtremeFilename(c.id)<<','<<spec.group<<','<<spec.variant<<','<<TypeName(c.truth.type)<<','<<c.rational<<','<<s.ControlPoints().size()*s.ControlPoints()[0].size()<<','<<s.DegreeU()<<','<<s.DegreeV()<<','<<c.collapsedMask<<','<<spec.length<<','<<spec.width<<','<<spec.length/spec.width<<','<<spec.sag<<','<<spec.translation<<','<<ulp<<','<<c.storedUpper<<','<<(c.storedUpper>=0)<<','<<lo<<','<<hi<<','<<change<<','<<spec.orientation<<','<<c.seed<<'\n';
    }
    Require(bool(manifest),"Manifest write failure");
}
}
