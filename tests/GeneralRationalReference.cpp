#include "GeneralRationalDataset.h"
#include <algorithm>

namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
BPoly C(I x){return BPoly(0,0,x);}
double Mid(I x){if(!x.finite())throw std::runtime_error("Nonfinite reference coefficient");return x.lo+(x.hi-x.lo)*.5;}
BPoly Map(const GeneralFixture &f,unsigned d,double u0,double u1,double v0,double v1)
{
    double a[4];for(unsigned k=0;k<4;++k)a[k]=d?f.corners[k].Y():f.corners[k].X();
    unsigned p=(a[0]==a[1]&&a[2]==a[3])?0:1,q=(a[0]==a[2]&&a[1]==a[3])?0:1;
    if(u0==u1)p=0;if(v0==v1)q=0;
    BPoly out(p,q);
    for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
    {
        I u(i?u1:u0),v(j?v1:v0);
        out.at(i,j)=I(a[0])*(I(1)-u)*(I(1)-v)+I(a[1])*u*(I(1)-v)+I(a[2])*(I(1)-u)*v+I(a[3])*u*v;
    }
    if(d&&f.bend)
    {
        BPoly u(u0==u1?0:1,0);u.c=u0==u1?std::vector<I>{I(u0)}:std::vector<I>{I(u0),I(u1)};
        out=out+u*(C(I(1))-u)*I(f.bend);
    }
    return out;
}
std::vector<BPoly> Basis(unsigned degree,const BPoly &x)
{
    std::vector<BPoly> b{C(I(1))};auto one=C(I(1))-x;
    for(unsigned k=1;k<=degree;++k)
    {
        std::vector<BPoly> next(k+1);
        for(unsigned j=0;j<=k;++j){if(j<k)next[j]=b[j]*one;if(j)next[j]=next[j]+b[j-1]*x;}
        b=std::move(next);
    }
    return b;
}
BPoly Compose(const BPoly &p,const BPoly &u,const BPoly &v)
{
    auto bu=Basis(p.u,u),bv=Basis(p.v,v);BPoly out;
    for(unsigned i=0;i<=p.u;++i)for(unsigned j=0;j<=p.v;++j)out=out+(bu[i]*bv[j])*p.at(i,j);
    return out;
}
struct Axis{BPoly c,s,w,t;};
Axis AxisReference(double a,double b,bool angular,bool u)
{
    if(!angular){BPoly t(u?1:0,u?0:1);t.c={I(a),I(b)};return {C(I(0)),C(I(0)),C(I(1)),t};}
    I mid=(I(a)+I(b))/I(2),half=(I(b)-I(a))/I(2);
    BPoly c(u?2:0,u?0:2),s=c,w=c;
    // Homogeneous quadratic circular arc; no inverse mapping or fitted controls.
    c.c={Cos(I(a)),Cos(mid),Cos(I(b))};s.c={Sin(I(a)),Sin(mid),Sin(I(b))};w.c={I(1),Cos(half),I(1)};
    if(w.range().lo<=0)throw std::runtime_error("Reference arc weight is not positive");
    return {c,s,w,C(I(0))};
}
double Eval(const BPoly &p,double u,double v)
{
    std::vector<double> rows(p.u+1);
    for(unsigned i=0;i<=p.u;++i)
    {
        std::vector<double> a(p.v+1);for(unsigned j=0;j<=p.v;++j)a[j]=Mid(p.at(i,j));
        for(unsigned n=p.v;n;--n)for(unsigned j=0;j<n;++j)a[j]=(1-v)*a[j]+v*a[j+1];rows[i]=a[0];
    }
    for(unsigned n=p.u;n;--n)for(unsigned i=0;i<n;++i)rows[i]=(1-u)*rows[i]+u*rows[i+1];return rows[0];
}
}
RationalReference::RationalReference(const GeneralFixture &f):spec(f)
{
    if(!f.rationalDomain)throw std::runtime_error("Retired analytic-UV dataset: rational-domain V2 required");
    bool angular[2]={f.type!=0,f.type==2||f.type==4};unsigned count[2];
    for(unsigned d=0;d<2;++d)
    {
        double span=f.chartDomain[2*d+1]-f.chartDomain[2*d];
        if(!(span>0&&std::isfinite(span))||(angular[d]&&span>2*Pi+1e-10))throw std::runtime_error("Invalid reference chart span");
        count[d]=angular[d]?std::max(1u,unsigned(std::ceil(span/(Pi/2)))):1;
    }
    const auto &a=f.truth;
    for(unsigned i=0;i<count[0];++i)for(unsigned j=0;j<count[1];++j)
    {
        auto endpoint=[&](unsigned d,unsigned k){return k==count[d]?f.chartDomain[2*d+1]:f.chartDomain[2*d]+(f.chartDomain[2*d+1]-f.chartDomain[2*d])*k/count[d];};
        auto u=AxisReference(endpoint(0,i),endpoint(0,i+1),angular[0],true);
        auto v=AxisReference(endpoint(1,j),endpoint(1,j+1),angular[1],false);
        BV p;BPoly w=u.w*v.w;
        switch(f.type)
        {
        case 0:p={u.t,v.t,C(I(0))};break;
        case 1:p={u.c*I(a.radius),u.s*I(a.radius),u.w*v.t};break;
        case 2:p={(u.c*v.s)*I(a.radius),(u.s*v.s)*I(a.radius),(u.w*v.c)*I(-a.radius)};break;
        case 3:{I slope=Sin(I(a.halfAngle))/Cos(I(a.halfAngle));p={(u.c*v.t)*slope,(u.s*v.t)*slope,u.w*v.t};break;}
        case 4:{auto radial=v.w*I(a.majorRadius)+v.c*I(a.radius);p={u.c*radial,u.s*radial,(u.w*v.s)*I(a.radius)};break;}
        default:throw std::runtime_error("Invalid reference mother type");
        }
        chart.push_back({p,double(i)/count[0],double(i+1)/count[0],double(j)/count[1],double(j+1)/count[1],w});
    }
}
V RationalReference::Value(double s,double t) const
{
    double a[2]{};double weights[]={(1-s)*(1-t),s*(1-t),(1-s)*t,s*t};
    for(unsigned i=0;i<4;++i){a[0]+=weights[i]*spec.corners[i].X();a[1]+=weights[i]*spec.corners[i].Y();}a[1]+=spec.bend*s*(1-s);
    const BezierCell *found=nullptr;
    for(const auto &c:chart)if(a[0]>=c.u0-1e-12&&a[0]<=c.u1+1e-12&&a[1]>=c.v0-1e-12&&a[1]<=c.v1+1e-12){found=&c;break;}
    if(!found)throw std::runtime_error("Reference point outside chart");const auto &c=*found;
    double u=(a[0]-c.u0)/(c.u1-c.u0),v=(a[1]-c.v0)/(c.v1-c.v0),w=Eval(c.w,u,v);
    V p{Eval(c.p.x,u,v)/w,Eval(c.p.y,u,v)/w,Eval(c.p.z,u,v)/w};
    V z=Unit(Get(spec.truth.axis)),x=Get(spec.truth.xDirection);x=Unit(x-z*Dot(x,z));
    return Get(spec.truth.origin)+x*p.x+Cross(z,x)*p.y+z*p.z;
}
std::vector<BezierCell> RationalReference::On(double u0,double u1,double v0,double v1,V anchor) const
{
    auto u=Map(spec,0,u0,u1,v0,v1),v=Map(spec,1,u0,u1,v0,v1);auto ur=u.range(),vr=v.range();std::vector<BezierCell> out;
    IV z=Normalize(Exact(Get(spec.truth.axis))),x=Exact(Get(spec.truth.xDirection));x=Normalize(x-z*Dot(x,z));IV y=Cross(z,x),o=Exact(Get(spec.truth.origin))-Exact(anchor);
    for(const auto &c:chart)
    {
        if(ur.hi<c.u0||ur.lo>c.u1||vr.hi<c.v0||vr.lo>c.v1)continue;
        auto U=(u-C(I(c.u0)))* (I(1)/(I(c.u1)-I(c.u0))),Vv=(v-C(I(c.v0)))*(I(1)/(I(c.v1)-I(c.v0)));
        auto px=Compose(c.p.x,U,Vv),py=Compose(c.p.y,U,Vv),pz=Compose(c.p.z,U,Vv),w=Compose(c.w,U,Vv);
        BV p={px*x.x+py*y.x+pz*z.x+w*o.x,px*x.y+py*y.y+pz*z.y+w*o.y,px*x.z+py*y.z+pz*z.z+w*o.z};
        out.push_back({p,u0,u1,v0,v1,w});
    }
    if(out.empty())throw std::runtime_error("Reference block has no chart coverage");return out;
}
std::vector<RationalReferencePiece> RationalReference::Cover(const BezierCell &source,V anchor) const
{
    auto targets=On(source.u0,source.u1,source.v0,source.v1,anchor);
    auto ur=Map(spec,0,source.u0,source.u1,source.v0,source.v1).range();
    auto vr=Map(spec,1,source.u0,source.u1,source.v0,source.v1).range();
    std::vector<RationalReferencePiece> out;size_t index=0;
    for(const auto &cell:chart)
    {
        if(ur.hi<cell.u0||ur.lo>cell.u1||vr.hi<cell.v0||vr.lo>cell.v1)continue;
        auto target=std::move(targets[index++]);
        std::array<double,4> domain{source.u0,source.u1,source.v0,source.v1};bool empty=false;
        for(unsigned d=0;d<2;++d)
        {
            // Tensor-aligned seams can be pulled back exactly with an outward
            // rounded affine division. Keep a small OVERLAP, never a seam gap.
            if(d&&spec.bend)continue;
            double a[4];for(unsigned k=0;k<4;++k)a[k]=d?spec.corners[k].Y():spec.corners[k].X();
            bool uOnly=a[0]==a[2]&&a[1]==a[3],vOnly=a[0]==a[1]&&a[2]==a[3];
            if(!uOnly&&!vOnly)continue;unsigned axis=uOnly?0:1;
            I slope=I(a[axis?2:1])-I(a[0]);if(slope.lo<=0&&slope.hi>=0)continue;
            I lo=(I(d?cell.v0:cell.u0)-I(a[0]))/slope;
            I hi=(I(d?cell.v1:cell.u1)-I(a[0]))/slope;
            if(slope.hi<0)std::swap(lo,hi);
            domain[2*axis]=std::max(domain[2*axis],lo.lo);
            domain[2*axis+1]=std::min(domain[2*axis+1],hi.hi);
            if(domain[2*axis]>domain[2*axis+1]){empty=true;break;}
        }
        if(!empty)out.push_back({std::move(target),domain});
    }
    if(out.empty())throw std::runtime_error("Reference seam coverage empty");return out;
}
std::array<std::vector<double>,2> RationalReference::TensorBreaks() const
{
    std::array<std::vector<double>,2> breaks{{{0,1},{0,1}}};
    for(unsigned d=0;d<2;++d)
    {
        double a[4];for(unsigned i=0;i<4;++i)a[i]=d?spec.corners[i].Y():spec.corners[i].X();
        bool uOnly=a[0]==a[2]&&a[1]==a[3],vOnly=a[0]==a[1]&&a[2]==a[3];
        for(const auto &c:chart)
        {
            double knot=d?c.v1:c.u1;if(knot>=1)continue;
            auto mm=std::minmax_element(a,a+4);if(knot<=*mm.first||knot>=*mm.second)continue;
            if(!uOnly&&!vOnly)throw std::runtime_error("Non-tensor chart seam requires constrained fitting");
            unsigned axis=uOnly?0:1;double f=(knot-a[0])/(a[axis?2:1]-a[0]);if(f>0&&f<1)breaks[axis].push_back(f);
        }
    }
    for(auto &b:breaks){std::sort(b.begin(),b.end());b.erase(std::unique(b.begin(),b.end()),b.end());}return breaks;
}
}
