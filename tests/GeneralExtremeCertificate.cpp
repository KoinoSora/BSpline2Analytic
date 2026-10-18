/**
 * @file GeneralExtremeCertificate.cpp
 * @brief 历史一般极端输入的独立来源认证工具。
 * 包含三角参数图与区间包络，源码保留不代表允许重跑归档数据。
 */
#include "GeneralExtremeDataset.h"
#include <functional>
#include <sstream>

namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
BPoly ConstantP(I x){return BPoly(0,0,x);}
// Taylor 多项式加向外舍入 Lagrange 余项，余项以区间常数表示；参考模型不调用 SDK 解析求值。
BPoly TrigP(const BPoly &angle,bool cosine)
{
    I range=angle.range();double center=range.lo+(range.hi-range.lo)*.5;
    BPoly d=angle-ConstantP(I(center)),power=ConstantP(I(1)),sum;
    I sin=Sin(I(center)),cos=Cos(I(center)),factorial(1);
    for(unsigned n=0;n<=9;++n)
    {
        I coefficient;
        switch((n+(cosine?1:0))%4){case 0:coefficient=sin;break;case 1:coefficient=cos;break;case 2:coefficient=-sin;break;default:coefficient=-cos;}
        sum=sum+power*(coefficient/factorial);
        power=power*d;factorial=factorial*I(n+1);
    }
    I remainder(1),radius(d.range().absMax());
    for(unsigned j=1;j<=10;++j)remainder=remainder*radius/I(j);
    return sum+ConstantP(I(-remainder.hi,remainder.hi));
}
BV TargetPolynomial(const GeneralFixture &c,const BezierCell &cell,V anchor)
{
    BPoly s(1,0),t(0,1);s.at(0,0)=I(cell.u0);s.at(1,0)=I(cell.u1);
    t.at(0,0)=I(cell.v0);t.at(0,1)=I(cell.v1);
    auto one=ConstantP(I(1));
    std::array<BPoly,4> weights={(one-s)*(one-t),s*(one-t),(one-s)*t,s*t};
    BPoly u,v;for(unsigned i=0;i<4;++i){u=u+weights[i]*I(c.corners[i].X());v=v+weights[i]*I(c.corners[i].Y());}
    if(c.bend!=0)v=v+s*(one-s)*I(c.bend);
    BV local;
    const auto &a=c.truth;
    switch(a.type)
    {
    case sggk::SurfaceType::Plane:local={u,v,ConstantP(I(0))};break;
    case sggk::SurfaceType::Cylinder:local={TrigP(u,true)*I(a.radius),TrigP(u,false)*I(a.radius),v};break;
    case sggk::SurfaceType::Cone:{I slope=Sin(I(a.halfAngle))/Cos(I(a.halfAngle));local={v*TrigP(u,true)*slope,v*TrigP(u,false)*slope,v};break;}
    case sggk::SurfaceType::Sphere:local={TrigP(v,false)*TrigP(u,true)*I(a.radius),TrigP(v,false)*TrigP(u,false)*I(a.radius),TrigP(v,true)*I(-a.radius)};break;
    case sggk::SurfaceType::Torus:{BPoly r=ConstantP(I(a.majorRadius))+TrigP(v,true)*I(a.radius);local={r*TrigP(u,true),r*TrigP(u,false),TrigP(v,false)*I(a.radius)};break;}
    default:throw std::runtime_error("Unsupported certificate mother");
    }
    IV z=Normalize(Exact(Get(a.axis))),x=Exact(Get(a.xDirection));
    x=Normalize(x-z*Dot(x,z));IV y=Cross(z,x),o=Exact(Get(a.origin))-Exact(anchor);
    return {ConstantP(o.x)+local.x*x.x+local.y*y.x+local.z*z.x,
            ConstantP(o.y)+local.x*x.y+local.y*y.y+local.z*z.y,
            ConstantP(o.z)+local.x*x.z+local.y*y.z+local.z*z.z};
}
bool ValidChart(const GeneralFixture &c)
{
    if(c.bend!=0)return false;
    const unsigned ends[4][2]={{0,2},{1,3},{0,1},{2,3}};
    std::vector<sggk::Point2D> unique;
    unsigned collapsed=0;
    for(auto edge:ends)
        collapsed+=c.corners[edge[0]].X()==c.corners[edge[1]].X()&&c.corners[edge[0]].Y()==c.corners[edge[1]].Y();
    for(auto p:c.corners)
    {
        bool exists=false;for(auto q:unique)exists|=p.X()==q.X()&&p.Y()==q.Y();
        if(!exists)unique.push_back(p);
    }
    if(collapsed==1&&unique.size()==3)
    {
        // 恰有一对相邻参数角点相等；双线性行列式为非零三角形行列式乘 s、1-s、t 或 1-t。
        I ax=I(unique[1].X())-I(unique[0].X()),ay=I(unique[1].Y())-I(unique[0].Y());
        I bx=I(unique[2].X())-I(unique[0].X()),by=I(unique[2].Y())-I(unique[0].Y());
        I det=ax*by-ay*bx;return det.lo>0||det.hi<0;
    }
    if(unique.size()!=4)return false;
    int sign=0;
    for(double u:{0.,1.})for(double v:{0.,1.})
    {
        I ux=(I(c.corners[1].X())-I(c.corners[0].X()))*I(1-v)+(I(c.corners[3].X())-I(c.corners[2].X()))*I(v);
        I uy=(I(c.corners[1].Y())-I(c.corners[0].Y()))*I(1-v)+(I(c.corners[3].Y())-I(c.corners[2].Y()))*I(v);
        I vx=(I(c.corners[2].X())-I(c.corners[0].X()))*I(1-u)+(I(c.corners[3].X())-I(c.corners[1].X()))*I(u);
        I vy=(I(c.corners[2].Y())-I(c.corners[0].Y()))*I(1-u)+(I(c.corners[3].Y())-I(c.corners[1].Y()))*I(u);
        I det=ux*vy-uy*vx;int current=det.lo>0?1:det.hi<0?-1:0;
        if(!current||(sign&&sign!=current))return false;sign=current;
    }
    return true;
}
}
FixtureCertificate CertifyGeneralExtremePair(const GeneralFixture &c,double tolerance,unsigned maxPieces)
{
    FixtureCertificate out;
    if(!IntervalEnvironmentSupported()||!c.input||!(tolerance>0)||!ValidChart(c))
    {out.message="Invalid environment or folded/degenerate UV chart";return out;}
    const auto &s=*c.input;
    if(s.MinParamU()!=0||s.MaxParamU()!=1||s.MinParamV()!=0||s.MaxParamV()!=1)
    {out.message="Certificate expects [0,1]^2";return out;}
    auto cp=s.ControlPoints();V anchor=Get(cp[0][0]);double shiftError=0;
    for(auto &row:cp)for(auto &point:row)
    {
        V value=Get(point)-anchor;
        IV e=Exact(Get(point))-Exact(anchor)-Exact(value);
        shiftError=std::max(shiftError,Sqrt(Sqr(e.x)+Sqr(e.y)+Sqr(e.z)).hi);
        point=Point(value);
    }
    sggk::BSplineSurfacePtr shifted;
    if(s.IsRational())shifted=std::make_shared<sggk::BSplineSurface>(s.DegreeU(),s.DegreeV(),cp,s.Weights(),s.KnotsU(),s.KnotsV(),s.MultsU(),s.MultsV());
    else shifted=std::make_shared<sggk::BSplineSurface>(s.DegreeU(),s.DegreeV(),cp,s.KnotsU(),s.KnotsV(),s.MultsU(),s.MultsV());
    double angular[2]{};
    for(unsigned axis=0;axis<2;++axis)
        for(unsigned k=0;k<2;++k)
        {
            unsigned a=axis?k:k*2,b=a+(axis?2:1);
            angular[axis]=std::max(angular[axis],std::abs(c.corners[a].X()-c.corners[b].X()));
            if(c.truth.type==sggk::SurfaceType::Sphere||c.truth.type==sggk::SurfaceType::Torus)
                angular[axis]=std::max(angular[axis],std::abs(c.corners[a].Y()-c.corners[b].Y()));
        }
    out.upper=0;
    std::function<bool(const BezierCell &,unsigned)> visit=[&](const BezierCell &cell,unsigned depth) {
        if(++out.pieces>maxPieces){out.message="Pair certificate piece budget exhausted";return false;}
        BV error=cell.p-TargetPolynomial(c,cell,anchor)*cell.w;
        IV range=Range(error);I norm=Sqrt(Sqr(range.x)+Sqr(range.y)+Sqr(range.z))/cell.w.range()+I(shiftError);
        if(norm.finite()&&norm.hi<tolerance){out.upper=std::max(out.upper,norm.hi);return true;}
        if(depth>=16){std::ostringstream m;m<<"Pair bound unresolved: "<<norm.hi;out.message=m.str();return false;}
        bool alongU=(cell.u1-cell.u0)*angular[0]>=(cell.v1-cell.v0)*angular[1];
        if(angular[0]==0&&angular[1]==0)alongU=depth%2==0;
        auto parts=Split(cell,alongU);return visit(parts.first,depth+1)&&visit(parts.second,depth+1);
    };
    auto cells=ExtractCells(*shifted);out.passed=!cells.empty();
    for(const auto &cell:cells)if(!visit(cell,0)){out.passed=false;break;}
    if(out.passed)out.message="Continuous recentered spline/analytic correspondence including four boundaries";
    else out.upper=INFINITY;
    return out;
}
}
