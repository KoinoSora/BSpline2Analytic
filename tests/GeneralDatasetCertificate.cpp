/**
 * @file GeneralDatasetCertificate.cpp
 * @brief 测试来源及四次多项式排除的独立证书实现。
 * 采用区间余项与秩证明，不由生产转换器的成功/失败决定输入标签。
 */
#include "GeneralDataset.h"
#include <functional>
#include <sstream>

namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
BPoly ConstantP(I x){return BPoly(0,0,x);}
// 以区间常数包围向外舍入的 Lagrange 余项并叠加 Taylor 多项式；此参考模型不调用 SDK 解析曲面求值。
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
BV TargetPolynomial(const GeneralFixture &c,const BezierCell &cell)
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
    x=Normalize(x-z*Dot(x,z));IV y=Cross(z,x),o=Exact(Get(a.origin));
    return {ConstantP(o.x)+local.x*x.x+local.y*y.x+local.z*z.x,
            ConstantP(o.y)+local.x*x.y+local.y*y.y+local.z*z.y,
            ConstantP(o.z)+local.x*x.z+local.y*y.z+local.z*z.z};
}
I Eval(const BPoly &p,I u,I v)
{
    std::vector<I> row(p.u+1);
    for(unsigned i=0;i<=p.u;++i)
    {
        std::vector<I> d(p.v+1);for(unsigned j=0;j<=p.v;++j)d[j]=p.at(i,j);
        for(unsigned n=p.v;n>0;--n)for(unsigned j=0;j<n;++j)d[j]=d[j]*(I(1)-v)+d[j+1]*v;
        row[i]=d[0];
    }
    for(unsigned n=p.u;n>0;--n)for(unsigned i=0;i<n;++i)row[i]=row[i]*(I(1)-u)+row[i+1]*u;
    return row[0];
}
IV EvalSurface(const std::vector<BezierCell> &cells,double u,double v)
{
    for(const auto &cell:cells)if(u>=cell.u0&&u<=cell.u1&&v>=cell.v0&&v<=cell.v1)
    {
        I s=(I(u)-I(cell.u0))/(I(cell.u1)-I(cell.u0)),t=(I(v)-I(cell.v0))/(I(cell.v1)-I(cell.v0));
        I w=Eval(cell.w,s,t);return {Eval(cell.p.x,s,t)/w,Eval(cell.p.y,s,t)/w,Eval(cell.p.z,s,t)/w};
    }
    throw std::runtime_error("Witness outside source domain");
}
using Exponent=std::array<unsigned,3>;
std::vector<Exponent> Exponents()
{
    std::vector<Exponent> out;
    for(unsigned degree=0;degree<=4;++degree)for(unsigned a=0;a<=degree;++a)for(unsigned b=0;b<=degree-a;++b)out.push_back({a,b,degree-a-b});
    return out;
}
std::vector<I> Monomials(IV p,const std::vector<Exponent> &exponents,int derivative=-1)
{
    std::array<std::array<I,5>,3> powers;
    std::array<I,3> x{p.x,p.y,p.z};
    for(unsigned j=0;j<3;++j){powers[j][0]=I(1);for(unsigned k=1;k<=4;++k)powers[j][k]=powers[j][k-1]*x[j];}
    std::vector<I> result;
    for(auto e:exponents)
    {
        I factor(1);
        if(derivative>=0){factor=I(e[derivative]);if(!e[derivative]){result.push_back(I(0));continue;}--e[derivative];}
        result.push_back(factor*powers[0][e[0]]*powers[1][e[1]]*powers[2][e[2]]);
    }
    return result;
}
}
FixtureCertificate CertifyAnalyticPair(const GeneralFixture &c,double tolerance,unsigned maxPieces)
{
    FixtureCertificate out;
    if(!IntervalEnvironmentSupported() || !(tolerance>0)){out.message="Interval environment/tolerance invalid";return out;}
    const auto &s=*c.input;
    if(s.MinParamU()!=0||s.MaxParamU()!=1||s.MinParamV()!=0||s.MaxParamV()!=1){out.message="Certificate expects complete [0,1]^2";return out;}
    // 直边测试数据的双线性参数图必须始终保持同一定向。
    if(c.bend==0)
    {
        int sign=0;
        for(double u:{0.,1.})for(double v:{0.,1.})
        {
            I ux=(I(c.corners[1].X())-I(c.corners[0].X()))*I(1-v)+(I(c.corners[3].X())-I(c.corners[2].X()))*I(v);
            I uy=(I(c.corners[1].Y())-I(c.corners[0].Y()))*I(1-v)+(I(c.corners[3].Y())-I(c.corners[2].Y()))*I(v);
            I vx=(I(c.corners[2].X())-I(c.corners[0].X()))*I(1-u)+(I(c.corners[3].X())-I(c.corners[1].X()))*I(u);
            I vy=(I(c.corners[2].Y())-I(c.corners[0].Y()))*I(1-u)+(I(c.corners[3].Y())-I(c.corners[1].Y()))*I(u);
            I det=ux*vy-uy*vx;int thisSign=det.lo>0?1:det.hi<0?-1:0;
            if(!thisSign || (sign && sign!=thisSign)){out.message="Folded/degenerate analytic UV chart";return out;}sign=thisSign;
        }
    }
    out.upper=0;
    std::function<bool(const BezierCell &,unsigned)> visit=[&](const BezierCell &cell,unsigned depth){
        if(++out.pieces>maxPieces)return false;
        BV reference=TargetPolynomial(c,cell),error=cell.p-reference*cell.w;
        IV range=Range(error);I norm=Sqrt(Sqr(range.x)+Sqr(range.y)+Sqr(range.z))/cell.w.range();
        if(norm.finite() && norm.hi<tolerance){out.upper=std::max(out.upper,norm.hi);return true;}
        if(depth>=12){std::ostringstream m;m<<"Pair bound unresolved at "<<cell.u0<<','<<cell.v0<<" bound="<<norm.hi;out.message=m.str();return false;}
        auto parts=Split(cell,depth%2==0);return visit(parts.first,depth+1)&&visit(parts.second,depth+1);
    };
    auto cells=ExtractCells(s);out.passed=!cells.empty();
    for(const auto &cell:cells)if(!visit(cell,0)){out.passed=false;break;}
    if(out.passed)out.message="Continuous same-parameter analytic/spline bound (including all four boundaries)";
    else{out.upper=INFINITY;if(out.message.empty())out.message="Pair certificate budget exhausted";}
    return out;
}
FixtureCertificate CertifyNoQuartic(const sggk::BSplineSurface &s,double radius)
{
    FixtureCertificate out;
    if(!IntervalEnvironmentSupported() || !(radius>0)){out.message="Invalid exclusion environment/radius";return out;}
    const auto exponents=Exponents();const unsigned n=35,grid=12,m=(grid+1)*(grid+1);
    auto cells=ExtractCells(s);std::vector<IV> points;
    for(unsigned i=0;i<=grid;++i)for(unsigned j=0;j<=grid;++j)points.push_back(EvalSurface(cells,double(i)/grid,double(j)/grid));
    out.witnesses=m;
    std::array<double,3> center{},scale{};
    for(unsigned d=0;d<3;++d)
    {
        double lo=INFINITY,hi=-INFINITY;
        for(auto p:points){I x=d==0?p.x:d==1?p.y:p.z;lo=std::min(lo,x.lo);hi=std::max(hi,x.hi);}
        center[d]=lo+(hi-lo)*.5;scale[d]=(hi-lo)*.5;
        if(!(scale[d]>0)){out.message="Witness coordinate rank deficiency";return out;}
    }
    std::vector<std::vector<I>> values(m);
    std::vector<IV> normalized;
    for(unsigned i=0;i<m;++i)
    {
        IV p{(points[i].x-I(center[0]))/I(scale[0]),(points[i].y-I(center[1]))/I(scale[1]),(points[i].z-I(center[2]))/I(scale[2])};
        normalized.push_back(p);values[i]=Monomials(p,exponents);
    }
    // 重正交 QR 只用于预条件处理，不是数值证书；最终秩下界由后续区间 Gram 矩阵独立证明。
    std::vector<std::vector<long double>> q(m,std::vector<long double>(n)),r(n,std::vector<long double>(n));
    for(unsigned j=0;j<n;++j)
    {
        for(unsigned i=0;i<m;++i)q[i][j]=((long double)values[i][j].lo+values[i][j].hi)/2;
        for(unsigned pass=0;pass<2;++pass)for(unsigned k=0;k<j;++k)
        {
            long double dot=0;for(unsigned i=0;i<m;++i)dot+=q[i][k]*q[i][j];r[k][j]+=dot;
            for(unsigned i=0;i<m;++i)q[i][j]-=dot*q[i][k];
        }
        long double norm=0;for(unsigned i=0;i<m;++i)norm+=q[i][j]*q[i][j];r[j][j]=std::sqrt(norm);
        if(!(r[j][j]>1e-17L)){out.message="Degree-four witness matrix deficient";return out;}
        for(unsigned i=0;i<m;++i)q[i][j]/=r[j][j];
    }
    std::vector<std::vector<double>> basis(n,std::vector<double>(n));
    for(unsigned j=0;j<n;++j)for(int i=int(j);i>=0;--i)
    {
        long double v=unsigned(i)==j?1:0;for(unsigned k=unsigned(i)+1;k<=j;++k)v-=r[i][k]*basis[k][j];
        basis[i][j]=double(v/r[i][i]);
    }
    auto transform=[&](const std::vector<I> &v){std::vector<I> a(n);for(unsigned j=0;j<n;++j)for(unsigned k=0;k<=j;++k)a[j]=a[j]+v[k]*I(basis[k][j]);return a;};
    std::vector<std::vector<I>> a(m);for(unsigned i=0;i<m;++i)a[i]=transform(values[i]);
    std::vector<std::vector<I>> gram(n,std::vector<I>(n));
    for(unsigned j=0;j<n;++j)for(unsigned k=j;k<n;++k)
    {
        I sum(0);for(unsigned i=0;i<m;++i)sum=sum+(j==k?Sqr(a[i][j]):a[i][j]*a[i][k]);
        gram[j][k]=gram[k][j]=sum;
    }
    double eigenLower=INFINITY;
    for(unsigned j=0;j<n;++j){I lower=gram[j][j];for(unsigned k=0;k<n;++k)if(k!=j)lower=lower-I(gram[j][k].absMax());eigenLower=std::min(eigenLower,lower.lo);}
    if(!(eigenLower>0)){out.message="Interval Gram rank bound unresolved";return out;}
    I squaredGradient(0);
    for(unsigned i=0;i<m;++i)
    {
        IV p=normalized[i];I *coordinates[]={&p.x,&p.y,&p.z};
        for(unsigned d=0;d<3;++d)*coordinates[d]=*coordinates[d]+I(-radius,radius)/I(scale[d]);
        for(unsigned d=0;d<3;++d)
        {
            auto derivative=transform(Monomials(p,exponents,int(d)));
            for(auto component:derivative)squaredGradient=squaredGradient+Sqr(component/I(scale[d]));
        }
    }
    I margin=Sqrt(I(eigenLower))-I(radius)*Sqrt(squaredGradient);
    out.rankMargin=margin.lo;out.passed=margin.finite()&&margin.lo>0;
    std::ostringstream message;message<<"169 source witnesses; sigma lower="<<Sqrt(I(eigenLower)).lo<<" perturbation upper="<<(I(radius)*Sqrt(squaredGradient)).hi<<" rank margin="<<margin.lo;
    out.message=message.str();return out;
}
}
