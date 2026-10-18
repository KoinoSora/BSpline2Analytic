/**
 * @file GeneralBoundaryError.cpp
 * @brief 有理边界位置误差与正比例切向残差的区间代数。
 * 源/目标统一到相同局部参数后构造误差；微小段豁免仍需完整覆盖、双侧直径和严格位置证书。
 */
#include "GeneralRationalBoundaryInternal.h"

namespace bspline2analytic::detail::rational_boundary
{
namespace
{
I Choose(unsigned n,unsigned k)
{
    I result(1);k=std::min(k,n-k);
    for(unsigned i=1;i<=k;++i)result=result*I(n-k+i)/I(i);
    return result;
}
bool Homogeneous1D(const BezierCell &c)
{
    return !c.p.x.v&&!c.p.y.v&&!c.p.z.v&&!c.w.v&&
        c.p.x.u==c.w.u&&c.p.y.u==c.w.u&&c.p.z.u==c.w.u;
}
I Length(IV a){return Sqrt(Sqr(a.x)+Sqr(a.y)+Sqr(a.z));}
I DerivativeRange(const BPoly &p)
{
    if(p.v)return I::Whole();
    if(!p.u)return I(0);
    I r=I(p.u)*(p.c[1]-p.c[0]);
    for(unsigned i=1;i<p.u;++i)
    {I x=I(p.u)*(p.c[i+1]-p.c[i]);r.lo=std::min(r.lo,x.lo);r.hi=std::max(r.hi,x.hi);}
    return r;
}
double Estimate(const BPoly &p,double t)
{
    if(p.v||p.c.size()!=p.u+1)return NAN;
    auto midpoint=[](I x){return x.lo*.5+x.hi*.5;};
    if(t==0)return midpoint(p.c.front());
    if(t==1)return midpoint(p.c.back());
    std::vector<double> values;values.reserve(p.c.size());
    for(I x:p.c)values.push_back(midpoint(x));
    for(size_t n=values.size();n>1;--n)for(size_t i=0;i+1<n;++i)
        values[i]=(1-t)*values[i]+t*values[i+1];
    return values.front();
}
V Estimate(const BV &p,double t)
{return {Estimate(p.x,t),Estimate(p.y,t),Estimate(p.z,t)};}
double Magnitude(V p){return std::max({std::abs(p.x),std::abs(p.y),std::abs(p.z)});}
I LengthSquared(IV p){return Sqr(p.x)+Sqr(p.y)+Sqr(p.z);}
}
// 源 S=Ps/Ws、目标 T=Pt/Wt 时，位置误差 T-S=(Pt*Ws-Ps*Wt)/(Ws*Wt)。
// 只有源/目标已限制到同一局部参数后才能相减。融合卷积一次累加三维分子和分母，
// 复用二项式权重；常权重或不匹配的齐次次数交回通用乘积路径。
CurveError CurveErrors::Build(const BezierCell &s,const BezierCell &t,paired::Products &products,bool &fused)
{
    fused=false;
    // 保留常权重路径；仅为合并乘积而将多项式权重升阶，会增加运算量并无谓放宽分母区间。
    if(!s.w.u||!t.w.u||!Homogeneous1D(s)||!Homogeneous1D(t))
        return {products.Subtract(products.MulShared(t.p,s.w),products.MulShared(s.p,t.w)),
                products.Mul(t.w,s.w)};
    unsigned p=s.w.u,q=t.w.u;
    auto key=std::make_pair(p,q);auto found=factors.find(key);
    if(found==factors.end())
    {
        std::vector<I> weights((p+1)*(q+1));
        for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
            weights[i*(q+1)+j]=Choose(p,i)*Choose(q,j)/Choose(p+q,i+j);
        found=factors.emplace(key,std::move(weights)).first;
    }
    CurveError out{{BPoly(p+q,0),BPoly(p+q,0),BPoly(p+q,0)},BPoly(p+q,0)};
    const auto &weights=found->second;
    for(unsigned i=0;i<=p;++i)for(unsigned j=0;j<=q;++j)
    {
        I f=weights[i*(q+1)+j],sw=s.w.c[i],tw=t.w.c[j];unsigned k=i+j;
        out.numerator.x.c[k]=out.numerator.x.c[k]+(t.p.x.c[j]*sw-s.p.x.c[i]*tw)*f;
        out.numerator.y.c[k]=out.numerator.y.c[k]+(t.p.y.c[j]*sw-s.p.y.c[i]*tw)*f;
        out.numerator.z.c[k]=out.numerator.z.c[k]+(t.p.z.c[j]*sw-s.p.z.c[i]*tw)*f;
        out.denominator.c[k]=out.denominator.c[k]+(sw*tw)*f;
    }
    fused=true;return out;
}
// 微小段合并不是删除短边：要求完整覆盖 [lo,hi]，源与目标各自空间直径均 <= 阈值，
// 且两包围盒任意点对的距离严格小于位置容差。只有切向检查被豁免，并返回具体区间。
// 目标跨接缝时逐段取并集包围盒，任何覆盖缺口都导致无法豁免。
std::optional<GeneralBoundaryTangentExemption> MicroSegmentBound(const BezierCell &source,
    const std::vector<BezierCell> &target,double lo,double hi,double diameterTolerance,double distanceTolerance)
{
    if(!(diameterTolerance>0&&lo<=hi))return {};
    auto box=[](const BezierCell &c)->std::optional<IV>
    {
        I w=c.w.range();if(!w.finite()||!(w.lo>0))return {};
        IV p=Range(c.p);p={p.x/w,p.y/w,p.z/w};
        if(!p.x.finite()||!p.y.finite()||!p.z.finite())return {};
        return p;
    };
    auto diameter=[](IV b){return Length({I(b.x.hi)-I(b.x.lo),I(b.y.hi)-I(b.y.lo),I(b.z.hi)-I(b.z.lo)}).hi;};
    auto s=box(source);if(!s)return {};
    double sd=diameter(*s);if(!(sd<=diameterTolerance))return {};
    std::optional<IV> t;double covered=lo;
    auto hull=[](I a,I b){return I(std::min(a.lo,b.lo),std::max(a.hi,b.hi));};
    for(const auto &part:target)
    {
        double first=std::max(lo,part.u0),last=std::min(hi,part.u1);
        if(first>last)continue;
        if(first>covered)return {}; // 微小分段也不能掩盖参数覆盖缺口。
        auto b=box(paired::Crop(part,first,last,0,1));if(!b)return {};
        if(t)t=IV{hull(t->x,b->x),hull(t->y,b->y),hull(t->z,b->z)};else t=b;
        covered=std::max(covered,last);
    }
    if(!t||covered<hi)return {};
    double td=diameter(*t),distance=Length(*s-*t).hi;
    if(!(td<=diameterTolerance&&distance<distanceTolerance))return {};
    GeneralBoundaryTangentExemption result;
    result.sourceFirst=source.u0;result.sourceLast=source.u1;result.targetFirst=lo;result.targetLast=hi;
    result.sourceDiameterUpper=sd;result.targetDiameterUpper=td;result.distanceUpper=distance;
    return result;
}
BV TangentNumerator(const BezierCell &c,paired::Products &products)
{
    if(!c.w.u&&!c.w.v)return Derivative(c.p,true)*c.w.c[0];
    return products.Subtract(products.MulShared(Derivative(c.p,true),c.w),
                             products.MulShared(c.p,Derivative(c.w,true)));
}
double ResidualTangentRatio(const CurveError &e,const BezierCell &s,const BV &st,double distanceUpper)
{
    I h=e.denominator.range(),w=s.w.range();
    if(!h.finite()||h.lo<=0||!w.finite()||w.lo<=0||
       !std::isfinite(distanceUpper)||distanceUpper<0)return INFINITY;
    // 由 E'=G'/H-E*H'/H 求误差导数，并复用已认证的位置上界控制 ||E||。
    I g1=Length({DerivativeRange(e.numerator.x),DerivativeRange(e.numerator.y),DerivativeRange(e.numerator.z)});
    I h1=DerivativeRange(e.denominator);
    I upper=(I(g1.hi)+I(distanceUpper)*I(h1.absMax()))/I(h.lo);
    // 系数区间盒提供源切向速度的下界；区间盒含零只能说明无法判定，不能推出曲线不正则。
    I speed=Length(Range(st))/Sqr(I(w.hi));
    if(!upper.finite()||!speed.finite()||!(upper.hi<speed.lo))return INFINITY;
    // 利用 T'=S'+E'，在源速度下界 m 大于误差导数上界 M 时，有 dot(S',T')>0 且 tan(theta)<=M/(m-M)。
    return (I(upper.hi)/(I(speed.lo)-I(upper.hi))).hi;
}
std::optional<BPoly> ChooseTangentScale(const BV &source,const BV &target)
{
    double ratio[3];
    for(unsigned i=0;i<3;++i)
    {
        V a=Estimate(source,.5*i),b=Estimate(target,.5*i);
        double ma=Magnitude(a),mb=Magnitude(b);
        if(!Finite(a)||!Finite(b)||!(ma>0&&mb>0))return {};
        a={a.x/ma,a.y/ma,a.z/ma};b={b.x/mb,b.y/mb,b.z/mb};
        ratio[i]=(ma/mb)*(Dot(a,b)/Dot(b,b));
        if(!(ratio[i]>0)||!std::isfinite(ratio[i]))return {};
    }
    // 存储的双精度值精确定义辅助比例 lambda；启发式选值的舍入不影响证明，因为后续对完整残差另行包络。
    double tolerance=64*std::numeric_limits<double>::epsilon()*std::max({ratio[0],ratio[1],ratio[2]});
    BPoly scale(2,0);scale.c={I(ratio[0]),I(2*ratio[1]-.5*ratio[0]-.5*ratio[2]),I(ratio[2])};
    if(std::max(std::abs(ratio[0]-ratio[1]),std::abs(ratio[2]-ratio[1]))<=tolerance||
       !scale.range().finite()||!(scale.range().lo>0))
        scale=BPoly(0,0,I(ratio[1]));
    else if(std::abs(ratio[1]-.5*ratio[0]-.5*ratio[2])<=tolerance)
    {scale=BPoly(1,0);scale.c={I(ratio[0]),I(ratio[2])};}
    return scale;
}
BV BuildScaledTangentResidual(const BV &source,const BV &target,const BPoly &scale,paired::Products &products)
{return products.Subtract(source,products.MulShared(target,scale));}
V TangentResidualReference(const BV &source)
{
    V e=Estimate(source,.5);double size=Magnitude(e);
    if(!Finite(e)||!(size>0))return {};
    // 按最大分量相除，避免 e 很小时先计算 1/size 引发溢出。
    return {e.x/size,e.y/size,e.z/size};
}
double ScaledTangentSineSquaredUpper(const BV &source,const BV &residual,const BPoly &scale)
{
    I lambda=scale.range();if(!lambda.finite()||!(lambda.lo>0))return INFINITY;
    V e=TangentResidualReference(source);I e2=LengthSquared(Exact(e));
    if(!e2.finite()||!(e2.lo>0))return INFINITY;
    I projection=Dot(source,Exact(e)).range();
    if(!projection.finite()||!(projection.lo>0))return INFINITY;
    I length2=Sqr(I(projection.lo))/e2;
    if(!length2.finite()||!(length2.lo>0))return INFINITY;
    double error2=0;
    if(residual.x.u==residual.y.u&&residual.x.u==residual.z.u&&
       residual.x.v==residual.y.v&&residual.x.v==residual.z.v)
    {
        // 利用范数的凸性直接对向量控制系数求界，不构造范数平方样条乘积。
        for(size_t i=0;i<residual.x.c.size();++i)
        {
            I norm2=LengthSquared({residual.x.c[i],residual.y.c[i],residual.z.c[i]});
            if(!norm2.finite())return INFINITY;error2=std::max(error2,norm2.hi);
        }
    }
    else
    {
        I norm2=LengthSquared(Range(residual));if(!norm2.finite())return INFINITY;error2=norm2.hi;
    }
    I q=I(error2)/I(length2.lo);
    // 残差上界 E<||A|| 且 lambda>0，可证明 B 非零且 A.B>0，排除反向平行切向场。
    return q.finite()&&q.hi<1?std::max(0.,q.hi):INFINITY;
}
}
