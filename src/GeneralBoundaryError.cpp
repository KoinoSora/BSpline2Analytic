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
}
CurveError CurveErrors::Build(const BezierCell &s,const BezierCell &t,paired::Products &products,bool &fused)
{
    fused=false;
    // Retain constant-weight paths: elevating a polynomial weight just to fuse
    // products would add work and widen its denominator unnecessarily.
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
    // E'=G'/H-E H'/H. Use the already certified position bound for ||E||.
    I g1=Length({DerivativeRange(e.numerator.x),DerivativeRange(e.numerator.y),DerivativeRange(e.numerator.z)});
    I h1=DerivativeRange(e.denominator);
    I upper=(I(g1.hi)+I(distanceUpper)*I(h1.absMax()))/I(h.lo);
    // A coefficient BOX gives a lower bound on the source speed. A box
    // containing zero is inconclusive, even when the curve is regular.
    I speed=Length(Range(st))/Sqr(I(w.hi));
    if(!upper.finite()||!speed.finite()||!(upper.hi<speed.lo))return INFINITY;
    // T'=S'+E': dot(S',T')>0 and tan(theta)<=M/(m-M).
    return (I(upper.hi)/(I(speed.lo)-I(upper.hi))).hi;
}
}
