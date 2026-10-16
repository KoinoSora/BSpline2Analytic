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
    // The stored doubles DEFINE lambda exactly. Rounding in this heuristic
    // cannot invalidate the proof; the complete residual is enclosed below.
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
    // Division by the largest component avoids overflow of 1/size for tiny e.
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
        // Convexity bounds the vector itself: no squared-norm spline products.
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
    // E<||A|| and lambda>0 prove B!=0 and A.B>0, ruling out antiparallel fields.
    return q.finite()&&q.hi<1?std::max(0.,q.hi):INFINITY;
}
}
