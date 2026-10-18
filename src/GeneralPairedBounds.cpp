/**
 * @file GeneralPairedBounds.cpp
 * @brief 历史一般面片的双线性参数对应充分界。
 * 用三角分解与连续复合求界，UV 坍缩与三维极点须分别处理；不替代当前兜底验收合同。
 */
#include "GeneralInternal.h"
#include <functional>
namespace bspline2analytic::detail
{
namespace
{
BPoly Constant(I x) { return BPoly(0,0,x); }
// 两个角参数同时变化时，用积化和差公式控制多项式次数；区间常数包含八阶 Lagrange 余项。
BPoly Trig(const BPoly &angle, bool cosine)
{
    I range=angle.range();
    if (!range.finite()) return Constant(I::Whole());
    double center=range.lo+(range.hi-range.lo)*.5;
    BPoly d=angle-Constant(I(center)), power=Constant(I(1)), out;
    I sn=Sin(I(center)),cs=Cos(I(center)),factorial(1);
    for (unsigned k=0;k<8;++k)
    {
        I coefficient;
        switch ((k+unsigned(cosine))%4)
        { case 0: coefficient=sn;break;case 1:coefficient=cs;break;case 2:coefficient=-sn;break;default:coefficient=-cs; }
        out=out+power*(coefficient/factorial);
        if(k<7) power=power*d;
        factorial=factorial*I(k+1);
    }
    I remainder(1);
    for(unsigned k=1;k<=8;++k) remainder=remainder*I(d.range().absMax())/I(k);
    return out+Constant(I(-remainder.hi,remainder.hi));
}
BV Reference(const sggk::BSplineSurface &s,const BezierCell &cell,const GeneralTrim &trim)
{
    BPoly u(1,1),v(1,1);
    for(unsigned i=0;i<2;++i) for(unsigned j=0;j<2;++j)
    {
        I x=(I(i?cell.u1:cell.u0)-I(s.MinParamU()))/(I(s.MaxParamU())-I(s.MinParamU()));
        I y=(I(j?cell.v1:cell.v0)-I(s.MinParamV()))/(I(s.MaxParamV())-I(s.MinParamV()));
        std::array<I,4> w{(I(1)-x)*(I(1)-y),x*(I(1)-y),(I(1)-x)*y,x*y};
        for(unsigned k=0;k<4;++k) {u.at(i,j)=u.at(i,j)+w[k]*I(trim.corners[k].X());v.at(i,j)=v.at(i,j)+w[k]*I(trim.corners[k].Y());}
    }
    const auto &a=trim.parameters;
    BV p;
    switch(a.type)
    {
    case sggk::SurfaceType::Plane:p={u,v,Constant(I(0))};break;
    case sggk::SurfaceType::Cylinder:p={Trig(u,true)*I(a.radius),Trig(u,false)*I(a.radius),v};break;
    case sggk::SurfaceType::Cone:{I t=Sin(I(a.halfAngle))/Cos(I(a.halfAngle));p={v*Trig(u,true)*t,v*Trig(u,false)*t,v};break;}
    case sggk::SurfaceType::Sphere:
        p={(Trig(v+u,false)+Trig(v-u,false))*I(a.radius*.5),
           (Trig(v-u,true)-Trig(v+u,true))*I(a.radius*.5),Trig(v,true)*I(-a.radius)};break;
    default:
        p={Trig(u,true)*I(a.majorRadius)+(Trig(u+v,true)+Trig(u-v,true))*I(a.radius*.5),
           Trig(u,false)*I(a.majorRadius)+(Trig(u+v,false)+Trig(u-v,false))*I(a.radius*.5),Trig(v,false)*I(a.radius)};
    }
    IV z=Normalize(Exact(Get(a.axis))),x=Exact(Get(a.xDirection));
    x=Normalize(x-z*Dot(x,z));IV y=Cross(z,x),o=Exact(Get(a.origin));
    return {Constant(o.x)+p.x*x.x+p.y*y.x+p.z*z.x,
            Constant(o.y)+p.x*x.y+p.y*y.y+p.z*z.y,
            Constant(o.z)+p.x*x.z+p.y*y.z+p.z*z.z};
}
bool RegularQuad(const GeneralTrim &trim)
{
    int sign=0;
    for(unsigned u=0;u<2;++u) for(unsigned v=0;v<2;++v)
    {
        const auto &a=trim.corners[2*v],&b=trim.corners[2*v+1];
        const auto &c=trim.corners[u],&d=trim.corners[u+2];
        I det=(I(b.X())-I(a.X()))*(I(d.Y())-I(c.Y()))-(I(b.Y())-I(a.Y()))*(I(d.X())-I(c.X()));
        int current=det.lo>0?1:det.hi<0?-1:0;
        if(!current || (sign && current!=sign)) return false;
        sign=current;
    }
    return true;
}
bool CollapsedSideTriangle(const GeneralTrim &trim)
{
    const auto &corners=trim.corners;
    auto equal=[&](unsigned a,unsigned b) {
        return corners[a].X()==corners[b].X() && corners[a].Y()==corners[b].Y();
    };
    // 统计的是 UV 坍缩，不是空间极点或短边；不能将接近的 UV 强行吸附。
    constexpr unsigned ends[4][2]={{0,2},{1,3},{0,1},{2,3}};
    unsigned collapsed=0;
    for(const auto &edge:ends) collapsed+=equal(edge[0],edge[1]);
    if(collapsed!=1) return false;
    std::array<unsigned,3> unique{};
    unsigned count=0;
    for(unsigned i=0;i<4;++i)
    {
        bool seen=false;
        for(unsigned j=0;j<count;++j) seen|=equal(i,unique[j]);
        if(seen) continue;
        if(count==unique.size()) return false;
        unique[count++]=i;
    }
    if(count!=3) return false;
    const auto &a=corners[unique[0]],&b=corners[unique[1]],&c=corners[unique[2]];
    I ax=I(b.X())-I(a.X()),ay=I(b.Y())-I(a.Y());
    I bx=I(c.X())-I(a.X()),by=I(c.Y())-I(a.Y());
    I det=ax*by-ay*bx;
    // 忽略角点定向符号后，双线性映射的 Jacobian 为 det 乘 s、1-s、t 或 1-t；内部非零，只在坍缩边上为零。
    return det.lo>0 || det.hi<0;
}
}
bool TryPairedTrimBound(const sggk::BSplineSurface &s,const std::vector<BezierCell> &cells,
                       const GeneralTrim &trim,const ValidationOptions &options,GeneralValidationReport &report)
{
    // 这只是充分认证条件；一般源面重参数化可能不通过该检查，却能通过历史投影验证器。
    // 此说明不授权重新运行已废弃测试。
    if(cells.empty() || (!RegularQuad(trim) && !CollapsedSideTriangle(trim))) return false;
    double upper=0;
    unsigned visited=0;
    std::function<bool(const BezierCell &,unsigned)> visit=[&](const BezierCell &cell,unsigned depth) {
        ++visited;++report.piecesVisited;
        if(report.piecesVisited>options.maxValidationPieces || visited>512) return false;
        I weight=cell.w.range();
        if(weight.lo<=0) return false;
        IV error=Range(cell.p-Reference(s,cell,trim)*cell.w);
        I bound=Sqrt(Sqr(error.x)+Sqr(error.y)+Sqr(error.z))/weight;
        // 到最近母面投影对应 PCurve 的距离上界为 2d。
        if(bound.finite() && (I(2)*bound).hi<options.distanceTolerance)
        {upper=std::max(upper,bound.hi);return true;}
        if(depth>=std::min(8u,options.maxSubdivisionDepth)) return false;
        auto halves=Split(cell,depth%2==0);
        return visit(halves.first,depth+1)&&visit(halves.second,depth+1);
    };
    for(const auto &cell:cells) if(!visit(cell,0)) return false;
    report.geometry=report.boundaryHausdorff=report.pcurveLinearity=CheckState::Passed;
    report.pcurveFailure=PCurveFailure::None;report.hasWitness=false;
    report.distanceUpperBound=report.boundaryHausdorffUpperBound=upper;
    report.message="Continuous bilinear trim correspondence: surface and bidirectional boundary bounds; projected boundary deviation <= twice the bound";
    return true;
}
}
