#include "GeneralRationalDataset.h"
#include <functional>
#include <sstream>

namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
double UpperNorm(IV p){return Sqrt(Sqr(p.x)+Sqr(p.y)+Sqr(p.z)).hi;}
double LowerNorm(IV p)
{
    auto lower=[](I a){return a.lo>0?a.lo:a.hi<0?-a.hi:0.;};
    return Sqrt(Sqr(I(lower(p.x)))+Sqr(I(lower(p.y)))+Sqr(I(lower(p.z)))).lo;
}
// Restrict AFTER differentiating. Tiny seam overlaps must not turn into
// ill-conditioned tiny-speed tangent comparisons.
BPoly RestrictAxis(const BPoly &p,bool u,I a,I b)
{
    if(a.lo<=0&&b.hi>=1)return p;
    unsigned n=u?p.u:p.v,m=u?p.v:p.u;BPoly out=p;
    I first(std::max(0.,a.lo)),last(std::min(1.,b.hi));
    for(unsigned k=0;k<=m;++k)
    {
        std::vector<I> row(n+1),right(n+1),left(n+1);
        for(unsigned j=0;j<=n;++j)row[j]=u?p.at(j,k):p.at(k,j);
        right[n]=row[n];
        for(unsigned r=1;r<=n;++r){for(unsigned j=0;j<=n-r;++j)row[j]=row[j]*(I(1)-first)+row[j+1]*first;right[n-r]=row[n-r];}
        I t=first.lo==1?I(0):(last-first)/(I(1)-first);left[0]=right[0];
        for(unsigned r=1;r<=n;++r){for(unsigned j=0;j<=n-r;++j)right[j]=right[j]*(I(1)-t)+right[j+1]*t;left[r]=right[0];}
        for(unsigned j=0;j<=n;++j)(u?out.at(j,k):out.at(k,j))=left[j];
    }
    return out;
}
struct Region
{
    I u0,u1,v0,v1;
    Region(const BezierCell &s,const std::array<double,4> &d)
    {
        auto range=[](double a,double b,double x){return a==b?I(0):(I(x)-I(a))/(I(b)-I(a));};
        u0=range(s.u0,s.u1,d[0]);u1=range(s.u0,s.u1,d[1]);v0=range(s.v0,s.v1,d[2]);v1=range(s.v0,s.v1,d[3]);
    }
    I operator()(const BPoly &p) const {return RestrictAxis(RestrictAxis(p,true,u0,u1),false,v0,v1).range();}
    IV operator()(const BV &p) const {return {(*this)(p.x),(*this)(p.y),(*this)(p.z)};}
};
double Position(const BezierCell &s,const RationalReferencePiece &piece)
{
    const auto &t=piece.target;Region range(s,piece.domain);
    auto h=s.w*t.w;if(range(h).lo<=0)return INFINITY;
    return (I(UpperNorm(range(s.p*t.w-t.p*s.w)))/I(range(h).lo)).hi;
}
double TangentRatio(const BezierCell &s,const RationalReferencePiece &piece,bool u)
{
    const auto &t=piece.target;Region range(s,piece.domain);
    auto h=s.w*t.w;auto g=t.p*s.w-s.p*t.w;
    I hr=range(h);if(hr.lo<=0)return INFINITY;
    auto numerator=Derivative(g,u)*h-g*Derivative(h,u);
    I error=I(UpperNorm(range(numerator)))/Sqr(I(hr.lo));
    auto j=Derivative(s.p,u)*s.w-s.p*Derivative(s.w,u);
    I speed=I(LowerNorm(range(j)))/Sqr(I(range(s.w).hi));
    if(speed.lo<=error.hi)return INFINITY;
    return (error/(I(speed.lo)-error)).hi;
}
bool ValidMapping(const GeneralFixture &c)
{
    if(c.bend)
    {
        if(c.group!=GeneralFixtureGroup::NonlinearPCurve)return false;
        BPoly x(1,1),y(1,1),bend(2,0);
        for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)
        {x.at(i,j)=I(c.corners[i+2*j].X());y.at(i,j)=I(c.corners[i+2*j].Y());}
        bend.at(1,0)=I(c.bend)/I(2);y=y+bend;
        auto jacobian=Derivative(x,true)*Derivative(y,false)-Derivative(x,false)*Derivative(y,true);
        int sign=0;
        std::function<bool(const BPoly&,unsigned)> prove=[&](const BPoly &p,unsigned depth)
        {
            I det=p.range();int current=det.lo>0?1:det.hi<0?-1:0;
            if(current){if(sign&&current!=sign)return false;sign=current;return true;}
            if(depth==8)return false;auto halves=Split(p,depth%2==0);
            return prove(halves.first,depth+1)&&prove(halves.second,depth+1);
        };
        return prove(jacobian,0);
    }
    int sign=0;unsigned zero=0;
    for(double u:{0.,1.})for(double v:{0.,1.})
    {
        auto p=c.corners;I ux=(I(p[1].X())-I(p[0].X()))*I(1-v)+(I(p[3].X())-I(p[2].X()))*I(v);
        I uy=(I(p[1].Y())-I(p[0].Y()))*I(1-v)+(I(p[3].Y())-I(p[2].Y()))*I(v);
        I vx=(I(p[2].X())-I(p[0].X()))*I(1-u)+(I(p[3].X())-I(p[1].X()))*I(u);
        I vy=(I(p[2].Y())-I(p[0].Y()))*I(1-u)+(I(p[3].Y())-I(p[1].Y()))*I(u);
        I det=ux*vy-uy*vx;int current=det.lo>0?1:det.hi<0?-1:0;
        if(!current)++zero;else if(sign&&sign!=current)return false;else sign=current;
    }
    if(!sign)return false;if(!zero)return true;
    // Only a genuine triangle with exactly one collapsed parameter side.
    unsigned collapsed=0;const unsigned ends[][2]={{0,2},{1,3},{0,1},{2,3}};
    for(auto &e:ends)collapsed+=c.corners[e[0]].X()==c.corners[e[1]].X()&&c.corners[e[0]].Y()==c.corners[e[1]].Y();
    return zero==2&&collapsed==1;
}
}
RationalSourceCertificate CertifyRationalDomain(const GeneralFixture &c,unsigned maxPieces,bool boundariesOnly)
{
    RationalSourceCertificate out;
    if(!c.rationalDomain||!c.input||!IntervalEnvironmentSupported()) {out.message="V2 source/environment required";return out;}
    if(c.input->MinParamU()!=0||c.input->MaxParamU()!=1||c.input->MinParamV()!=0||c.input->MaxParamV()!=1)
    {out.message="Source parameter domain must be the unit square";return out;}
    bool negative=c.group==GeneralFixtureGroup::Infeasible;
    if(negative)
    {auto n=CertifyNoQuartic(*c.input,c.exclusionRadius);static_cast<FixtureCertificate&>(out)=n;if(!n.passed)return out;out.passed=false;}
    if(!ValidMapping(c)){out.message="Invalid rational-domain mapping";return out;}
    for(auto q:c.corners)if(q.X()<0||q.X()>1||q.Y()<0||q.Y()>1){out.message="Corner outside reference chart";return out;}
    auto cp=c.input->ControlPoints();V anchor=Get(cp.front().front());double shiftError=0;
    for(auto &row:cp)for(auto &p:row)
    {
        V x=Get(p)-anchor;shiftError=std::max(shiftError,UpperNorm(Exact(Get(p))-Exact(anchor)-Exact(x)));p=Point(x);
    }
    const auto &s=*c.input;sggk::BSplineSurfacePtr shifted;
    if(s.IsRational())shifted=std::make_shared<sggk::BSplineSurface>(s.DegreeU(),s.DegreeV(),cp,s.Weights(),s.KnotsU(),s.KnotsV(),s.MultsU(),s.MultsV());
    else shifted=std::make_shared<sggk::BSplineSurface>(s.DegreeU(),s.DegreeV(),cp,s.KnotsU(),s.KnotsV(),s.MultsU(),s.MultsV());
    auto reference=c;if(reference.type==5)reference.type=0;
    RationalReference ref(reference);out.upper=0;double maxRatio=0;
    // Interval extraction from the ACTUAL stored controls, including recentering.
    auto cells=ExtractCells(*shifted);
    for(auto &cell:cells)
    {
        I e(-shiftError,shiftError);cell.p.x=cell.p.x+cell.w*e;cell.p.y=cell.p.y+cell.w*e;cell.p.z=cell.p.z+cell.w*e;
    }
    const double ratioLimit=(Sin(I(.01)*IPi()/I(180))/Cos(I(.01)*IPi()/I(180))).lo;
    std::function<bool(const BezierCell&,unsigned)> visit=[&](const BezierCell &cell,unsigned depth)
    {
        if(++out.pieces>maxPieces){out.message="Rational source certificate budget exhausted";return false;}
        auto targets=ref.Cover(cell,anchor);double bound=0;
        for(const auto &t:targets)bound=std::max(bound,Position(cell,t));
        if(bound<c.pairTolerance){out.upper=std::max(out.upper,bound);return true;}
        if(depth>=24){std::ostringstream m;m<<"Source distance unresolved: "<<bound<<" tolerance="<<c.pairTolerance;out.message=m.str();return false;}
        bool u=(cell.u1-cell.u0)>=(cell.v1-cell.v0);
        auto halves=Split(cell,u);return visit(halves.first,depth+1)&&visit(halves.second,depth+1);
    };
    if(!negative&&!boundariesOnly)for(const auto &cell:cells)if(!visit(cell,0))return out;
    // Separate 1D proof avoids evaluating boundary derivatives on 2D interiors.
    for(unsigned edge=0;edge<4;++edge)
    {
        if(c.collapsedMask&(1u<<edge))continue;bool fixedU=edge<2,high=edge%2;
        std::function<bool(const BezierCell&,unsigned)> tangent=[&](const BezierCell &cell,unsigned depth)
        {
            if(++out.pieces>maxPieces){out.message="Source tangent budget exhausted";return false;}
            auto ts=ref.Cover(cell,anchor);double ratio=0;
            for(const auto &t:ts)ratio=std::max(ratio,TangentRatio(cell,t,!fixedU));
            if(ratio<ratioLimit){maxRatio=std::max(maxRatio,ratio);return true;}
            if(depth>=24){out.message="Source boundary tangent unresolved";return false;}
            auto halves=Split(cell,!fixedU);return tangent(halves.first,depth+1)&&tangent(halves.second,depth+1);
        };
        for(const auto &cell:cells)
        {
            if((fixedU?(high?cell.u1:cell.u0):(high?cell.v1:cell.v0))!=(high?1.:0.))continue;
            auto side=BoundaryCell(cell,static_cast<BoundaryId>(edge));
            if((negative||boundariesOnly)&&!visit(side,0))return out;
            if(!tangent(side,0))return out;
        }
    }
    out.tangentDeg=(Atan(I(maxRatio))*I(180)/IPi()).hi;out.passed=true;
    out.boundaryUpper=out.upper;if(negative)out.upper=INFINITY;
    out.message="Independent rational-domain source position and boundary-tangent certificate";return out;
}
}
