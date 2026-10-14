#include "GeneralInternal.h"
namespace bspline2analytic
{
const char *GeneralSeedFailureName(GeneralSeedFailure failure)
{
    switch (failure)
    {
    case GeneralSeedFailure::None: return "None";
    case GeneralSeedFailure::InsufficientPoints: return "InsufficientPoints";
    case GeneralSeedFailure::RankDeficient: return "RankDeficient";
    case GeneralSeedFailure::AmbiguousAxis: return "AmbiguousAxis";
    case GeneralSeedFailure::InvalidParameters: return "InvalidParameters";
    case GeneralSeedFailure::NonFinite: return "NonFinite";
    default: return "ConstructionFailed";
    }
}
}
namespace bspline2analytic::detail
{
namespace
{
using Row = std::array<double,3>;
using Rows = std::vector<Row>;
constexpr double Eps = std::numeric_limits<double>::epsilon();
double RankFloor(size_t rows) { return 64 * Eps * double(std::max(size_t(3),rows)); }
void Reject(GeneralSeedGroupDiagnostic &d, GeneralSeedFailure why)
{
    d.failure=why;
    throw Failure(Status::IllConditioned,GeneralSeedFailureName(why));
}
struct SVD3
{
    std::array<double,3> singular{}; // Descending.
    std::array<V,3> right{};
};
// One-sided Jacobi directly on the n x 3 matrix. Never form A^T A: doing so
// would square the conditioning of the short-arc direction/normal plane fit.
bool Decompose(Rows a,SVD3 &out)
{
    if(a.size()<3) return false;
    double magnitude=0;
    for(const auto &r:a) for(double x:r)
    {
        if(!std::isfinite(x)) return false;
        magnitude=std::max(magnitude,std::abs(x));
    }
    if(!(magnitude>0)) return false;
    for(auto &r:a) for(double &x:r) x/=magnitude;
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    bool converged=false;
    for(unsigned sweep=0;sweep<48;++sweep)
    {
        bool changed=false;
        for(unsigned p=0;p<2;++p) for(unsigned q=p+1;q<3;++q)
        {
            double alpha=0,beta=0,gamma=0;
            for(const auto &r:a)
            {alpha+=r[p]*r[p];beta+=r[q]*r[q];gamma+=r[p]*r[q];}
            if(std::abs(gamma)<=16*Eps*std::sqrt(alpha)*std::sqrt(beta)) continue;
            double tau=(beta-alpha)/(2*gamma);
            double t=std::copysign(1.0,tau)/(std::abs(tau)+std::hypot(1.0,tau));
            double c=1/std::hypot(1.0,t),s=c*t;
            for(auto &r:a)
            {double x=r[p],y=r[q];r[p]=c*x-s*y;r[q]=s*x+c*y;}
            for(auto &r:v)
            {double x=r[p],y=r[q];r[p]=c*x-s*y;r[q]=s*x+c*y;}
            changed=true;
        }
        if(!changed) {converged=true;break;}
    }
    if(!converged) return false;
    std::array<double,3> lengths{};
    for(const auto &r:a) for(unsigned k=0;k<3;++k) lengths[k]=std::hypot(lengths[k],r[k]);
    std::array<unsigned,3> order{0,1,2};
    std::sort(order.begin(),order.end(),[&](unsigned i,unsigned j){return lengths[i]>lengths[j];});
    for(unsigned k=0;k<3;++k)
    {
        unsigned j=order[k];out.singular[k]=lengths[j]*magnitude;
        out.right[k]=Unit({v[0][j],v[1][j],v[2][j]});
    }
    return out.singular[0]>0 && std::isfinite(out.singular[0]);
}
V Solve(const Rows &rows,const std::vector<double> &rhs,GeneralSeedGroupDiagnostic &d)
{
    SVD3 svd;
    if(!Decompose(rows,svd)) Reject(d,GeneralSeedFailure::RankDeficient);
    double rank=svd.singular[2]/svd.singular[0];
    d.relativeRank=std::min(d.relativeRank,rank);
    std::array<double,3> x;
    if(!(rank>RankFloor(rows.size())) || !LeastSquares3(rows,rhs,RankFloor(rows.size()),x))
        Reject(d,GeneralSeedFailure::RankDeficient);
    V result{x[0],x[1],x[2]};
    if(!Finite(result)) Reject(d,GeneralSeedFailure::NonFinite);
    for(size_t i=0;i<rows.size();++i)
    {
        double residual=std::abs(Dot({rows[i][0],rows[i][1],rows[i][2]},result)-rhs[i]);
        if(!std::isfinite(residual)) Reject(d,GeneralSeedFailure::NonFinite);
        d.constructionResidual=std::max(d.constructionResidual,residual);
    }
    return result;
}
V PlaneNormal(const Rows &rows,GeneralSeedGroupDiagnostic &d)
{
    SVD3 svd;
    if(!Decompose(rows,svd)) Reject(d,GeneralSeedFailure::AmbiguousAxis);
    double rank=svd.singular[1]/svd.singular[0];
    double gap=(svd.singular[1]-svd.singular[2])/svd.singular[0];
    d.relativeRank=std::min(d.relativeRank,rank);d.axisSeparation=gap;
    if(!(rank>RankFloor(rows.size()) && gap>8*RankFloor(rows.size())))
        Reject(d,GeneralSeedFailure::AmbiguousAxis);
    for(const auto &r:rows)
        d.constructionResidual=std::max(d.constructionResidual,
            std::abs(Dot({r[0],r[1],r[2]},svd.right[2])));
    return svd.right[2];
}
}
AnalyticParameters ConditionedGeneralSeed(sggk::SurfaceType type,
    const std::vector<SurfaceSample> &points,double scale,GeneralSeedGroupDiagnostic &diagnostic)
{
    auto &d=diagnostic;
    d.sampleCount=static_cast<unsigned>(points.size());
    if(points.size()<4) Reject(d,GeneralSeedFailure::InsufficientPoints);
    if(!(scale>0) || !std::isfinite(scale)) Reject(d,GeneralSeedFailure::NonFinite);
    // Subtract a nearby anchor before averaging, including for large world coordinates.
    V anchor=points.front().point,mean{};
    for(const auto &p:points)
    {
        if(!Finite(p.point) || !Finite(p.normal) || !p.normalValid)
            Reject(d,GeneralSeedFailure::NonFinite);
        mean=mean+(p.point-anchor)/scale/double(points.size());
    }
    std::vector<V> local,normals;
    for(const auto &p:points)
    {local.push_back((p.point-anchor)/scale-mean);normals.push_back(Unit(p.normal));}
    AnalyticParameters a;a.type=type;
    V center{},axis{0,0,1};
    Rows rows;std::vector<double> rhs;
    if(type==sggk::SurfaceType::Sphere)
    {
        for(size_t i=0;i<points.size();++i)
        {
            V x=Perpendicular(normals[i]),y=Cross(normals[i],x);
            for(V v:{x,y}) {rows.push_back({v.x,v.y,v.z});rhs.push_back(Dot(v,local[i]));}
        }
        center=Solve(rows,rhs,d);
        double radius=0;
        for(V p:local) radius+=Norm(p-center)/double(points.size());
        a.radius=radius*scale;
    }
    else if(type==sggk::SurfaceType::Cylinder)
    {
        for(V n:normals) rows.push_back({n.x,n.y,n.z});
        axis=PlaneNormal(rows,d);
        V x=Perpendicular(axis),y=Cross(axis,x);rows.clear();
        for(size_t i=0;i<points.size();++i)
        {
            V n=normals[i]-axis*Dot(normals[i],axis);
            if(!(Norm(n)>RankFloor(points.size()))) Reject(d,GeneralSeedFailure::AmbiguousAxis);
            n=Unit(n);
            rows.push_back({1,0,Dot(n,x)});rhs.push_back(Dot(local[i],x));
            rows.push_back({0,1,Dot(n,y)});rhs.push_back(Dot(local[i],y));
        }
        V circle=Solve(rows,rhs,d);
        center=x*circle.x+y*circle.y;a.radius=std::abs(circle.z)*scale;
    }
    else if(type==sggk::SurfaceType::Cone)
    {
        for(size_t i=0;i<points.size();++i)
        {V n=normals[i];rows.push_back({n.x,n.y,n.z});rhs.push_back(Dot(n,local[i]));}
        center=Solve(rows,rhs,d);
        std::vector<V> directions;
        for(V p:local)
        {
            if(!(Norm(p-center)>RankFloor(points.size()))) Reject(d,GeneralSeedFailure::InvalidParameters);
            directions.push_back(Unit(p-center));
        }
        V first=directions.front(),offset{};
        for(V v:directions) offset=offset+(v-first)/double(directions.size());
        rows.clear();
        for(V v:directions) {v=(v-first)-offset;rows.push_back({v.x,v.y,v.z});}
        axis=PlaneNormal(rows,d);
        if(Dot(axis,first+offset)<0) axis=axis*-1;
        double sine=0,cosine=0;
        for(V v:directions) {sine+=Norm(Cross(axis,v));cosine+=Dot(axis,v);}
        a.halfAngle=std::atan2(sine,cosine);
        if(!(a.halfAngle>0 && a.halfAngle<Pi/2)) Reject(d,GeneralSeedFailure::InvalidParameters);
    }
    else Reject(d,GeneralSeedFailure::InvalidParameters);
    V world=anchor+(mean+center)*scale;
    if(!Finite(world) || !Finite(axis) || !std::isfinite(a.radius)) Reject(d,GeneralSeedFailure::NonFinite);
    if(type!=sggk::SurfaceType::Cone && !(a.radius>0)) Reject(d,GeneralSeedFailure::InvalidParameters);
    a.origin=Point(world);a.axis=Direction(axis);a.xDirection=Direction(Perpendicular(axis));
    return a;
}
}
