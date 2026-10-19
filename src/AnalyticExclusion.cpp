/** @file AnalyticExclusion.cpp
 * @brief 用 14 维隐式多项式空间及经区间验证的左逆证明五类解析母面均超差。
 * QR 只生成证书候选；所有排除结论必须由向外舍入的残差和距离下界支持。
 */
#include "AnalyticExclusionInternal.h"
#include "IntervalBounds.h"
#include <numeric>

namespace bspline2analytic::detail
{
namespace
{
constexpr unsigned N = 14;
using Row = std::array<double,N>;
using IRow = std::array<I,N>;
struct LocalLimit {};
struct ClockGuard
{
    AnalyticExclusionReport &report;
    const ExecutionDeadline &deadline;
    double limit;
    ExecutionDeadline::Clock::time_point start = ExecutionDeadline::Clock::now();
    double Ms() const { return std::chrono::duration<double,std::milli>(ExecutionDeadline::Clock::now()-start).count(); }
    void Check() const
    {
        deadline.Check("fallback_exclusion");
        if (Ms() >= limit) throw LocalLimit{};
    }
    ~ClockGuard() { report.totalMs = Ms(); }
};
IRow Features(IV p)
{
    I x=p.x,y=p.y,z=p.z,xx=Sqr(x),yy=Sqr(y),zz=Sqr(z),r=xx+yy+zz;
    return {I(1),x,y,z,xx,x*y,x*z,yy,y*z,zz,x*r,y*r,z*r,Sqr(r)};
}
// 带列主元、二次正交化的 QR。病态或秩不足时保守退出；不以浮点秩估计作排除。
bool LeftInverse(const std::vector<IRow> &a, std::vector<Row> &inverse, const ClockGuard &clock)
{
    const size_t m=a.size();
    std::array<std::vector<double>,N> v,q;
    std::array<Row,N> r{};
    std::array<unsigned,N> order{};
    std::iota(order.begin(),order.end(),0u);
    for (unsigned j=0;j<N;++j)
    {
        v[j].resize(m);q[j].resize(m);
        for(size_t i=0;i<m;++i)v[j][i]=a[i][j].lo*.5+a[i][j].hi*.5;
    }
    double initialNorm=0;
    for(unsigned k=0;k<N;++k)
    {
        clock.Check();
        unsigned pivot=k;double norm2=-1;
        for(unsigned j=k;j<N;++j)
        {
            double t=std::inner_product(v[j].begin(),v[j].end(),v[j].begin(),0.);
            if(t>norm2){norm2=t;pivot=j;}
        }
        if(k==0)initialNorm=std::sqrt(norm2);
        double norm=std::sqrt(norm2);
        if(!std::isfinite(norm)||!(norm>initialNorm*1e-12))return false;
        std::swap(v[k],v[pivot]);std::swap(order[k],order[pivot]);
        for(unsigned i=0;i<k;++i)std::swap(r[i][k],r[i][pivot]);
        r[k][k]=norm;
        for(size_t i=0;i<m;++i)q[k][i]=v[k][i]/norm;
        for(unsigned j=k+1;j<N;++j)
            for(unsigned pass=0;pass<2;++pass)
            {
                double t=std::inner_product(q[k].begin(),q[k].end(),v[j].begin(),0.);
                r[k][j]+=t;
                for(size_t i=0;i<m;++i)v[j][i]-=t*q[k][i];
            }
    }
    // inverse[i][j] 表示 B 的第 j 行第 i 列，方便连续访问每个见证。
    inverse.resize(m);
    for(size_t i=0;i<m;++i)
    {
        Row x{};
        for(int j=N-1;j>=0;--j)
        {
            double t=q[j][i];
            for(unsigned k=j+1;k<N;++k)t-=r[j][k]*x[k];
            x[j]=t/r[j][j];
            if(!std::isfinite(x[j]))return false;
            inverse[i][order[j]]=x[j];
        }
    }
    return true;
}
void Compute(const sggk::BSplineSurface &s,double tolerance,AnalyticExclusionReport &report,const ClockGuard &clock)
{
    if(!IntervalEnvironmentSupported()||!(tolerance>0&&std::isfinite(tolerance)))
    {report.message="Unsupported interval environment or tolerance";return;}
    const auto &cp=s.ControlPoints();
    // 限制不可分割的提取工作量；超出支持规模时交回原优化流程。
    if(cp.empty()||cp.front().empty()||cp.size()*cp.front().size()>512||s.DegreeU()>5||s.DegreeV()>5)
    {report.message="Exclusion complexity limit";return;}
    auto spans=[](const auto &knots,double lo,double hi) {
        unsigned n=0;
        for(size_t i=1;i<knots.size();++i)
            if(knots[i-1]<knots[i]&&knots[i]>lo&&knots[i-1]<hi)++n;
        return n;
    };
    if(size_t(spans(s.KnotsU(),s.MinParamU(),s.MaxParamU()))*
       spans(s.KnotsV(),s.MinParamV(),s.MaxParamV())>32)
    {report.message="Exclusion cell limit";return;}
    std::array<double,3> lo{INFINITY,INFINITY,INFINITY},hi{-INFINITY,-INFINITY,-INFINITY},center{};
    for(const auto &row:cp)for(const auto &p:row)
    {
        const std::array<double,3> xyz{p.X(),p.Y(),p.Z()};
        for(unsigned j=0;j<3;++j)
        {
            if(!std::isfinite(xyz[j])){report.message="Nonfinite control point";return;}
            lo[j]=std::min(lo[j],xyz[j]);hi[j]=std::max(hi[j],xyz[j]);
        }
    }
    double extent=0;
    for(unsigned j=0;j<3;++j){center[j]=lo[j]*.5+hi[j]*.5;extent=std::max(extent,hi[j]-lo[j]);}
    if(!(extent>0&&std::isfinite(extent))){report.message="Degenerate normalization";return;}
    int exponent=0;std::frexp(extent,&exponent);
    const double scale=std::ldexp(1.,exponent-1);
    if(!(scale>0&&std::isfinite(scale))){report.message="Unsupported normalization";return;}
    report.scale=scale;
    clock.Check();
    auto cells=ExtractCells(s);
    clock.Check();
    if(cells.empty()||cells.size()>32){report.message="Exclusion cell limit";return;}
    std::vector<IRow> a;a.reserve(cells.size()*25);
    double box=0;
    for(const auto &cell:cells)
        for(unsigned u=0;u<=4;++u)for(unsigned v=0;v<=4;++v)
        {
            clock.Check();
            I w=Evaluate(cell.w,u*.25,v*.25).c[0];
            if(!w.finite()||w.lo<=0){report.message="Uncertified rational denominator";return;}
            IV p{(Evaluate(cell.p.x,u*.25,v*.25).c[0]/w-I(center[0]))/I(scale),
                 (Evaluate(cell.p.y,u*.25,v*.25).c[0]/w-I(center[1]))/I(scale),
                 (Evaluate(cell.p.z,u*.25,v*.25).c[0]/w-I(center[2]))/I(scale)};
            auto row=Features(p);
            for(I value:row)if(!value.finite()){report.message="Nonfinite feature interval";return;}
            box=std::max({box,p.x.absMax(),p.y.absMax(),p.z.absMax()});
            a.push_back(row);
        }
    report.witnesses=static_cast<unsigned>(a.size());
    std::vector<Row> b;
    if(a.size()<N||!LeftInverse(a,b,clock))
    {report.message="Rank deficient or ill-conditioned feature matrix";return;}
    double rho=0,norm=0;
    for(unsigned j=0;j<N;++j)
    {
        clock.Check();
        I rowResidual(0),rowNorm(0);
        for(size_t i=0;i<a.size();++i)rowNorm=rowNorm+I(std::abs(b[i][j]));
        for(unsigned k=0;k<N;++k)
        {
            I dot(0);
            for(size_t i=0;i<a.size();++i)dot=dot+I(b[i][j])*a[i][k];
            rowResidual=rowResidual+I((I(j==k?1.:0.)-dot).absMax());
        }
        rho=std::max(rho,rowResidual.hi);norm=std::max(norm,rowNorm.hi);
    }
    report.residualNormUpper=rho;report.inverseNormUpper=norm;
    if(!(rho<1&&norm>0&&std::isfinite(norm)))
    {report.message="Left inverse residual not certified below one";return;}
    // 任一系数无穷范数为 1 的隐式多项式，在离见证至多 1 的参数盒内的梯度 L1 上界。
    I boxBound=I(box)+I(1);
    I lipschitz=I(3)+I(12)*boxBound+I(27)*Sqr(boxBound)+I(36)*Sqr(boxBound)*boxBound;
    report.gradientBound=lipschitz.hi;
    I lower=(I(1)-I(rho))/(I(norm)*I(lipschitz.hi));
    lower=I(scale)*I(std::min(1.,lower.lo));
    if(lower.finite())report.distanceLowerBound=std::max(0.,lower.lo);
    clock.Check();
    if(report.distanceLowerBound>tolerance)
    {report.status=AnalyticExclusionStatus::CertifiedExcluded;report.message="All five analytic mother types excluded by certified distance lower bound";}
    else report.message="Certified distance lower bound does not exceed tolerance";
}
}
void CertifyAnalyticExclusion(const sggk::BSplineSurface &s,double tolerance,const ExecutionDeadline &deadline,
                             AnalyticExclusionReport &report,double localBudgetMs)
{
    report={};report.attempts=1;report.status=AnalyticExclusionStatus::Inconclusive;
    ClockGuard clock{report,deadline,localBudgetMs};
    try
    {
        if(!(localBudgetMs>0&&std::isfinite(localBudgetMs)))throw LocalLimit{};
        Compute(s,tolerance,report,clock);
    }
    catch(const LocalLimit &)
    {report.status=AnalyticExclusionStatus::LocalBudgetExceeded;report.message="Exclusion local budget exhausted; continue original search";}
    catch(const DeadlineExceeded &){throw;}
    catch(const std::exception &e)
    {report.status=AnalyticExclusionStatus::Inconclusive;report.message=std::string("Exclusion inconclusive: ")+e.what();}
}
}
