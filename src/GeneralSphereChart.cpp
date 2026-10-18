/**
 * @file GeneralSphereChart.cpp
 * @brief 历史球面线性裁剪参数图的方向搜索。
 * 在固定母面上尝试不同球轴和参数图，处理周期及奇异性；当前有理边界路径采用自身恢复策略。
 */
#include "GeneralInternal.h"
namespace bspline2analytic::detail
{
namespace
{
using Points=std::array<std::vector<V>,4>;
struct ChartScore { AnalyticParameters chart; double cost=INFINITY,worst=INFINITY; };
AnalyticParameters Move(const AnalyticParameters &a,double dx,double dy)
{
    auto f=FrameOf(a);
    V z=Unit(f.z+f.x*dx+f.y*dy),x=Unit(f.x-z*Dot(f.x,z));
    auto b=a;b.axis=Direction(z);b.xDirection=Direction(x);return b;
}
std::vector<double> Residual(const AnalyticParameters &a,const Points &points)
{
    auto f=FrameOf(a);
    auto unit=a;unit.radius=1;
    std::vector<double> out;
    for(const auto &edge:points)
    {
        std::vector<std::array<double,2>> uv;
        std::vector<bool> regular;
        for(V q:edge)
        {
            V p{Dot(q,f.x),Dot(q,f.y),Dot(q,f.z)};
            auto t=MotherParameters(unit,p);
            if(!uv.empty()) t[0]=Unwrap(t[0],uv.back()[0]);
            uv.push_back(t);regular.push_back(std::hypot(p.x,p.y)>1e-7);
        }
        for(size_t end:{size_t(0),uv.size()-1}) if(!regular[end])
        {
            size_t first=uv.size(),second=uv.size();
            for(size_t k=0;k<uv.size();++k)
            {
                size_t j=end==0?k:uv.size()-1-k;
                if(regular[j]) {if(first==uv.size()) first=j;else {second=j;break;}}
            }
            if(second<uv.size() && std::abs(uv[second][1]-uv[first][1])>1e-12)
                uv[end][0]=uv[first][0]+(uv[second][0]-uv[first][0])*
                    (uv[end][1]-uv[first][1])/(uv[second][1]-uv[first][1]);
        }
        double du=uv.back()[0]-uv.front()[0],dv=uv.back()[1]-uv.front()[1],den=du*du+dv*dv;
        for(size_t i=0;i<uv.size();++i)
        {
            double t=den>1e-24?std::clamp(((uv[i][0]-uv.front()[0])*du+(uv[i][1]-uv.front()[1])*dv)/den,0.0,1.0):0;
            V p=CanonicalPoint(unit,uv.front()[0]+t*du,uv.front()[1]+t*dv);
            V q=f.x*p.x+f.y*p.y+f.z*p.z;
            V r=q-edge[i];out.push_back(r.x);out.push_back(r.y);out.push_back(r.z);
        }
    }
    return out;
}
ChartScore Score(const AnalyticParameters &a,const Points &points)
{
    auto residual=Residual(a,points);ChartScore s{a,0,0};
    for(size_t i=0;i<residual.size();i+=3)
    {
        double norm=residual[i]*residual[i]+residual[i+1]*residual[i+1]+residual[i+2]*residual[i+2];
        s.cost+=norm;s.worst=std::max(s.worst,std::sqrt(norm)*a.radius);
    }
    return s;
}
}
std::vector<AnalyticParameters> RefinedSphereCharts(SurfaceSampleContext &context,
    const GeneralConversionOptions &o,const AnalyticParameters &a,const std::vector<AnalyticParameters> &starts,unsigned &iterations)
{
    Points points;
    const auto &s=context.source;
    for(unsigned edge=0;edge<4;++edge) for(unsigned i=0;i<=16;++i)
    {
        double t=double(i)/16;
        double u=edge<2?(edge?s.MaxParamU():s.MinParamU()):s.MinParamU()+(s.MaxParamU()-s.MinParamU())*t;
        double v=edge<2?s.MinParamV()+(s.MaxParamV()-s.MinParamV())*t:(edge==2?s.MinParamV():s.MaxParamV());
        points[edge].push_back(Unit(context.At(u,v,sggk::DerivType::NONE,sggk::DerivType::NONE,false).point-Get(a.origin)));
    }
    std::vector<ChartScore> scores;
    auto add=[&](const AnalyticParameters &b) {
        for(const auto &old:scores) if(std::abs(Dot(Get(old.chart.axis),Get(b.axis)))>1-1e-10) return;
        scores.push_back(Score(b,points));
    };
    for(const auto &b:starts) add(b);
    // 通过确定性的方向覆盖搜索，避免依赖边界必须是圆或等参线。
    for(unsigned i=0;i<32;++i)
    {
        double z=1-2*(i+.5)/32,angle=i*2.399963229728653;
        V axis{std::sqrt(1-z*z)*std::cos(angle),std::sqrt(1-z*z)*std::sin(angle),z};
        auto b=a;b.axis=Direction(axis);b.xDirection=Direction(Perpendicular(axis));add(b);
    }
    std::stable_sort(scores.begin(),scores.end(),[](const auto &x,const auto &y){return x.cost<y.cost;});
    if(scores.size()>8) scores.resize(8);
    for(auto &score:scores)
    {
        double lambda=1e-8;
        for(unsigned iteration=0;iteration<60 && score.worst>o.fitting.validation.distanceTolerance*.005;++iteration)
        {
            ++iterations;
            auto r=Residual(score.chart,points);
            std::vector<std::vector<double>> jac(r.size(),std::vector<double>(2));
            constexpr double h=1e-5;
            for(unsigned j=0;j<2;++j)
            {
                auto plus=Residual(Move(score.chart,j?0:h,j?h:0),points);
                auto minus=Residual(Move(score.chart,j?0:-h,j?-h:0),points);
                for(size_t i=0;i<r.size();++i) jac[i][j]=(plus[i]-minus[i])/(2*h);
            }
            std::vector<double> rhs;for(double x:r) rhs.push_back(-x);
            for(unsigned j=0;j<2;++j) {std::vector<double> row(2);row[j]=std::sqrt(lambda);jac.push_back(row);rhs.push_back(0);}
            std::vector<double> step;
            if(!SolveLeastSquares(jac,rhs,step)) break;
            double length=std::hypot(step[0],step[1]);
            if(length>.35) for(auto &x:step) x*=.35/length;
            bool improved=false;
            for(unsigned line=0;line<8;++line)
            {
                auto trial=Score(Move(score.chart,step[0],step[1]),points);
                if(trial.cost<score.cost) {score=trial;improved=true;break;}
                for(auto &x:step) x*=.5;
            }
            lambda=improved?std::max(1e-14,lambda*.2):lambda*10;
            if(lambda>1e4) break;
        }
    }
    std::stable_sort(scores.begin(),scores.end(),[](const auto &x,const auto &y){return x.worst<y.worst;});
    std::vector<AnalyticParameters> out;for(const auto &score:scores) out.push_back(score.chart);
    return out;
}
}
