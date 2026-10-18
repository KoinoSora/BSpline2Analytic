/**
 * @file GeneralRationalBoundaryValidation.cpp
 * @brief 当前一般兜底的完整验证流程。
 * 内部采样筛查、图恢复、SDK 边界构造及连续位置/切向认证共享预算；未通过或超时不能发布成功输出。
 */
#include "GeneralRationalBoundaryInternal.h"
#include <chrono>

namespace bspline2analytic::detail::rational_boundary
{
namespace
{
using Clock=std::chrono::steady_clock;
double Ms(Clock::time_point t){return std::chrono::duration<double,std::milli>(Clock::now()-t).count();}
struct Measure{double &out;Clock::time_point start=Clock::now();explicit Measure(double &x):out(x){}~Measure(){out+=Ms(start);}};
struct Budget:std::runtime_error{using std::runtime_error::runtime_error;};
class Certifier
{
    SurfaceSampleContext &samples;const Trim &trim;const Target &target;
    const GeneralRationalBoundaryOptions &options;GeneralRationalBoundaryReport &report;
    paired::Products products;CurveErrors errors;CurveEnvelopes envelopes;
    double distance=0,ratio=0,sineSquared=0,limit;I sineSquaredLimit;
    std::array<bool,4> envelopeOnEdge{true,true,true,true};
    std::vector<GeneralBoundaryTangentExemption> exemptions;
    bool MergeMicro(const BezierCell &s,unsigned edge,double lo,double hi)
    {
        auto proof=MicroSegmentBound(s,target.edges[edge],lo,hi,options.microSegmentDiameter,
                                    options.general.fitting.validation.distanceTolerance);
        if(!proof)return false;
        proof->edge=edge;exemptions.push_back(*proof);++report.microMergedSegments;++report.boundarySegments;
        distance=std::max(distance,proof->distanceUpper);return true;
    }
    double ScaledTangent(const BV &s,const BV &t,std::optional<BPoly> &scale,bool envelope)
    {
        ++report.tangentResidualChecks;
        if(!envelope)++report.tangentResidualIntervalFallbacks;
        if(!scale)scale=ChooseTangentScale(s,t);
        if(!scale)return INFINITY;
        auto residual=envelope?envelopes.TangentResidual(s,t,*scale):
            std::optional<BV>(BuildScaledTangentResidual(s,t,*scale,products));
        double q=residual?ScaledTangentSineSquaredUpper(s,*residual,*scale):INFINITY;
        if(q<sineSquaredLimit.lo)++report.tangentResidualAccepted;
        return q;
    }
    bool Accept(double bound,double angleRatio,double angleSquared)
    {
        if(!(bound<options.general.fitting.validation.distanceTolerance))return false;
        // 保留原快速判据，只有新的后备残差路径使用正弦平方界；反三角函数留到最终报告组装时计算。
        if(angleRatio<limit)ratio=std::max(ratio,angleRatio);
        else if(angleSquared<sineSquaredLimit.lo)sineSquared=std::max(sineSquared,angleSquared);
        else return false;
        distance=std::max(distance,bound);++report.boundarySegments;return true;
    }
    // 节点状态：源段 s 的原参数域对应目标 [lo,hi]。先对齐目标接缝，再依次尝试
    // 双精度保守包络、完整区间位置界、误差导数/正比例切向界；这些路径共享节点预算。
    // 未获通过时尝试合法微段豁免，否则细分并可修正中点对应。修正后残差必须重建。
    bool Visit(const BezierCell &s,unsigned edge,double lo,double hi,unsigned depth)
    {
        samples.CheckDeadline("boundary_certificate");
        if(++report.boundaryNodes>options.maxBoundaryNodes)throw Budget("Rational boundary node budget exhausted");
        if(!(lo<hi))return lo==hi&&MergeMicro(s,edge,lo,hi);
        for(const auto &part:target.edges[edge])if(part.u1>lo&&part.u1<hi)
        {
            double cut=s.u0+(s.u1-s.u0)*(part.u1-lo)/(hi-lo);
            if(!(cut>s.u0&&cut<s.u1))
            {
                // 完整保留目标接缝片段。
                // 先通过空间微小段证书将其与源端点配对，再按普通规则认证剩余对应区间，不能删除接缝附近的覆盖。
                if(cut<=s.u0&&MergeMicro(paired::Crop(s,s.u0,s.u0,0,1),edge,lo,part.u1))
                    return Visit(s,edge,part.u1,hi,depth);
                if(cut>=s.u1&&MergeMicro(paired::Crop(s,s.u1,s.u1,0,1),edge,part.u1,hi))
                    return Visit(s,edge,lo,part.u1,depth);
                return false;
            }
            return Visit(paired::Crop(s,s.u0,cut,0,1),edge,lo,part.u1,depth)&&
                   Visit(paired::Crop(s,cut,s.u1,0,1),edge,part.u1,hi,depth);
        }
        const BezierCell *part=nullptr;
        for(const auto &p:target.edges[edge])if(p.u0<=lo&&hi<=p.u1){part=&p;break;}
        if(!part)return false;
        auto t=paired::Crop(*part,lo,hi,0,1);
        // 辅助比例只对当前精确对应关系有效；子节点可能改变映射端点，必须重新选择 lambda 并重建残差。
        std::optional<BPoly> tangentScale;
        // 对同一源/目标段先用双精度系数和保守舍入包络过滤。
        // 无结论不额外消耗节点，也不改变对应关系或验收阈值；纯多项式差没有可节省的有理卷积，因此跳过此过滤。
        bool useEnvelope=envelopeOnEdge[edge]&&(s.w.u||t.w.u);
        bool positionExceeded=false,tangentExceeded=false,fastPositionPassed=false;
        if(useEnvelope)++report.envelopeChecks;
        if(auto fast=useEnvelope?envelopes.Build(s,t):std::nullopt)
        {
            double bound=paired::PositionBoxUpper(fast->numerator,fast->denominator);
            if(bound<options.general.fitting.validation.distanceTolerance)
            {
                fastPositionPassed=true;
                double angleRatio=trim.collapsed[edge]?0:INFINITY,angleSquared=INFINITY;
                if(!trim.collapsed[edge])
                    if(auto st=envelopes.Tangent(s))
                    {
                        angleRatio=ResidualTangentRatio(*fast,s,*st,bound);
                        if(!(angleRatio<limit))
                            if(auto tt=envelopes.Tangent(t))
                            {
                                tangentExceeded=envelopes.TangentExceeded(*st,*tt,sineSquaredLimit.hi);
                                if(!tangentExceeded)angleSquared=ScaledTangent(*st,*tt,tangentScale,true);
                            }
                    }
                if(Accept(bound,angleRatio,angleSquared))
                {
                    ++report.envelopeAccepted;return true;
                }
            }
            else positionExceeded=envelopes.PositionExceeded(*fast,options.general.fitting.validation.distanceTolerance);
        }
        if(positionExceeded){++report.envelopePositionWitnesses;++report.positionOnlyNodes;}
        else if(tangentExceeded)++report.envelopeTangentWitnesses;
        else
        {
            if(useEnvelope)++report.envelopeFallbacks;
            bool fused=false;auto error=errors.Build(s,t,products,fused);
            report.fusedErrorBuilds+=fused;
            double bound=paired::PositionBoxUpper(error.numerator,error.denominator);
            if(bound>=options.general.fitting.validation.distanceTolerance)
                bound=std::min(bound,paired::PositionControlUpper(error.numerator,error.denominator,products));
            double angleRatio=INFINITY,angleSquared=INFINITY;
            if(!(bound<options.general.fitting.validation.distanceTolerance))++report.positionOnlyNodes;
            else if(trim.collapsed[edge])angleRatio=0;
            else
            {
                auto st=TangentNumerator(s,products);
                ++report.tangentFastChecks;
                angleRatio=ResidualTangentRatio(error,s,st,bound);
                if(angleRatio<limit)++report.tangentFastAccepted;
                else
                {
                    ++report.tangentFallbacks;
                    auto tt=TangentNumerator(t,products);
                    angleSquared=ScaledTangent(st,tt,tangentScale,false);
                }
            }
            if(Accept(bound,angleRatio,angleSquared))return true;
            // 部分边的位置已通过，切向仍需细分收紧。
            // 若两个切向求界内核均无结论，则关闭该边重复过滤以节省开销；其他边及位置超限见证路径仍保留过滤。
            // 这只影响成本，不改变验收条件。
            if(fastPositionPassed&&bound<options.general.fitting.validation.distanceTolerance)
                envelopeOnEdge[edge]=false;
        }
        if(MergeMicro(s,edge,lo,hi))return true;
        if(depth>=options.maxBoundaryDepth)return false;
        double mid=s.u0+(s.u1-s.u0)*.5,mapped=lo+(hi-lo)*.5;
        if(!(mid>s.u0&&mid<s.u1))return false;
        if(!trim.collapsed[edge])
        {
            const auto &src=samples.source;bool u=edge<2;
            auto p=samples.At(u?(edge==0?src.MinParamU():src.MaxParamU()):mid,
                u?mid:(edge==2?src.MinParamV():src.MaxParamV()),sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
            double candidate=Fraction(trim,edge,p,mapped);
            if(candidate>lo&&candidate<hi)mapped=candidate;
        }
        if(!(mapped>lo&&mapped<hi))return false;
        ++report.mappingRefinements;auto halves=Split(s,true);
        return Visit(halves.first,edge,lo,mapped,depth+1)&&Visit(halves.second,edge,mapped,hi,depth+1);
    }
public:
    const std::vector<GeneralBoundaryTangentExemption> &Exemptions() const {return exemptions;}
    Certifier(SurfaceSampleContext &s,const Trim &t,const Target &dst,const GeneralRationalBoundaryOptions &o,
        GeneralRationalBoundaryReport &r):samples(s),trim(t),target(dst),options(o),report(r)
    {
        I angle=I(o.general.fitting.validation.boundaryTangentToleranceDeg)*IPi()/I(180);
        limit=(Sin(angle)/Cos(angle)).lo;
        sineSquaredLimit=Sqr(Sin(angle));
    }
    bool Run(const std::array<std::vector<BezierCell>,4> &source)
    {
        for(unsigned e=0;e<4;++e)
        {
            if(source[e].empty())return false;
            double begin=source[e].front().u0,end=source[e].back().u1,previous=0;
            for(size_t i=0;i<source[e].size();++i)
            {
                const auto &cell=source[e][i];double next=i+1==source[e].size()?1:(cell.u1-begin)/(end-begin);
                if(i+1<source[e].size()&&!trim.collapsed[e])
                {
                    const auto &s=samples.source;bool u=e<2;
                    auto p=samples.At(u?(e==0?s.MinParamU():s.MaxParamU()):cell.u1,
                        u?cell.u1:(e==2?s.MinParamV():s.MaxParamV()),sggk::DerivType::NONE,sggk::DerivType::NONE,false).point;
                    next=Fraction(trim,e,p,next);
                }
                if(!(next>=previous&&next<=1)||!Visit(cell,e,previous,next,0))return false;
                previous=next;
            }
        }
        report.boundaryDistanceUpperBound=distance;
        report.boundaryTangentUpperBoundDeg=(Atan(I(ratio))*I(180)/IPi()).hi;
        if(sineSquared>0)
        {
            I residualAngle=Atan(Sqrt(I(sineSquared)/(I(1)-I(sineSquared))));
            report.boundaryTangentUpperBoundDeg=std::max(report.boundaryTangentUpperBoundDeg,
                                                       (residualAngle*I(180)/IPi()).hi);
        }
        return report.boundaryTangentUpperBoundDeg<options.general.fitting.validation.boundaryTangentToleranceDeg;
    }
};
std::vector<AnalyticParameters> Charts(SurfaceSampleContext &samples,const AnalyticParameters &a,
                                     const GeneralRationalBoundaryOptions &o)
{
    std::vector<AnalyticParameters> out{a};
    if(a.type!=sggk::SurfaceType::Sphere||o.chartDomain)return out;
    auto add=[&](V axis)
    {
        if(Norm(axis)<1e-12||out.size()>=o.maxChartAttempts)return;
        AnalyticParameters p=a;p.axis=Direction(axis);p.xDirection=Direction(Perpendicular(Unit(axis)));out.push_back(p);
    };
    // 仅生成有界数量的几何候选，不调用目标函数假设解析 UV 直线的旧优化器。
    for(const auto &edge:samples.Boundaries())
    {
        if(edge.size()<3)continue;
        add(Cross(edge[edge.size()/2].point-edge.front().point,edge.back().point-edge.front().point));
        if(Norm(edge.back().point-edge.front().point)<1e-12)add(edge.front().point-Get(a.origin));
    }
    add({1,0,0});add({0,1,0});add({0,0,1});return out;
}
}
std::vector<AnalyticParameters> RecoveryCharts(SurfaceSampleContext &samples,const AnalyticParameters &a,
                                              const GeneralRationalBoundaryOptions &o)
{return Charts(samples,a,o);}
// 当前合同的总编排：输入与环境检查 -> 内部采样 -> 候选参数图恢复 -> 有限域采样
// -> 实际 SDK 输出边界 -> 连续边界证书。内部两次采样都不能推出内部连续域已通过。
// 候选图失败可以换图，输入错误与预算耗尽由外层统一处理；成功时才提交 Target.patch。
GeneralRationalBoundaryReport Validate(SurfaceSampleContext &samples,const AnalyticParameters &a,
                                      const GeneralRationalBoundaryOptions &o)
{
    GeneralRationalBoundaryReport r;auto start=Clock::now();
    try
    {
        ExecutionDeadline localDeadline(o.maxTotalMs,start);
        DeadlineScope deadlineScope(samples.deadline,localDeadline);
        const auto &s=samples.source;ValidateSurface(s);const auto &v=o.general.fitting.validation;
        if(!(v.distanceTolerance>0&&std::isfinite(v.distanceTolerance))||
           !(v.boundaryTangentToleranceDeg>0&&v.boundaryTangentToleranceDeg<90)||
           !(o.microSegmentDiameter>=0&&o.microSegmentDiameter<=1e-8)||
           !o.maxBoundaryNodes||o.maxBoundaryDepth>20||!o.maxChartAttempts||o.maxChartAttempts>64||
           o.general.boundarySamplesPerSpan<3||o.general.boundarySamplesPerSpan>64||!o.general.fitting.sampling.maxSamples)
            throw Failure(Status::InvalidInput,"Invalid rational boundary options");
        if(o.chartDomain)for(unsigned d=0;d<2;++d)
            if(!(std::isfinite((*o.chartDomain)[2*d])&&std::isfinite((*o.chartDomain)[2*d+1])&&
                 (*o.chartDomain)[2*d]<(*o.chartDomain)[2*d+1]))throw Failure(Status::InvalidInput,"Invalid fixed chart box");
        if(!IntervalEnvironmentSupported())throw Failure(Status::Unsupported,"Directed interval environment unavailable");
        using T=sggk::SurfaceType;
        if(a.type!=T::Plane&&a.type!=T::Cylinder&&a.type!=T::Cone&&a.type!=T::Sphere&&a.type!=T::Torus)
            throw Failure(Status::Unsupported,"Unsupported rational mother");
        if(a.type!=T::Plane&&a.type!=T::Cone&&!(a.radius>0&&std::isfinite(a.radius)))
            throw Failure(Status::InvalidInput,"Invalid mother radius");
        if(a.type==T::Torus&&!(a.majorRadius>a.radius&&std::isfinite(a.majorRadius)))
            throw Failure(Status::Unsupported,"Only regular ring torus supported");
        if(a.type==T::Cone&&!(a.halfAngle>0&&a.halfAngle<Pi/2))throw Failure(Status::InvalidInput,"Invalid cone angle");
        auto frame=FrameOf(a);if(!Finite(frame.o)||!Finite(frame.x)||!Finite(frame.z))throw Failure(Status::InvalidInput,"Invalid mother frame");
        samples.CheckDeadline("validation_input");
        if(s.IsRational())for(const auto &row:s.Weights())for(double w:row)
            if(!(w>0))throw Failure(Status::Unsupported,"Positive source weights required");
        for(bool u:{false,true})
        {
            const auto &m=u?s.MultsU():s.MultsV();unsigned degree=u?s.DegreeU():s.DegreeV();
            for(size_t i=1;i+1<m.size();++i)if(m[i]>degree)throw Failure(Status::Unsupported,"Discontinuous source knots");
        }
        {
            Measure timer(r.samplingMs);
            for(const auto &p:samples.Greville(false))
            {++r.samples;r.interiorSampleMaxDistance=std::max(r.interiorSampleMaxDistance,Norm(p.point-ProjectToMother(a,p.point)));}
        }
        if(!r.samples||!(r.interiorSampleMaxDistance<v.distanceTolerance))
            throw Failure(Status::VerificationFailed,"Interior sampled mother distance exceeded");
        std::vector<AnalyticParameters> charts;
        {Measure timer(r.recoveryMs);charts=Charts(samples,a,o);}
        std::optional<std::array<std::vector<BezierCell>,4>> source;
        for(const auto &chart:charts)
        {
            samples.CheckDeadline("chart_recovery");
            ++r.chartAttempts;Trim trim;bool recovered=false;
            try{Measure timer(r.recoveryMs);recovered=Recover(samples,chart,o,trim);}
            catch(const Failure &e){r.message=e.what();continue;}
            r.recovery3DRefinements+=trim.refinementCalls;
            if(!recovered){r.message="NURBS-domain straight-boundary recovery unresolved";continue;}
            {
                Measure timer(r.samplingMs);double bound=0;
                for(const auto &p:samples.Greville(false))
                {samples.CheckDeadline("finite_patch_sampling");bound=std::max(bound,FiniteDistance(trim,p.point));}
                r.interiorSampleMaxFiniteDistance=bound;
                if(!(bound<v.distanceTolerance)){r.message="Finite-patch sample witness exceeded";continue;}
            }
            Target target;
            {
                Measure timer(r.preparationMs);
                if(!source)source=ExtractHybridBoundaries(s,Get(a.origin));
                samples.CheckDeadline("source_preparation");
                target=BuildTarget(trim,v.distanceTolerance*.1,r,samples.deadline);
                target.patch.sourceDomain={s.MinParamU(),s.MaxParamU(),s.MinParamV(),s.MaxParamV()};
            }
            bool passed;
            {
                Measure timer(r.boundaryMs);Certifier validator(samples,trim,target,o,r);passed=validator.Run(*source);
                if(passed)target.patch.tangentExemptions=validator.Exemptions();
            }
            samples.CheckDeadline("validation_complete");
            if(passed)
            {
                r.outputMotherDistanceUpperBound=target.motherBound;r.patch=std::move(target.patch);
                r.status=HybridValidationStatus::Passed;
                r.message=r.patch->tangentExemptions.empty()?
                    "Sampled interior; continuous exported NURBS boundary position/tangent; no analytic PCurve":
                    "Sampled interior; continuous boundary position; tangent excludes documented spatial micro-segments";
                break;
            }
            r.message="Rational boundary correspondence not certified";
        }
    }
    catch(const DeadlineExceeded &e)
    {r.patch.reset();r.status=HybridValidationStatus::BudgetExceeded;r.timeBudgetExceeded=true;r.timeoutStage=e.stage;r.message=e.what();}
    catch(const std::invalid_argument &e){r.status=HybridValidationStatus::InvalidInput;r.message=e.what();}
    catch(const Budget &e){r.status=HybridValidationStatus::BudgetExceeded;r.message=e.what();}
    catch(const Failure &e)
    {r.status=e.status==Status::InvalidInput?HybridValidationStatus::InvalidInput:e.status==Status::Unsupported?
        HybridValidationStatus::Unsupported:HybridValidationStatus::NotCertified;r.message=e.what();}
    catch(const std::exception &e){r.status=HybridValidationStatus::NotCertified;r.message=e.what();}
    if(!r.Passed())
    {r.patch.reset();r.boundaryDistanceUpperBound=r.boundaryTangentUpperBoundDeg=r.outputMotherDistanceUpperBound=INFINITY;}
    r.totalMs=Ms(start);
    if(o.maxTotalMs>0&&r.totalMs>o.maxTotalMs)
    {
        r.patch.reset();r.status=HybridValidationStatus::BudgetExceeded;r.timeBudgetExceeded=true;
        r.timeoutStage="validation_completion_overrun";r.message="Validation wall time exceeded hard limit";
        r.boundaryDistanceUpperBound=r.boundaryTangentUpperBoundDeg=r.outputMotherDistanceUpperBound=INFINITY;
    }
    return r;
}
}
namespace bspline2analytic
{
GeneralRationalBoundaryReport ValidateGeneralRationalBoundaries(const sggk::BSplineSurface &s,
    const AnalyticParameters &a,const GeneralRationalBoundaryOptions &o)
{
    auto start=std::chrono::steady_clock::now();detail::SurfaceSampleContext samples(s,o.general.fitting.sampling.maxSamples);
    auto r=detail::rational_boundary::Validate(samples,a,o);
    r.totalMs=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();return r;
}
}
