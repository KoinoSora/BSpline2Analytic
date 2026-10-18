/**
 * @file SequentialConversion.cpp
 * @brief 当前串联转换入口及两条分支的调度实现。
 * 等参分支先行，失败才进入一般候选与母面优化；总计时包含失败尝试，超时清除全部可发布几何。
 */
#include "SequentialConversion.h"
#include "MotherDistanceInternal.h"
#include "GeneralRationalBoundaryInternal.h"

namespace bspline2analytic
{
namespace
{
using namespace detail;
double Elapsed(ExecutionDeadline::Clock::time_point t)
{ return std::chrono::duration<double,std::milli>(ExecutionDeadline::Clock::now()-t).count(); }
struct Timer
{
    double &ms;
    ExecutionDeadline::Clock::time_point start = ExecutionDeadline::Clock::now();
    ~Timer() { ms = Elapsed(start); }
};
// 等参分支只接收源几何和配置。有理面保留精确识别合同；多项式面逐个恢复候选，
// 交给母面距离与联合切向验证。一次候选失败只推进搜索，不直接判定源面不可转换。
// deadline 为等参子预算；超时必须向上传递，交由入口决定是否仍有兜底时间。
void Iso(const sggk::BSplineSurface &s, const SequentialConversionOptions &o,
         const ExecutionDeadline &deadline, SequentialConversionResult &r)
{
    Timer timer{r.isoMs};
    r.isoAttempted = true;
    r.path = s.IsRational() ? "iso_exact" : "iso_mother";
    try
    {
        deadline.Check("iso_start");
        if (s.IsRational())
        {
            SurfaceSampleContext samples(s);
            samples.deadline = &deadline;
            r.exact = RecognizeWithSamples(s,o.exact,samples);
            r.accepted = r.exact.status == Status::Success;
            r.message = r.exact.message;
        }
        else
        {
            paired::AnalysisContext context(s);
            context.samples.deadline = &deadline;
            context.classification.emplace();
            ClassificationDiagnostics diagnostics;
            IsoClassificationSearch search(context.samples,o.iso,diagnostics);
            while (auto candidate = search.Next(*context.classification))
            {
                deadline.Check("iso_candidate");
                ++r.isoCandidates;
                r.mother = mother::ValidatePrepared(context,*candidate,o.mother);
                r.message = r.mother.message;
                if (r.mother.Passed())
                { r.isoPatch = candidate->patch; r.accepted = true; break; }
            }
            if (r.message.empty()) r.message = "Iso search exhausted without a candidate";
        }
        deadline.Check("iso_complete");
    }
    catch (const DeadlineExceeded &) { throw; }
    catch (const std::exception &e) { r.accepted = false; r.message = e.what(); }
    r.isoMessage = r.message;
}
// 一般兜底复用源采样，但每个母面候选必须重做有理参数域恢复和连续边界认证。
// 先验证初值；失败且允许优化时，只优化母面自由度，再用相同合同独立复验。
// 候选池、残差评价和 Jacobian 全部沿用入口的总截止时间，不能按候选重新计满预算。
void Fallback(const sggk::BSplineSurface &s, const SequentialConversionOptions &o,
              const ExecutionDeadline &deadline, SequentialConversionResult &r)
{
    Timer timer{r.fallbackMs};
    r.fallbackAttempted = true;
    r.path = "general_rational";
    SurfaceSampleContext samples(s,o.fallback.general.fitting.sampling.maxSamples);
    samples.deadline = &deadline;
    deadline.Check("fallback_start");
    GeneralClassificationDiagnostics diagnostics;
    auto pool = InitializeGeneralCandidates(samples,o.fallback.general,diagnostics);
    ConversionStatistics statistics;
    std::optional<SampleCache> base;
    for (;;)
    {
        auto batch = NextGeneralCandidateBatch(pool,samples,o.fallback.general,diagnostics);
        if (batch.empty()) break;
        for (const auto &seed : batch)
        {
            deadline.Check("fallback_candidate");
            ++r.fallbackCandidates;
            AnalyticPatch candidate;
            candidate.parameters = seed.parameters;
            // 这是母面参数更新所需的有效临时域，不是恢复出的裁剪域，也不作为验收真值。
            candidate.domain = {0,1,.25,1};
            // 捕获 candidate 的当前参数；优化原地更新后再次调用，必须验证新的母面。
            auto validate = [&] {
                r.fallback = rational_boundary::Validate(samples,candidate.parameters,o.fallback);
                if (r.fallback.timeBudgetExceeded) throw DeadlineExceeded("fallback_validation");
            };
            validate();
            const auto &fit = o.fallback.general.fitting;
            if (!r.fallback.Passed() && fit.refinement.maxRounds && fit.refinement.maxIterations)
            {
                deadline.Check("fallback_optimization");
                if (!base) base = GeneralSamples(samples);
                EvaluationWorkspace workspace;
                auto problem = MotherRefinementProblem(*base,fit);
                auto evaluate = problem.evaluate;
                problem.evaluate = [&](const AnalyticPatch &p, EvaluationWorkspace &w, ResidualBuffer &b) {
                    deadline.Check("fallback_optimization"); evaluate(p,w,b);
                };
                auto jacobian = problem.jacobian;
                problem.jacobian = [&](const AnalyticPatch &p, const ActiveParameters &a, double scale,
                                       EvaluationWorkspace &w, RefinementDiagnostics &d) {
                    deadline.Check("fallback_jacobian"); return jacobian(p,a,scale,w,d);
                };
                auto refined = RefineCandidateLevel(candidate,*base,fit,statistics,
                    ParameterLevel::MotherSurface,fit.refinement.maxIterations,true,workspace,&problem);
                r.optimizationIterations = statistics.optimizationIterations;
                if (refined.improved) validate();
            }
            r.message = r.fallback.message;
            if (r.fallback.Passed()) { r.accepted = true; return; }
        }
    }
    if (r.message.empty()) r.message = "No fallback candidate certified";
}
}
SequentialConversionResult ConvertIsoThenFallback(const sggk::BSplineSurface &s,
                                                  const SequentialConversionOptions &o)
{
    const auto start = detail::ExecutionDeadline::Clock::now();
    SequentialConversionResult r;
    try
    {
        if (!(o.maxTotalMs > 0 && o.maxTotalMs <= 50 && o.maxIsoMs > 0 && std::isfinite(o.maxIsoMs)))
            throw std::invalid_argument("Sequential total budget must be positive and <=50 ms; iso budget must be positive and finite");
        // 两个截止时间共用 start。等参分支花掉的时间直接占用总预算，兜底不能重置它。
        detail::ExecutionDeadline deadline(o.maxTotalMs,start);
        detail::ExecutionDeadline isoDeadline(std::min(o.maxTotalMs,o.maxIsoMs),start);
        try { Iso(s,o,isoDeadline,r); }
        catch (const detail::DeadlineExceeded &e)
        {
            r.accepted=false; r.isoTimeBudgetExceeded=true;
            r.isoMessage=std::string("Iso sub-budget exhausted: ")+e.stage;
            // 子预算耗尽即撤销等参候选，即使临近超时曾有局部检查通过，也不能发布。
            r.isoPatch.reset(); r.exact.surface.reset(); r.exact.domain.reset();
        }
        // 成功直接短路；失败才进入兜底。若总预算已耗尽，兜底的首个检查立即取消。
        if (!r.accepted) Fallback(s,o,deadline,r);
        deadline.Check("sequential_complete");
    }
    catch (const detail::DeadlineExceeded &e)
    { r.accepted=false; r.timeBudgetExceeded=true; r.timeoutStage=e.stage; r.message=e.what(); }
    catch (const std::exception &e) { r.accepted=false; r.message=e.what(); }
    // 协作式检查无法中断正在执行的不可分割操作，必须用实际完成时间再次兜住超时。
    // 不截断为 50 ms：报告真实墙钟，便于调用方区分超时与几何未认证。
    r.totalMs = Elapsed(start);
    if (r.totalMs > o.maxTotalMs && o.maxTotalMs > 0)
    { r.accepted=false; r.timeBudgetExceeded=true; r.timeoutStage="completion_overrun"; r.message="Actual wall time exceeded case limit"; }
    // 失败结果保留诊断，但移除所有分支的可用几何，防止调用方误用失败尝试的候选。
    if (!r.Passed()) { r.isoPatch.reset(); r.exact.surface.reset(); r.exact.domain.reset(); r.fallback.patch.reset(); }
    return r;
}
}
