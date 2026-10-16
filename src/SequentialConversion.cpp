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
            candidate.domain = {0,1,.25,1};
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
        detail::ExecutionDeadline deadline(o.maxTotalMs,start);
        detail::ExecutionDeadline isoDeadline(std::min(o.maxTotalMs,o.maxIsoMs),start);
        try { Iso(s,o,isoDeadline,r); }
        catch (const detail::DeadlineExceeded &e)
        {
            r.accepted=false; r.isoTimeBudgetExceeded=true;
            r.isoMessage=std::string("Iso sub-budget exhausted: ")+e.stage;
            r.isoPatch.reset(); r.exact.surface.reset(); r.exact.domain.reset();
        }
        if (!r.accepted) Fallback(s,o,deadline,r);
        deadline.Check("sequential_complete");
    }
    catch (const detail::DeadlineExceeded &e)
    { r.accepted=false; r.timeBudgetExceeded=true; r.timeoutStage=e.stage; r.message=e.what(); }
    catch (const std::exception &e) { r.accepted=false; r.message=e.what(); }
    r.totalMs = Elapsed(start);
    if (r.totalMs > o.maxTotalMs && o.maxTotalMs > 0)
    { r.accepted=false; r.timeBudgetExceeded=true; r.timeoutStage="completion_overrun"; r.message="Actual wall time exceeded case limit"; }
    if (!r.Passed()) { r.isoPatch.reset(); r.exact.surface.reset(); r.exact.domain.reset(); r.fallback.patch.reset(); }
    return r;
}
}
