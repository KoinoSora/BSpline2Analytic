/**
 * @file ValidationTiming.h
 * @brief 验证器五阶段互斥计时与统计归并。
 * 嵌套阶段切换时扣除外层重复时间，返回报告前显式结束计时；输出日志由调用方另计。
 */
#pragma once
#include "ApproximateConversion.h"
#include <chrono>
#include <utility>

namespace bspline2analytic::detail
{
inline void AddValidationTiming(ValidationTiming &out, const ValidationTiming &in)
{
    for (auto pair : {std::make_pair(&out.preparation, &in.preparation),
                      std::make_pair(&out.surfaceBounds, &in.surfaceBounds),
                      std::make_pair(&out.subdivision, &in.subdivision),
                      std::make_pair(&out.boundaryChecks, &in.boundaryChecks),
                      std::make_pair(&out.reporting, &in.reporting)})
    {
        pair.first->elapsedMs += pair.second->elapsedMs;
        pair.first->calls += pair.second->calls;
    }
    out.totalMs += in.totalMs;
}

// 相邻两次读钟之间的时间只归入一个活动阶段；嵌套 Scope 退出后恢复调用方阶段，避免将嵌套工作重复计入外层。
class ValidationTimer
{
    using Clock = std::chrono::steady_clock;
    ValidationTiming &output;
    StageTiming *active;
    Clock::time_point start = Clock::now(), last = start;

    void Charge(Clock::time_point now)
    {
        active->elapsedMs += std::chrono::duration<double, std::milli>(now - last).count();
        last = now;
    }

  public:
    explicit ValidationTimer(ValidationTiming &out) : output(out), active(&out.preparation)
    {
        ++active->calls;
    }
    ValidationTimer(const ValidationTimer &) = delete;
    ValidationTimer &operator=(const ValidationTimer &) = delete;
    void Switch(StageTiming &next, bool entering = true)
    {
        Charge(Clock::now());
        active = &next;
        if (entering) ++active->calls;
    }
    class Scope
    {
        ValidationTimer &timer;
        StageTiming *previous;
      public:
        Scope(ValidationTimer &owner, StageTiming &next) : timer(owner), previous(owner.active)
        {
            timer.Switch(next);
        }
        Scope(const Scope &) = delete;
        Scope &operator=(const Scope &) = delete;
        ~Scope() { timer.Switch(*previous, false); }
    };
    // 返回或移动报告前显式结束计时，不依赖具名返回值优化 NRVO。
    void Finish()
    {
        const auto now = Clock::now();
        Charge(now);
        output.totalMs = std::chrono::duration<double, std::milli>(now - start).count();
    }
};
} // namespace bspline2analytic::detail
