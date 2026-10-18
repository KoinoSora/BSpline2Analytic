/**
 * @file StageTimer.h
 * @brief 转换层互斥阶段的 RAII 计时工具。
 * 用于分类、优化和额外验证阶段，不放入逐点循环，不重复累加嵌套细分明细。
 */
#pragma once
#include "ApproximateConversion.h"
#include <chrono>
namespace bspline2analytic::detail
{
// 只在互斥的流程阶段边界使用，避免在逐点循环内引入计时开销。
class StageTimer
{
    using Clock = std::chrono::steady_clock;
    StageTiming &output;
    StageTiming *aggregate;
    Clock::time_point start = Clock::now();

  public:
    explicit StageTimer(StageTiming &t, StageTiming *phase = nullptr) : output(t), aggregate(phase)
    {
        ++output.calls;
        if (aggregate) ++aggregate->calls;
    }
    StageTimer(const StageTimer &) = delete;
    StageTimer &operator=(const StageTimer &) = delete;
    ~StageTimer() noexcept
    {
        const double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
        output.elapsedMs += elapsed;
        if (aggregate) aggregate->elapsedMs += elapsed;
    }
};
} // namespace bspline2analytic::detail
