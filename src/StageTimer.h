#pragma once
#include "ApproximateConversion.h"
#include <chrono>
namespace bspline2analytic::detail
{
// Use only at disjoint orchestration boundaries, never in point loops.
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
