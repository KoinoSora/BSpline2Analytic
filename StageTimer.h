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
    Clock::time_point start = Clock::now();

  public:
    explicit StageTimer(StageTiming &t) : output(t)
    {
        ++output.calls;
    }
    StageTimer(const StageTimer &) = delete;
    StageTimer &operator=(const StageTimer &) = delete;
    ~StageTimer() noexcept
    {
        output.elapsedMs += std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    }
};
} // namespace bspline2analytic::detail
