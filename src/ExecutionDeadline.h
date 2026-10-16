#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <stdexcept>

namespace bspline2analytic::detail
{
// Deliberately not a runtime_error: geometric candidate catches must not swallow cancellation.
struct DeadlineExceeded : std::exception
{
    const char *stage;
    explicit DeadlineExceeded(const char *where) : stage(where) {}
    const char *what() const noexcept override { return "Case time budget exhausted"; }
};
class ExecutionDeadline
{
public:
    using Clock = std::chrono::steady_clock;
    explicit ExecutionDeadline(double milliseconds, Clock::time_point start = Clock::now())
        : enabled(milliseconds > 0)
    {
        if (!std::isfinite(milliseconds) || milliseconds < 0 || milliseconds > 60000)
            throw std::invalid_argument("Invalid case time budget");
        // Cooperative cancellation: reserve 10% (at most 5ms) for the current
        // indivisible operation and stack cleanup. This is not a real-time OS guarantee.
        if (enabled) stop = start + std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double, std::milli>(milliseconds - std::min(5., milliseconds * .1)));
    }
    void Check(const char *stage) const
    { if (enabled && Clock::now() >= stop) throw DeadlineExceeded(stage); }
private:
    bool enabled;
    Clock::time_point stop{};
};
class DeadlineScope
{
    const ExecutionDeadline *&slot;
    const ExecutionDeadline *previous;
public:
    DeadlineScope(const ExecutionDeadline *&slot, const ExecutionDeadline &local)
        : slot(slot), previous(slot) { if (!slot) slot = &local; }
    ~DeadlineScope() { slot = previous; }
};
}
