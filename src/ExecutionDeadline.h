/**
 * @file ExecutionDeadline.h
 * @brief 转换调用共享的单调时钟截止时间与作用域传递。
 * 使用协作式检查传播取消；内层不覆盖外层预算，返回时由入口复核真实总耗时。
 */
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <stdexcept>

namespace bspline2analytic::detail
{
// 特意不继承 runtime_error，防止候选几何失败的捕获分支吞掉超时取消信号。
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
        // 协作式取消：提前预留预算的 10%（至多 5 ms），供不可分割操作完成和栈清理使用。
        // 这不构成实时操作系统级别的硬中断保证，调用返回时仍须复核真实耗时。
        if (enabled) stop = start + std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double, std::milli>(milliseconds - std::min(5., milliseconds * .1)));
    }
    void Check(const char *stage) const
    { if (enabled && Clock::now() >= stop) throw DeadlineExceeded(stage); }
private:
    bool enabled;
    Clock::time_point stop{};
};
// 将局部预算临时挂入共享上下文。已有外层预算优先，不能被每候选局部预算延长。
// 析构恢复原指针，保证异常退出后上下文不会悬挂到已经销毁的局部截止时间。
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
