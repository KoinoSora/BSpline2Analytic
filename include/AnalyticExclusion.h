/** @file AnalyticExclusion.h
 * @brief 五类母面的整体距离排除报告；未证明不表示可以转换。
 */
#pragma once
#include <string>
namespace bspline2analytic
{
enum class AnalyticExclusionStatus { NotAttempted, Inconclusive, CertifiedExcluded, LocalBudgetExceeded };
struct AnalyticExclusionReport
{
    AnalyticExclusionStatus status = AnalyticExclusionStatus::NotAttempted;
    unsigned attempts = 0, witnesses = 0;
    double totalMs = 0, distanceLowerBound = 0;
    double residualNormUpper = 0, inverseNormUpper = 0, gradientBound = 0, scale = 0;
    std::string message;
    bool Excluded() const noexcept { return status == AnalyticExclusionStatus::CertifiedExcluded; }
};
}
