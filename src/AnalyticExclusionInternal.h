/** @file AnalyticExclusionInternal.h
 * @brief 从源曲面实际区间点值构造排除证书，不接收候选或测试真值。
 */
#pragma once
#include "AnalyticExclusion.h"
#include "ExecutionDeadline.h"
#include <Geometry/3D/Surface/BSplineSurface.h>
namespace bspline2analytic::detail
{
void CertifyAnalyticExclusion(const sggk::BSplineSurface &, double distanceTolerance,
                             const ExecutionDeadline &, AnalyticExclusionReport &,
                             double localBudgetMs = 8);
}
