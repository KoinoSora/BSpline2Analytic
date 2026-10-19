/** @file AnalyticExclusionInternal.h
 * @brief 从源曲面实际区间点值构造排除证书，不接收候选或测试真值。
 */
#pragma once
#include "AnalyticExclusion.h"
#include "ExecutionDeadline.h"
#include "IntervalBounds.h"
#include <Geometry/3D/Surface/BSplineSurface.h>
namespace bspline2analytic::detail
{
// 两方向次数须 <=5。输出顺序为 U 外层、V 内层；固定四分点处用精确二进制基函数求区间包络。
std::array<I,25> ExclusionQuarterGrid(const BPoly &);
void CertifyAnalyticExclusion(const sggk::BSplineSurface &, double distanceTolerance,
                             const ExecutionDeadline &, AnalyticExclusionReport &,
                             double localBudgetMs = 8);
}
