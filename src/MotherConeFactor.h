/**
 * @file MotherConeFactor.h
 * @brief 圆锥位置后备：P=gQ+R 的非负线性因子及区间距离界。
 * 分解属于当前候选；误差作为距离扰动补回，不除以可能为零的 g。
 */
#pragma once
#include "MotherDistanceInternal.h"

namespace bspline2analytic::detail::mother
{
struct ConeFactorCell
{
    BPoly residual, radial, z, factor;
    BV remainder;
    I cosine, sine;
    double enclosure = 0;
};
struct ConeFactorEstimate
{
    double upper = INFINITY, perturbation = INFINITY, denominator = 0;
};
// 所有长度均使用当前候选的归一化局部坐标。
std::optional<ConeFactorCell> TryConeFactor(const std::array<Poly, 3> &, double enclosure,
    double angle, double tolerance, paired::Products &, MotherDistanceValidationReport &,
    const ExecutionDeadline *);
ConeFactorEstimate BoundConeFactor(const ConeFactorCell &);
std::pair<ConeFactorCell, ConeFactorCell> SplitConeFactor(const ConeFactorCell &, bool alongU);
} // namespace bspline2analytic::detail::mother
