/**
 * @file MotherTangentValidation.h
 * @brief 边界两参数方向及联合切向认证的内部接口。
 * 只需要边界行与相邻控制行，不构造源面法向；联合条件不是直接法向测量。
 */
#pragma once
#include "MotherDistanceValidation.h"
#include "PairedValidationInternal.h"

namespace bspline2analytic::detail::mother
{
struct TangentSegment
{
    double begin = 0, end = 1;
    // 按源 U/V 对齐的方向，包括目标图的交换与反向。
    // 压入单个 Bernstein 变量；可省略不影响方向的正参数尺度和正权重尺度。
    std::array<BV, 2> directions;
};
using BoundaryTangents = std::array<std::vector<TangentSegment>, 4>;
// 只提取最外控制行及相邻行，不构造完整二维导数场。
std::array<BV, 2> BoundaryPartials(const BV &, BoundaryId);
BoundaryTangents PrepareTargetTangents(const paired::Chart &, paired::Products &);
// 通过 A-lambda*B、lambda>0 及源切向长度下界认证；不构造样条叉乘或法向。
// 无法判定时，两方向共享一次细分预算。
bool ValidateJointTangents(const BoundaryTangents &, const BoundaryTangents &,
                          const MotherDistanceValidationOptions &, PairedValidationReport &,
                          paired::Products &, MotherDistanceValidationReport &, const ExecutionDeadline * = nullptr);
} // namespace bspline2analytic::detail::mother
