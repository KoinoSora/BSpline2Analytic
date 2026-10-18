/**
 * @file PairedSurfaceBounds.h
 * @brief 测试来源面与输入面的同参数距离认证接口。
 * 覆盖完整公共域，只用于输入资格和独立参考，不参与生产候选生成。
 */
#pragma once
#include "BSpline2Analytic.h"
namespace bspline2analytic::tests
{
struct PairBound
{
    bool passed = false;
    double upper = std::numeric_limits<double>::infinity();
    unsigned pieces = 0;
};
// 完整公共参数域上的同参数距离界，与识别算法独立。
PairBound BoundPair(const sggk::BSplineSurface &a, const sggk::BSplineSurface &b, double tolerance,
                    unsigned maxPieces = 50000);
} // namespace bspline2analytic::tests
