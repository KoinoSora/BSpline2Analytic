/**
 * @file ReferenceCases.h
 * @brief 精确识别参考案例的数据结构与生成声明。
 * 案例几何和预期类型供外部审计使用，不注入生产识别路径。
 */
#pragma once
#include "BSpline2Analytic.h"
#include <vector>
namespace bspline2analytic::tests
{
struct ReferenceCase
{
    std::string name;
    sggk::SurfacePtr analytic;
    sggk::BSplineSurfacePtr spline;
};
std::vector<ReferenceCase> MakeCases();
sggk::BSplineSurfacePtr ElevateBezierPatch(const sggk::BSplineSurface &surface);
} // namespace bspline2analytic::tests
