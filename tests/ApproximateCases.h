/**
 * @file ApproximateCases.h
 * @brief 近似测试几何与期望信息的声明。
 * 仅为测试使用；共享生成器可复用，是否运行套件须遵守其验证器状态。
 */
#pragma once
#include "ApproximateConversion.h"
#include "ReferenceCases.h"
namespace bspline2analytic::tests
{
sggk::BSplineSurfacePtr RefitPolynomial(const sggk::BSplineSurface &, unsigned subdivisions);
struct PerturbationCase
{
    std::string name;
    sggk::BSplineSurfacePtr spline;
    AnalyticPatch initial, feasible;
};
PerturbationCase CurvedPlateauFixture(const ReferenceCase &);
sggk::BSplineSurfacePtr PlateauFixture(const sggk::BSplineSurface &);
} // namespace bspline2analytic::tests
