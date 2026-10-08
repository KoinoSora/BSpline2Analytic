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
