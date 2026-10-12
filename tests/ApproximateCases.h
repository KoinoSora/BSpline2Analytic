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
