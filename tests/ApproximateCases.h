#pragma once
#include "ApproximateConversion.h"
#include "ReferenceCases.h"
namespace bspline2analytic::tests
{
sggk::BSplineSurfacePtr RefitPolynomial(const sggk::BSplineSurface &, unsigned subdivisions);
sggk::BSplineSurfacePtr PlateauFixture(const sggk::BSplineSurface &);
} // namespace bspline2analytic::tests
