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
// Same-parameter distance on the complete common domain; independent of recognition.
PairBound BoundPair(const sggk::BSplineSurface &a, const sggk::BSplineSurface &b, double tolerance,
                    unsigned maxPieces = 50000);
} // namespace bspline2analytic::tests
