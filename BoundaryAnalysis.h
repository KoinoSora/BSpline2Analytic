#pragma once
#include "BSpline2Analytic.h"
namespace bspline2analytic::detail
{
BoundaryInfo AnalyzeCurve(const sggk::Curve3DPtr &curve, BoundaryId id,
                          const RecognitionOptions &o);
}
