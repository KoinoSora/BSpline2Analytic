#pragma once
#include "AnalyticReconstruction.h"
#include "SurfaceSampleContext.h"
namespace bspline2analytic::detail
{
ValidationReport ValidateCandidate(const sggk::BSplineSurface &s, const Candidate &c,
                                   const std::array<BoundaryInfo, 4> &b,
                                   const RecognitionOptions &o, SurfaceSampleContext * = nullptr);
}
