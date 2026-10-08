#pragma once
#include "AnalyticReconstruction.h"
namespace bspline2analytic::detail
{
void RecoverPatch(const sggk::BSplineSurface &input, const Candidate &candidate,
                  const std::array<BoundaryInfo, 4> &boundaries, const RecognitionOptions &options,
                  sggk::UVRange &range, ParameterCorrespondence &mapping);
}
