#pragma once
#include "GeneralHybridValidation.h"
#include "GeneralInternal.h"
namespace bspline2analytic::detail
{
// Cache is owned by the caller and refers to an unchanged source. Samples are
// reused, but every candidate's distances and certificates are recomputed.
GeneralHybridReport ValidateGeneralHybridWithSamples(const sggk::BSplineSurface &,
    const AnalyticParameters &, const GeneralHybridOptions &, SurfaceSampleContext &);
// One-dimensional homogeneous extraction; also exposed to white-box tests.
std::array<std::vector<BezierCell>,4> ExtractHybridBoundaries(const sggk::BSplineSurface &, V origin);
}
