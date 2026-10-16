#pragma once
#include "MotherDistanceValidation.h"
#include "PairedValidationInternal.h"

namespace bspline2analytic::detail::mother
{
struct TangentSegment
{
    double begin = 0, end = 1;
    // Source U/V directions, including target chart swaps/reversals. Packed in
    // one Bernstein variable; positive parameter/weight scales may be omitted.
    std::array<BV, 2> directions;
};
using BoundaryTangents = std::array<std::vector<TangentSegment>, 4>;
// Extract only the outer row and its adjacent row, never a full 2D derivative.
std::array<BV, 2> BoundaryPartials(const BV &, BoundaryId);
BoundaryTangents PrepareTargetTangents(const paired::Chart &, paired::Products &);
// Certify with A-lambda*B, lambda>0, and a source-length lower bound.
// No spline cross products or normals; inconclusive bounds share one split budget.
bool ValidateJointTangents(const BoundaryTangents &, const BoundaryTangents &,
                          const MotherDistanceValidationOptions &, PairedValidationReport &,
                          paired::Products &, MotherDistanceValidationReport &);
} // namespace bspline2analytic::detail::mother
