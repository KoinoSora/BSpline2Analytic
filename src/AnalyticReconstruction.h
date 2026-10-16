#pragma once
#include "BSpline2Analytic.h"
#include <vector>
namespace bspline2analytic::detail
{
struct Candidate
{
    AnalyticParameters parameters;
    sggk::SurfacePtr surface;
};
std::vector<Candidate> Reconstruct(const std::array<BoundaryInfo, 4> &b,
                                   const RecognitionOptions &o, bool materialize = true);
double Distance(const AnalyticParameters &p, const sggk::Point3D &point);
} // namespace bspline2analytic::detail
