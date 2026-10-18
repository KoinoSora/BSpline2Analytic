/**
 * @file AnalyticReconstruction.h
 * @brief 精确识别内部的解析母面重建声明。
 * 承接边界分类，不对重建候选直接作最终验收。
 */
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
