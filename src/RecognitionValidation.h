/**
 * @file RecognitionValidation.h
 * @brief 精确识别候选验证的内部接口。
 * 验证等级与误差上界随报告返回，调用方不得把采样级结果当作连续证书。
 */
#pragma once
#include "AnalyticReconstruction.h"
#include "SurfaceSampleContext.h"
namespace bspline2analytic::detail
{
ValidationReport ValidateCandidate(const sggk::BSplineSurface &s, const Candidate &c,
                                   const std::array<BoundaryInfo, 4> &b,
                                   const RecognitionOptions &o, SurfaceSampleContext * = nullptr);
}
