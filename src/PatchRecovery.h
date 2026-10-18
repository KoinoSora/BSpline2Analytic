/**
 * @file PatchRecovery.h
 * @brief 精确识别内部的有限面片恢复接口。
 * 承接解析参数并生成有限域与方向信息，后续仍须验证。
 */
#pragma once
#include "AnalyticReconstruction.h"
namespace bspline2analytic::detail
{
void RecoverPatch(const sggk::BSplineSurface &input, const Candidate &candidate,
                  const std::array<BoundaryInfo, 4> &boundaries, const RecognitionOptions &options,
                  sggk::UVRange &range, ParameterCorrespondence &mapping);
}
