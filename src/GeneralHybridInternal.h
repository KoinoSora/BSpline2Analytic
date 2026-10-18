/**
 * @file GeneralHybridInternal.h
 * @brief 历史混合验证器共享的源边提取与内部声明。
 * 当前有理边界路径复用提取工具，不因此调用已废弃 ValidateGeneralHybrid。
 */
#pragma once
#include "GeneralHybridValidation.h"
#include "GeneralInternal.h"
namespace bspline2analytic::detail
{
// 缓存由调用方持有，要求源面不变。
// 可复用采样值，但每个候选的距离与证书必须重新计算。
GeneralHybridReport ValidateGeneralHybridWithSamples(const sggk::BSplineSurface &,
    const AnalyticParameters &, const GeneralHybridOptions &, SurfaceSampleContext &);
// 一维齐次 Bézier 提取，同时向白盒测试开放。
std::array<std::vector<BezierCell>,4> ExtractHybridBoundaries(const sggk::BSplineSurface &, V origin);
}
