/**
 * @file BoundaryAnalysis.h
 * @brief 精确边界分析的内部入口声明。
 * 为解析重建提供边界几何与失败原因，保持源 Umin/Umax/Vmin/Vmax 顺序。
 */
#pragma once
#include "BSpline2Analytic.h"
namespace bspline2analytic::detail
{
BoundaryInfo AnalyzeCurve(const sggk::Curve3DPtr &curve, BoundaryId id,
                          const RecognitionOptions &o);
}
