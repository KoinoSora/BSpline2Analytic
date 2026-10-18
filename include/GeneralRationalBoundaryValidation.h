/**
 * @file GeneralRationalBoundaryValidation.h
 * @brief 当前一般兜底验证的配置、边界输出与报告接口。
 * 约束位于母面 NURBS 参数域；内部仅采样，输出三维边界连续认证，不提供全域内部或直接法向证书。
 */
#pragma once
#include "GeneralHybridValidation.h"
#include <Geometry/3D/Curve/BSplineCurve3D.h>
#include <vector>

namespace bspline2analytic
{
struct GeneralRationalBoundaryOptions
{
    GeneralConversionOptions general;
    // 可选的固定解析参数盒，用于定义母面的有理参数图；不构造解析 PCurve。
    // 未指定时，依据有序边界轨迹推断参数盒。
    std::optional<std::array<double,4>> chartDomain;
    unsigned maxChartAttempts = 12, maxBoundaryNodes = 4096, maxBoundaryDepth = 14;
    // 包含所有验证阶段的墙钟预算，0 仅用于显式关闭诊断调用的限时。
    // 分类调用方应让所有候选共享同一个截止时间。
    double maxTotalMs = 50;
    // 只有源段和目标段各自的三维直径上界均不超过此值时才可合并；仍须满足 distanceTolerance 位置约束。
    // 对应切向明确记为豁免，不报告为角度认证。
    // 合法范围为 [0,1e-8]，0 表示关闭该策略。
    double microSegmentDiameter = 1e-8;
};
struct GeneralBoundaryTangentExemption
{
    unsigned edge = 0;
    // 源区间使用原曲面边界上的变化参数，目标区间使用 [0,1]。
    double sourceFirst = 0, sourceLast = 0, targetFirst = 0, targetLast = 0;
    double sourceDiameterUpper = 0, targetDiameterUpper = 0, distanceUpper = 0;
};
struct GeneralBoundaryPatch
{
    AnalyticParameters mother;
    // 边界存储顺序为 Umin、Umax、Vmin、Vmax，各自参数域为 [0,1]。
    // 闭合边环的正向遍历为 Vmin、Umax、-Vmax、-Umin。
    std::array<sggk::BSplineCurve3DPtr,4> boundaries;
    std::array<sggk::Point3D,4> vertices; // 角点顺序为 00、10、01、11，各输出曲线精确共享这些端点。
    std::array<bool,4> tangentNotApplicable{};
    std::vector<GeneralBoundaryTangentExemption> tangentExemptions;
    // 用于区域选择和构造追溯的元数据，不是解析参数域中的线性 PCurve。
    std::array<double,4> chartDomain{},sourceDomain{};
    std::array<sggk::Point2D,4> rationalCorners;
    sggk::Point3D interiorReference;
};
struct GeneralRationalBoundaryReport
{
    static constexpr bool certifiesInteriorContinuum = false;
    HybridValidationStatus status = HybridValidationStatus::NotCertified;
    std::optional<GeneralBoundaryPatch> patch;
    unsigned samples = 0, chartAttempts = 0, boundaryNodes = 0, boundarySegments = 0;
    unsigned mappingRefinements = 0, targetSegments = 0, maxTargetDegree = 0;
    unsigned microMergedSegments = 0, recovery3DRefinements = 0;
    unsigned fusedErrorBuilds = 0, positionOnlyNodes = 0;
    unsigned tangentFastChecks = 0, tangentFastAccepted = 0, tangentFallbacks = 0;
    // 正比例切向残差尝试次数，包括区间重算。
    unsigned tangentResidualChecks = 0, tangentResidualAccepted = 0;
    // 包括双精度包络关闭时直接采用区间计算的次数。
    unsigned tangentResidualIntervalFallbacks = 0;
    unsigned envelopeChecks = 0, envelopeAccepted = 0, envelopeFallbacks = 0;
    unsigned envelopePositionWitnesses = 0;
    unsigned envelopeTangentWitnesses = 0;
    double interiorSampleMaxDistance = 0, interiorSampleMaxFiniteDistance = 0;
    // 切向上界仅覆盖非坍缩且不在 patch.tangentExemptions 中的区间。
    double boundaryDistanceUpperBound = INFINITY, boundaryTangentUpperBoundDeg = INFINITY;
    // 输出曲线相对母面上理想有理轨迹的已认证距离上界。
    double outputMotherDistanceUpperBound = INFINITY;
    double samplingMs = 0, recoveryMs = 0, preparationMs = 0, boundaryMs = 0, totalMs = 0;
    bool timeBudgetExceeded = false;
    std::string timeoutStage;
    std::string message;
    bool Passed() const noexcept {return status==HybridValidationStatus::Passed&&patch.has_value();}
};
// 内部仅检查采样距离；输出边界进行连续位置认证，并在明确记录的微小段豁免之外进行连续切向认证。
// 不构造曲面法向系数，不检查法向转角。
/** @brief 对给定解析母面恢复有理参数域直线边界，并认证实际导出的三维 NURBS 边界。
 * @return 成功时含四边与母面参数；失败或超时不发布 patch，具体原因见 status/message。
 * @note 内部只提供采样证据；tangentExemptions 中的微小段仍认证位置，但不认证切向。
 */
GeneralRationalBoundaryReport ValidateGeneralRationalBoundaries(const sggk::BSplineSurface &,
    const AnalyticParameters &,const GeneralRationalBoundaryOptions & = {});
}
