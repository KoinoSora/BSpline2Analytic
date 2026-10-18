/**
 * @file GeneralHybridValidation.h
 * @brief 历史混合验证器及当前复用的共享状态类型。
 * ValidateGeneralHybrid 已废弃，禁止作为回归或性能对照；当前边界验证见 GeneralRationalBoundaryValidation.h。
 */
#pragma once
#include "GeneralConversion.h"

namespace bspline2analytic
{
// 历史独立实验验证器：内部采用采样验收，不提供全曲面的连续证书；自身不运行优化器或更早的验证器。
// 当前调用状态以本文件末尾的废弃说明为准。
enum class HybridValidationStatus { Passed, NotCertified, InvalidInput, Unsupported, BudgetExceeded };
struct GeneralHybridOptions
{
    GeneralConversionOptions general;
    unsigned maxBoundaryPieces = 2048;
    unsigned maxBoundaryDepth = 12;
};
struct GeneralHybridReport
{
    static constexpr bool certifiesInteriorContinuum = false;
    HybridValidationStatus status = HybridValidationStatus::NotCertified;
    CheckState interiorSamples = CheckState::NotChecked;
    CheckState boundaryPosition = CheckState::NotChecked;
    CheckState boundaryTangent = CheckState::NotChecked;
    // 成功验证时，true 表示源边界精确坍缩为一点；其切向无定义，因此豁免切向检查，不能解释为已认证零角误差。
    std::array<bool,4> tangentNotApplicable{};
    // 特意避免命名为 distanceUpperBound，因为该值不对内部连续域作出保证。
    double interiorSampleMaxDistance = 0;
    double interiorSampleMaxFiniteDistance = 0;
    double boundaryDistanceUpperBound = INFINITY;
    double boundaryTangentUpperBoundDeg = INFINITY;
    unsigned interiorSampleCount = 0, boundaryNodes = 0, boundarySegments = 0;
    unsigned mappingRefinements = 0, collapsedEdges = 0, chartAttempts = 0, chartIterations = 0;
    unsigned failedBoundary = 0;
    double samplingMs = 0, recoveryMs = 0, preparationMs = 0, boundaryMs = 0, totalMs = 0;
    std::optional<GeneralAnalyticPatch> patch;
    std::string message;
    bool Passed() const noexcept
    {
        return status == HybridValidationStatus::Passed && patch.has_value() &&
               interiorSamples == CheckState::Passed && boundaryPosition == CheckState::Passed &&
               boundaryTangent == CheckState::Passed;
    }
};
// 已于 2026-10-15 废弃。
// 保留实现，不运行测试或性能基准。
// 当前验证器为 ValidateGeneralRationalBoundaries；共享状态类型仍可使用。
GeneralHybridReport ValidateGeneralHybrid(const sggk::BSplineSurface &, const AnalyticParameters &,
                                         const GeneralHybridOptions & = {});
const char *HybridValidationStatusName(HybridValidationStatus) noexcept;
}
