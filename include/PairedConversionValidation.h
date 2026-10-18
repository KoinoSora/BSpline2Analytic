/**
 * @file PairedConversionValidation.h
 * @brief 独立配对验证器与可复用源面缓存接口。
 * 用有限参数图、分段单调对应及区间界验证位置、边界切向和法向；与母面距离验证合同不同。
 */
#pragma once
#include "ApproximateConversion.h"
#include <memory>

namespace bspline2analytic
{
// 独立实验认证接口，不调用投影验证器或优化器。
// 未成功认证不能推出几何上不存在可行解。
enum class PairedValidationStatus
{
    Passed,
    NotCertified,
    Unsupported,
    InvalidInput,
    BudgetExceeded,
    NumericalFailure
};
struct PairedValidationOptions
{
    ValidationOptions validation;
    bool enablePiecewiseMapping = true;
    unsigned maxMappingSegments = 512;
    unsigned maxMappingDepth = 12;
    unsigned maxCommonCells = 16384;
};
struct PairedParameterMap
{
    // 目标坐标属于经过方向调整和归一化的有理参数图，不是对解析角参数作线性插值。
    std::vector<double> source, target;
};
struct PairedValidationReport
{
    PairedValidationStatus status = PairedValidationStatus::NotCertified;
    // 此处 Exceeded 只针对报告中的参数对应关系，不排除其他对应关系可行。
    ConversionValidationReport certificate;
    PairedParameterMap uMap, vMap;
    unsigned mappingAttempts = 0, mappingPointEvaluations = 0;
    unsigned commonCells = 0, surfaceNodes = 0, boundaryNodes = 0;
    unsigned sourceExtractions = 0, targetPatches = 0;
    unsigned sourceCacheHits = 0, uExtractionMatrices = 0, vExtractionMatrices = 0;
    // certificate.timing.preparation 的子项，不可再次加到 totalMs。
    double nurbsMs = 0, mappingMs = 0;
    double sourceExtractionMs = 0, sourceCacheMs = 0, localTransformMs = 0;
    unsigned fastBoundPasses = 0, tightBoundChecks = 0;
    // certificate.timing.surfaceBounds 的子项，不包括细分耗时。
    double surfaceComposeMs = 0, surfaceProductMs = 0, surfaceDifferenceMs = 0;
    double surfaceFastBoundMs = 0, surfaceTightBoundMs = 0;
    // 计时结构版本 2：源/目标准备属于 preparation；边界复合、代数、求界与细分属于 boundaryChecks。
    unsigned preparationSchema = 2;
    unsigned sourceBoundaryBuilds = 0, boundarySegments = 0;
    unsigned recoveryTraceHits = 0, recoveryTraceBuilds = 0;
    double sourceCenterMs = 0, sourceBoundaryMs = 0, targetTransformMs = 0;
    double targetBoundaryMs = 0, edgeComposeMs = 0, edgeAlgebraMs = 0;
    double edgeBoundsMs = 0, edgeSubdivisionMs = 0;
    bool Passed() const noexcept
    {
        return status == PairedValidationStatus::Passed && certificate.Passed();
    }
};
// 同一源面连续验证多个候选时可复用的缓存。
// 逐项精确比较源面内容，改变时失效；不缓存候选几何、参数域、映射或证书。
// 非线程安全，并发调用须各用独立上下文。
class PairedValidationContext
{
  public:
    PairedValidationContext();
    ~PairedValidationContext();
    PairedValidationContext(PairedValidationContext &&) noexcept;
    PairedValidationContext &operator=(PairedValidationContext &&) noexcept;
    PairedValidationContext(const PairedValidationContext &) = delete;
    PairedValidationContext &operator=(const PairedValidationContext &) = delete;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    friend PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &,
                                                           const AnalyticPatch &,
                                                           PairedValidationContext &,
                                                           const PairedValidationOptions &);
};
/// 无缓存重载：为本次源面/候选建立独立上下文，返回有限参数对应下的完整配对证书。
PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &, const AnalyticPatch &,
                                                const PairedValidationOptions & = {});
/// 缓存重载：仅复用精确匹配的源面准备；候选、参数图及对应关系仍在每次调用中重新验证。
PairedValidationReport ValidatePairedConversion(const sggk::BSplineSurface &, const AnalyticPatch &,
                                                PairedValidationContext &,
                                                const PairedValidationOptions & = {});
const char *PairedValidationStatusName(PairedValidationStatus) noexcept;
} // namespace bspline2analytic
