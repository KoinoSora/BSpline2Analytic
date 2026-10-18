/**
 * @file SurfaceConversion.h
 * @brief 保留的历史统一转换入口与分支报告类型。
 * 历史 Convert 仍连接旧近似验证链；当前共享 50 ms 预算入口为 ConvertIsoThenFallback，不能混用验收合同。
 */
#pragma once
#include "ApproximateConversion.h"
#include "GeneralConversion.h"
#include <variant>

namespace bspline2analytic
{
enum class ConversionPath
{
    RationalRecognition,
    PolynomialApproximation,
    GeneralFallback
};
enum class SurfaceConversionStatus { Success, InvalidInput, ConversionFailed };
struct ConversionAttempt
{
    ConversionPath path;
    bool skipped = false;
    bool succeeded = false;
    ClassificationFailureReason classificationFailure = ClassificationFailureReason::None;
    std::string message;
    // 即使 details 被一般兜底结果替换，仍保留专用分支的拟合诊断。
    std::vector<IsoCurveFitDiagnostic> isoCurveFits;
};
struct SurfaceConversionOptions
{
    RecognitionOptions rational;
    ConversionOptions polynomial;
    GeneralConversionOptions general;
};
// 各分支的报告保留在 details；patch 是含显式线性 PCurve 的统一有限面片结果，仅在最终成功时发布。
struct SurfaceConversionResult
{
    ConversionPath path = ConversionPath::RationalRecognition;
    std::variant<RecognitionResult, ConversionResult, GeneralConversionResult> details;
    SurfaceConversionStatus status = SurfaceConversionStatus::ConversionFailed;
    std::vector<ConversionAttempt> attempts;
    ConversionTiming timing;
    std::optional<GeneralAnalyticPatch> patch;
    GeneralValidationReport pcurveValidation;
    bool Succeeded() const;
    bool IsDegenerate() const;
    const char *StatusText() const;
    sggk::SurfaceType Type() const;
    const ConversionTiming &Timing() const;
    std::string Message() const;
};
// 历史统一入口先尝试适用的专用分支，失败后转一般拟合；各采样视图在首次请求时才构造。
// 当前限时串联入口另见 SequentialConversion.h。
SurfaceConversionResult Convert(const sggk::BSplineSurface &input,
                                const SurfaceConversionOptions &options = {});
const char *ConversionPathName(ConversionPath path);
} // namespace bspline2analytic
