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
    // Retain specialized fits even when details is replaced by general fallback.
    std::vector<IsoCurveFitDiagnostic> isoCurveFits;
};
struct SurfaceConversionOptions
{
    RecognitionOptions rational;
    ConversionOptions polynomial;
    GeneralConversionOptions general;
};
// Stage-specific reports stay in details; patch is the common finite result with
// explicit linear PCurves, published only on final success.
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
// Applicable specialized paths are attempted first; failures continue to general
// fitting. Sampling views are constructed only when first requested by a path.
SurfaceConversionResult Convert(const sggk::BSplineSurface &input,
                                const SurfaceConversionOptions &options = {});
const char *ConversionPathName(ConversionPath path);
} // namespace bspline2analytic
