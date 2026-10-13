#pragma once
#include "ApproximateConversion.h"
#include <variant>

namespace bspline2analytic
{
enum class ConversionPath
{
    RationalRecognition,
    PolynomialApproximation
};
struct SurfaceConversionOptions
{
    RecognitionOptions rational;
    ConversionOptions polynomial;
};
// Preserve the two verification contracts rather than fabricating a shared certificate.
struct SurfaceConversionResult
{
    ConversionPath path = ConversionPath::RationalRecognition;
    std::variant<RecognitionResult, ConversionResult> details;
    bool Succeeded() const;
    bool IsDegenerate() const;
    const char *StatusText() const;
    sggk::SurfaceType Type() const;
    const ConversionTiming &Timing() const;
    std::string Message() const;
};
// Routing uses the SDK representation, not fixture provenance. Unit-weight planes
// therefore use the polynomial path. Neither path silently falls back to the other.
SurfaceConversionResult Convert(const sggk::BSplineSurface &input,
                                const SurfaceConversionOptions &options = {});
const char *ConversionPathName(ConversionPath path);
} // namespace bspline2analytic
