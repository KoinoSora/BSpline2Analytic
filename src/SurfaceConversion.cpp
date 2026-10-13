#include "SurfaceConversion.h"
namespace bspline2analytic
{
SurfaceConversionResult Convert(const sggk::BSplineSurface &s, const SurfaceConversionOptions &o)
{
    SurfaceConversionResult r;
    if (s.IsRational())
    {
        r.path = ConversionPath::RationalRecognition;
        r.details = Recognize(s, o.rational);
    }
    else
    {
        r.path = ConversionPath::PolynomialApproximation;
        r.details = ConvertApproximate(s, o.polynomial);
    }
    return r;
}
bool SurfaceConversionResult::Succeeded() const
{
    if (const auto *r = std::get_if<RecognitionResult>(&details))
        return r->status == Status::Success;
    return std::get<ConversionResult>(details).status == ConversionStatus::Success;
}
bool SurfaceConversionResult::IsDegenerate() const
{
    if (const auto *r = std::get_if<RecognitionResult>(&details))
        return r->status == Status::Degenerate;
    return std::get<ConversionResult>(details).status == ConversionStatus::Degenerate;
}
const char *SurfaceConversionResult::StatusText() const
{
    if (const auto *r = std::get_if<RecognitionResult>(&details))
        return StatusName(r->status);
    return ConversionStatusName(std::get<ConversionResult>(details).status);
}
sggk::SurfaceType SurfaceConversionResult::Type() const
{
    if (const auto *r = std::get_if<RecognitionResult>(&details))
        return r->parameters.type;
    const auto &r = std::get<ConversionResult>(details);
    return r.patch ? r.patch->parameters.type : sggk::SurfaceType::Unknown;
}
const ConversionTiming &SurfaceConversionResult::Timing() const
{
    return std::visit([](const auto &r) -> const ConversionTiming & { return r.timing; }, details);
}
std::string SurfaceConversionResult::Message() const
{
    return std::visit([](const auto &r) { return r.message; }, details);
}
const char *ConversionPathName(ConversionPath p)
{
    return p == ConversionPath::RationalRecognition ? "rational" : "polynomial";
}
} // namespace bspline2analytic
