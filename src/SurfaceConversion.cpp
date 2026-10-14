#include "SurfaceConversion.h"
#include "GeneralInternal.h"
#include "StageTimer.h"
namespace bspline2analytic
{
SurfaceConversionResult Convert(const sggk::BSplineSurface &s, const SurfaceConversionOptions &o)
{
    using namespace detail;
    const auto start=std::chrono::steady_clock::now();
    SurfaceConversionResult result;
    auto addTiming=[&](const ConversionTiming &timing) {
        for (auto pair:{std::make_pair(&result.timing.classification,&timing.classification),
                       std::make_pair(&result.timing.optimization,&timing.optimization),
                       std::make_pair(&result.timing.validation,&timing.validation)})
        { pair.first->elapsedMs+=pair.second->elapsedMs;pair.first->calls+=pair.second->calls; }
    };
    auto finish=[&] {
        result.timing.totalMs=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        result.timing.otherMs=result.timing.totalMs-result.timing.classification.elapsedMs-
            result.timing.optimization.elapsedMs-result.timing.validation.elapsedMs;
        if (result.status!=SurfaceConversionStatus::Success) result.patch.reset();
    };
    SurfaceSampleContext samples(s,std::max(o.polynomial.sampling.maxSamples,o.general.fitting.sampling.maxSamples));
    auto acceptSpecial=[&](const AnalyticParameters &parameters,double tolerance) {
        GeneralConversionOptions checks=o.general;
        checks.fitting.validation.distanceTolerance=tolerance;
        GeneralTrim trim;
        {
            StageTimer timer(result.timing.validation);
            result.pcurveValidation=ValidateGeneral(s,parameters,checks,samples,&trim);
        }
        if (!result.pcurveValidation.Passed()) return false;
        result.patch.emplace(trim.parameters,
            std::array<double,4>{s.MinParamU(),s.MaxParamU(),s.MinParamV(),s.MaxParamV()},trim.corners,trim.edges);
        result.status=SurfaceConversionStatus::Success;
        return true;
    };
    try
    {
        ValidateSurface(s);
        if (s.IsRational())
        {
            result.path=ConversionPath::RationalRecognition;
            auto exact=RecognizeWithSamples(s,o.rational,samples);
            addTiming(exact.timing);
            bool success=exact.status==Status::Success && acceptSpecial(exact.parameters,o.rational.distanceTolerance);
            if (exact.status==Status::Success && !success)
            { exact.status=Status::VerificationFailed;exact.message=result.pcurveValidation.message; }
            result.attempts.push_back({result.path,false,success,exact.classificationFailure,exact.message});
            bool invalid=exact.status==Status::InvalidInput;
            result.details=std::move(exact);
            if (success || invalid)
            {
                if (invalid) result.status=SurfaceConversionStatus::InvalidInput;
                finish();return result;
            }
            result.attempts.push_back({ConversionPath::PolynomialApproximation,true,false,
                ClassificationFailureReason::None,"Non-rational isocurve module is inapplicable to rational input"});
        }
        else
        {
            result.attempts.push_back({ConversionPath::RationalRecognition,true,false,
                ClassificationFailureReason::None,"Exact rational path skipped for polynomial representation"});
            result.path=ConversionPath::PolynomialApproximation;
            auto approximate=ConvertApproximateWithSamples(s,o.polynomial,samples);
            addTiming(approximate.timing);
            bool success=approximate.status==ConversionStatus::Success && approximate.patch &&
                acceptSpecial(approximate.patch->parameters,o.polynomial.validation.distanceTolerance);
            if (approximate.status==ConversionStatus::Success && !success)
            { approximate.status=ConversionStatus::VerificationFailed;approximate.message=result.pcurveValidation.message; }
            result.attempts.push_back({result.path,false,success,approximate.classificationFailure,approximate.message});
            result.attempts.back().isoCurveFits = approximate.statistics.classification.isoCurveFits;
            bool invalid=approximate.status==ConversionStatus::InvalidInput;
            result.details=std::move(approximate);
            if (success || invalid)
            {
                if (invalid) result.status=SurfaceConversionStatus::InvalidInput;
                finish();return result;
            }
        }
        result.path=ConversionPath::GeneralFallback;
        auto general=ConvertGeneralWithSamples(s,o.general,samples);
        addTiming(general.timing);
        bool success=general.status==ConversionStatus::Success;
        result.attempts.push_back({result.path,false,success,ClassificationFailureReason::None,general.message});
        result.pcurveValidation=general.validation;result.patch=general.patch;
        result.status=success ? SurfaceConversionStatus::Success :
            general.status==ConversionStatus::InvalidInput ? SurfaceConversionStatus::InvalidInput : SurfaceConversionStatus::ConversionFailed;
        result.details=std::move(general);
    }
    catch (const Failure &e)
    {
        result.status=e.status==Status::InvalidInput ? SurfaceConversionStatus::InvalidInput : SurfaceConversionStatus::ConversionFailed;
        GeneralConversionResult failure;
        failure.status=result.status==SurfaceConversionStatus::InvalidInput ? ConversionStatus::InvalidInput : ConversionStatus::VerificationFailed;
        failure.message=e.what();result.details=std::move(failure);
    }
    catch (const std::exception &e)
    {
        result.status=SurfaceConversionStatus::InvalidInput;
        GeneralConversionResult failure;failure.status=ConversionStatus::InvalidInput;failure.message=e.what();result.details=std::move(failure);
    }
    finish();return result;
}
bool SurfaceConversionResult::Succeeded() const { return status==SurfaceConversionStatus::Success && patch.has_value(); }
bool SurfaceConversionResult::IsDegenerate() const
{
    if (const auto *r=std::get_if<RecognitionResult>(&details)) return r->status==Status::Degenerate;
    if (const auto *r=std::get_if<ConversionResult>(&details)) return r->status==ConversionStatus::Degenerate;
    return std::get<GeneralConversionResult>(details).status==ConversionStatus::Degenerate;
}
const char *SurfaceConversionResult::StatusText() const
{
    return status==SurfaceConversionStatus::Success ? "Success" :
        status==SurfaceConversionStatus::InvalidInput ? "InvalidInput" : "ConversionFailed";
}
sggk::SurfaceType SurfaceConversionResult::Type() const
{
    return patch ? patch->Parameters().type : sggk::SurfaceType::Unknown;
}
const ConversionTiming &SurfaceConversionResult::Timing() const { return timing; }
std::string SurfaceConversionResult::Message() const
{
    return std::visit([](const auto &r) { return r.message; },details);
}
const char *ConversionPathName(ConversionPath p)
{
    switch (p)
    {
    case ConversionPath::RationalRecognition:return "rational";
    case ConversionPath::PolynomialApproximation:return "polynomial";
    case ConversionPath::GeneralFallback:return "general";
    }
    return "unknown";
}
}
