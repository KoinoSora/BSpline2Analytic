#include "BSpline2Analytic.h"
#include "AnalyticReconstruction.h"
#include "Numerics.h"
#include "PatchRecovery.h"
#include "RecognitionValidation.h"

namespace bspline2analytic
{
const char *StatusName(Status s)
{
    switch (s)
    {
    case Status::Success:
        return "Success";
    case Status::InvalidInput:
        return "InvalidInput";
    case Status::Degenerate:
        return "Degenerate";
    case Status::IllConditioned:
        return "IllConditioned";
    case Status::InconsistentBoundary:
        return "InconsistentBoundary";
    case Status::Ambiguous:
        return "Ambiguous";
    case Status::VerificationFailed:
        return "VerificationFailed";
    case Status::Unsupported:
        return "Unsupported";
    }
    return "Unknown";
}
const char *TypeName(sggk::SurfaceType t)
{
    switch (t)
    {
    case sggk::SurfaceType::Plane:
        return "Plane";
    case sggk::SurfaceType::Cylinder:
        return "Cylinder";
    case sggk::SurfaceType::Cone:
        return "Cone";
    case sggk::SurfaceType::Sphere:
        return "Sphere";
    case sggk::SurfaceType::Torus:
        return "Torus";
    default:
        return "Unknown";
    }
}
static bool Equivalent(const AnalyticParameters &a, const AnalyticParameters &b,
                       const RecognitionOptions &o)
{
    using namespace detail;
    if (a.type != b.type)
        return false;
    V delta = Get(a.origin) - Get(b.origin), axis = Get(a.axis);
    if (a.type == sggk::SurfaceType::Plane)
        return Parallel(axis, Get(b.axis), o.angularTolerance) &&
               std::abs(Dot(delta, axis)) <= o.distanceTolerance;
    if (a.type == sggk::SurfaceType::Cylinder)
        delta = delta - axis * Dot(delta, axis);
    return Norm(delta) <= o.distanceTolerance &&
           std::abs(a.radius - b.radius) <= o.distanceTolerance &&
           std::abs(a.majorRadius - b.majorRadius) <= o.distanceTolerance &&
           std::abs(a.halfAngle - b.halfAngle) <= o.angularTolerance &&
           (a.type == sggk::SurfaceType::Sphere || Parallel(axis, Get(b.axis), o.angularTolerance));
}
RecognitionResult Recognize(const sggk::BSplineSurface &s, const RecognitionOptions &o)
{
    RecognitionResult result;
    auto analysis = AnalyzeBoundaries(s, o);
    result.boundaries = analysis.boundaries;
    if (analysis.status != Status::Success)
    {
        result.status = analysis.status;
        result.message = analysis.message;
        return result;
    }
    try
    {
        auto candidates = detail::Reconstruct(result.boundaries, o);
        if (candidates.empty())
            throw detail::Failure(Status::InconsistentBoundary,
                                  "No analytic surface satisfies the boundary relationships");
        bool found = false;
        Status failureStatus = Status::VerificationFailed;
        std::string failures;
        for (const auto &c : candidates)
        {
            try
            {
                auto validation = detail::ValidateCandidate(s, c, result.boundaries, o);
                sggk::UVRange range;
                ParameterCorrespondence mapping;
                detail::RecoverPatch(s, c, result.boundaries, o, range, mapping);
                if (found)
                {
                    if (!Equivalent(result.parameters, c.parameters, o))
                        throw detail::Failure(
                            Status::Ambiguous,
                            "More than one inequivalent analytic candidate passed verification");
                    continue;
                }
                result.parameters = c.parameters;
                result.surface = c.surface;
                result.domain = range;
                result.mapping = mapping;
                result.validation = validation;
                found = true;
            }
            catch (const detail::Failure &e)
            {
                if (e.status == Status::Ambiguous)
                    throw;
                failureStatus = e.status;
                failures += std::string(TypeName(c.parameters.type)) + ": " + e.what() + "; ";
            }
        }
        if (!found)
            throw detail::Failure(failureStatus, failures);
        result.status = Status::Success;
        result.message = "Analytic surface and rectangular parameter range recovered";
    }
    catch (const detail::Failure &e)
    {
        result.status = e.status;
        result.message = e.what();
        result.surface.reset();
        result.domain.reset();
    }
    catch (const std::exception &e)
    {
        result.status = Status::InvalidInput;
        result.message = e.what();
        result.surface.reset();
        result.domain.reset();
    }
    return result;
}
} // namespace bspline2analytic
