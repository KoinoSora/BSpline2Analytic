#pragma once
#include "ApproximateConversion.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bspline2analytic::tests
{
inline void CheckValidationTiming(const ValidationTiming &t)
{
    double sum = 0;
    for (const auto *stage : {&t.preparation, &t.surfaceBounds, &t.subdivision,
                             &t.boundaryChecks, &t.reporting})
    {
        if (!std::isfinite(stage->elapsedMs) || stage->elapsedMs < 0 ||
            (!stage->calls && stage->elapsedMs != 0))
            throw std::runtime_error("Invalid validation step timing");
        sum += stage->elapsedMs;
    }
    if (!std::isfinite(t.totalMs) || t.totalMs < 0 ||
        std::abs(sum - t.totalMs) > 1e-8 * std::max(1.0, t.totalMs))
        throw std::runtime_error("Validation step timing partition mismatch");
}
inline void CheckValidationTimings(const ConversionResult &r)
{
    CheckValidationTiming(r.validation.timing);
    const auto &initial = r.statistics.initialValidationTiming;
    const auto &refined = r.statistics.refinementValidationTiming;
    CheckValidationTiming(initial);
    CheckValidationTiming(refined);
    if (initial.preparation.calls != r.timing.classificationValidation.calls ||
        refined.preparation.calls != r.timing.optimizationValidation.calls ||
        initial.totalMs > r.timing.classificationValidation.elapsedMs + 1e-8 ||
        refined.totalMs > r.timing.optimizationValidation.elapsedMs + 1e-8 ||
        r.validation.timing.totalMs > initial.totalMs + refined.totalMs + 1e-8)
        throw std::runtime_error("Validation step aggregation mismatch");
}
// Check real converter timings, including failed/zero-iteration paths.
inline void CheckSpecializedTiming(const ConversionTiming &t)
{
    const auto close = [](double a, double b) {
        return std::abs(a - b) <= 1e-8 * std::max({1.0, std::abs(a), std::abs(b)});
    };
    for (const auto *stage : {&t.classification, &t.optimization, &t.validation,
                             &t.classificationCompute, &t.classificationValidation,
                             &t.optimizationCompute, &t.optimizationValidation})
        if (!std::isfinite(stage->elapsedMs) || stage->elapsedMs < 0 ||
            (!stage->calls && stage->elapsedMs != 0))
            throw std::runtime_error("Invalid specialized timing value");
    if (!close(t.classification.elapsedMs,
               t.classificationCompute.elapsedMs + t.classificationValidation.elapsedMs) ||
        !close(t.optimization.elapsedMs,
               t.optimizationCompute.elapsedMs + t.optimizationValidation.elapsedMs) ||
        t.classification.calls != t.classificationCompute.calls + t.classificationValidation.calls ||
        t.optimization.calls != t.optimizationCompute.calls + t.optimizationValidation.calls ||
        t.validation.calls != 0 || t.validation.elapsedMs != 0 ||
        !std::isfinite(t.totalMs) || !std::isfinite(t.otherMs) || t.otherMs < -1e-8 ||
        !close(t.totalMs, t.classification.elapsedMs + t.optimization.elapsedMs + t.otherMs))
        throw std::runtime_error("Specialized timing partition mismatch");
}
}
