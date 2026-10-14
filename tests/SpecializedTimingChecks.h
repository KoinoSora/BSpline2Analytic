#pragma once
#include "BSpline2Analytic.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bspline2analytic::tests
{
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
