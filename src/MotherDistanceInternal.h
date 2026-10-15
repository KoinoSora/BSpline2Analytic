#pragma once
#include "MotherDistanceValidation.h"
#include "PairedValidationInternal.h"
#include "ValidationTiming.h"

namespace bspline2analytic::detail::mother
{
struct Measure
{
    double &elapsed;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    explicit Measure(double &v) : elapsed(v) {}
    ~Measure()
    {
        elapsed += std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - start).count();
    }
};
bool ValidatePosition(const paired::PreparedSource &, const AnalyticParameters &,
                      const RecognitionOptions &, MotherDistanceValidationReport &, ValidationTimer &);
bool ValidateEdgePositions(const paired::PreparedSource &, const paired::Chart &,
                           const PairedValidationOptions &, PairedValidationReport &,
                           paired::Products &, MotherDistanceValidationReport &);
MotherDistanceValidationReport ValidatePrepared(paired::AnalysisContext &,
    const RecoveredCandidate &, const MotherDistanceValidationOptions & = {});
} // namespace bspline2analytic::detail::mother
