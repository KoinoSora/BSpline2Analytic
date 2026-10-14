#pragma once
#include "Numerics.h"
#include <map>
#include <tuple>

namespace bspline2analytic::detail
{
struct SurfaceSample
{
    double u = 0, v = 0, area = 0;
    V point, du, dv, normal;
    bool normalValid = false;
    bool derivativesReady = false;
};
// Stack-owned per conversion. The source must not change during the conversion.
// Construction performs no SDK evaluation; each view and parameter pair is lazy.
class SurfaceSampleContext
{
public:
    explicit SurfaceSampleContext(const sggk::BSplineSurface &source, unsigned budget = 20000)
        : source(source), budget(budget) {}
    const SurfaceSample &At(double u, double v,
        sggk::DerivType sideU = sggk::DerivType::NONE,
        sggk::DerivType sideV = sggk::DerivType::NONE, bool needDerivatives = true);
    const std::vector<SurfaceSample> &Greville(bool needNormals = true);
    const std::vector<size_t> &Seeds();
    const std::vector<size_t> &Scores();
    const std::vector<size_t> &CoarseScores(unsigned count);
    const std::array<std::vector<SurfaceSample>, 4> &Boundaries();
    unsigned Evaluations() const { return static_cast<unsigned>(values.size()); }
    const sggk::BSplineSurface &source;
private:
    unsigned budget;
    bool grevilleReady = false, grevilleNormalsReady = false, boundariesReady = false;
    std::map<std::tuple<double, double, int, int>, SurfaceSample> values;
    std::vector<SurfaceSample> greville;
    std::vector<size_t> seeds, scores;
    std::map<unsigned, std::vector<size_t>> coarseScores;
    std::vector<size_t> seedPartition;
    std::vector<double> grevilleWeights;
    std::array<std::vector<SurfaceSample>, 4> boundaries;
};
RecognitionResult RecognizeWithSamples(const sggk::BSplineSurface &, const RecognitionOptions &,
                                      SurfaceSampleContext &);
}
