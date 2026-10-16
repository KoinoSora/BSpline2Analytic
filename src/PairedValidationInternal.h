#pragma once
#include "ApproximateInternal.h"
#include "PairedConversionValidation.h"
#include <map>

namespace bspline2analytic::detail::mother { struct SourceCache; }

namespace bspline2analytic::detail::paired
{
struct TargetCell
{
    BezierCell surface;
    BV normal; // Canonical normal numerator, with patch.normalSign applied.
    // Analytically constructed factors in oriented source axes; never inferred
    // from approximate rank of a coefficient matrix.
    std::array<BPoly, 2> weightFactors;
    bool separableWeight = false;
};
struct Chart
{
    std::array<std::vector<double>, 2> angles;
    std::array<bool, 2> angular{};
    std::vector<TargetCell> cells;
};
// Private coefficient cache: no changes to shared interval arithmetic.
class Products
{
    std::map<std::pair<unsigned, unsigned>, std::vector<I>> factors;
    const std::vector<I> &Weights(unsigned, unsigned);

  public:
    BPoly Mul(const BPoly &, const BPoly &);
    BV Mul(const BV &, const BPoly &);
    BPoly Dot(const BV &, const BV &);
    BV Cross(const BV &, const BV &);
    BPoly Elevate(const BPoly &, unsigned, unsigned);
    BV Subtract(BV, const BV &);
    BV MulShared(const BV &, const BPoly &);
    BV MulWeight(const BV &, const TargetCell &);
};
double PositionBoxUpper(const BV &, const BPoly &);
double PositionControlUpper(const BV &, const BPoly &, Products &);
BPoly Restrict(const BPoly &, bool alongU, I lo, I hi);
BV Restrict(const BV &, bool alongU, I lo, I hi);
BezierCell Crop(const BezierCell &, double, double, double, double);
Chart BuildChart(const AnalyticPatch &);
void PlaceChart(Chart &, const AnalyticParameters &, V sourceOrigin);
void CheckPatch(const AnalyticPatch &);
double RationalCoordinate(const Chart &, unsigned dimension, double analyticParameter);
struct LocalFrame
{
    IV origin, x, y, z;
};
LocalFrame MakeLocalFrame(const AnalyticParameters &);
void Localize(BV &, const LocalFrame &);
std::vector<BezierCell> ExtractPolynomialCells(const sggk::BSplineSurface &,
                                               PairedValidationReport &);
struct BoundarySegment
{
    double begin = 0, end = 1;
    BV tangent, normal; // Packed into U; physical source parameter derivatives.
};
using Boundaries = std::array<std::vector<BoundarySegment>, 4>;
struct PreparedSource
{
    V origin;
    std::vector<BezierCell> cells;
    Boundaries boundaries;
};
PreparedSource PrepareSource(const sggk::BSplineSurface &, PairedValidationReport &);
struct BoundaryResult
{
    bool passed = false, tangentExceeded = false, normalExceeded = false;
    double tangentRatio = 0, normalRatio = 0;
};
// Maps are already finalized (including chart breaks); this never changes them.
Boundaries PrepareTargetBoundaries(const Chart &, Products &);
BoundaryResult ValidateBoundaries(const Boundaries &, const Boundaries &,
                                  const PairedParameterMap &, const PairedParameterMap &,
                                  const PairedValidationOptions &, PairedValidationReport &,
                                  Products &);
struct Budget : std::runtime_error
{
    using std::runtime_error::runtime_error;
};
// Source-owned immutable coefficients. No pointer-identity or lossy hash keys.
struct SourceCache
{
    unsigned degreeU = 0, degreeV = 0;
    std::array<double, 4> domain{};
    sggk::RealArray knotsU, knotsV;
    sggk::UIntArray multsU, multsV;
    std::vector<std::vector<V>> controls;
    PreparedSource prepared;
    bool Matches(const sggk::BSplineSurface &) const;
    void Store(const sggk::BSplineSurface &, PreparedSource &&);
};
// Stack-owned bridge between classification/recovery and paired validation.
// The source must remain unchanged for this object's lifetime.
struct AnalysisContext
{
    explicit AnalysisContext(const sggk::BSplineSurface &s, unsigned budget = 50000)
        : samples(s, budget)
    {
    }
    SurfaceSampleContext samples;
    std::optional<ClassificationSamples> classification;
    SourceCache source;
    // Independent polynomial cache used only by the mother-distance validator.
    std::shared_ptr<mother::SourceCache> motherSource;
};
PairedValidationReport ValidatePrepared(AnalysisContext &, const RecoveredCandidate &,
                                        const PairedValidationOptions & = {});
// Shared correspondence construction; no geometry certification is performed.
void CompleteChartMaps(PairedValidationReport &, const Chart &);
void RefineBoundaryMaps(AnalysisContext &, const AnalyticPatch &, const Chart &,
                        const PairedValidationOptions &, PairedValidationReport &,
                        const RecoveredCandidate *);
} // namespace bspline2analytic::detail::paired
