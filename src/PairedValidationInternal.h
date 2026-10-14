#pragma once
#include "ApproximateInternal.h"
#include "PairedConversionValidation.h"
#include <map>

namespace bspline2analytic::detail::paired
{
struct TargetCell
{
    BezierCell surface;
    BV normal; // Canonical normal numerator, with patch.normalSign applied.
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
};
BPoly Restrict(const BPoly &, bool alongU, I lo, I hi);
BV Restrict(const BV &, bool alongU, I lo, I hi);
BezierCell Crop(const BezierCell &, double, double, double, double);
Chart BuildChart(const AnalyticPatch &);
void CheckPatch(const AnalyticPatch &);
double RationalCoordinate(const Chart &, unsigned dimension, double analyticParameter);
BV LocalSource(const BV &, const AnalyticParameters &);
} // namespace bspline2analytic::detail::paired
