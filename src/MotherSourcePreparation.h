#pragma once
#include "MotherDistanceValidation.h"
#include "PairedValidationInternal.h"
#include "MotherTangentValidation.h"

namespace bspline2analytic::detail::mother
{
// Positive scalar error arithmetic. Coordinates themselves use ordinary doubles.
double ErrorAdd(double, double);
double ErrorMul(double, double);
double Gamma(unsigned);
double RoundedError(double);
double ErrorNorm(V);
struct SourceCell
{
    Poly x, y, z;
    V error;
    std::array<double, 4> domain{};
    std::optional<BezierCell> tight;
};
struct ExtractionSpan
{
    unsigned first = 0, degree = 0;
    double lo = 0, hi = 1;
    std::vector<I> enclosure;
    std::vector<double> weights, rowNorm, rowError;
};
struct SourceStrip
{
    std::vector<V> points, error;
};
struct SourceCache
{
    unsigned p = 0, q = 0;
    std::array<double, 4> domain{};
    sggk::RealArray knotsU, knotsV;
    sggk::UIntArray multsU, multsV;
    std::vector<std::vector<V>> controls;
    std::vector<V> centered;
    V origin, initialError, magnitude;
    std::vector<ExtractionSpan> us, vs;
    std::vector<std::optional<SourceStrip>> strips;
    std::vector<std::optional<SourceCell>> cells;
    std::optional<paired::PreparedSource> boundary;
    BoundaryTangents boundaryTangents;
    bool Matches(const sggk::BSplineSurface &) const;
    void Initialize(const sggk::BSplineSurface &, MotherDistanceValidationReport &);
    const SourceCell &Get(size_t, MotherDistanceValidationReport &);
    const SourceCell &Tighten(size_t, MotherDistanceValidationReport &);
    const paired::PreparedSource &PrepareBoundary(MotherDistanceValidationReport &,
                                                  bool tight = false);
};
} // namespace bspline2analytic::detail::mother
