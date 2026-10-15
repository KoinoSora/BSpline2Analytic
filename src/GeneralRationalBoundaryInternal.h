#pragma once
#include "GeneralRationalBoundaryValidation.h"
#include "GeneralHybridInternal.h"
#include "PairedValidationInternal.h"

namespace bspline2analytic::detail::rational_boundary
{
inline constexpr unsigned Ends[4][2]={{0,2},{1,3},{0,1},{2,3}};
struct Trim
{
    AnalyticParameters mother;
    std::array<double,4> domain{};
    paired::Chart chart;
    GeneralTrim polygon; // Only corners used: these are rational-chart coordinates.
    std::array<std::vector<std::pair<double,std::array<double,2>>>,4> walks;
    std::array<bool,4> collapsed{};
};
struct Target
{
    GeneralBoundaryPatch patch;
    std::array<std::vector<BezierCell>,4> edges; // Actual exported curves, centered.
    double motherBound=0;
};
struct CurveError
{
    BV numerator;
    BPoly denominator;
};
// One-dimensional homogeneous error, in a COMMON local parameter [0,1].
// Rebuild after changing the source/target correspondence.
class CurveErrors
{
    std::map<std::pair<unsigned,unsigned>,std::vector<I>> factors;
public:
    CurveError Build(const BezierCell &source,const BezierCell &target,paired::Products &,bool &fused);
};
BV TangentNumerator(const BezierCell &,paired::Products &);
// Sufficient same-direction angle bound from E'=T'-S'. INFINITY is inconclusive.
double ResidualTangentRatio(const CurveError &,const BezierCell &source,
                            const BV &sourceTangentNumerator,double distanceUpper);
// Test-visible, exact tensor-to-line composition in homogeneous Bernstein form.
BezierCell ComposeLine(const paired::TargetCell &,sggk::Point2D,sggk::Point2D,
                       double,double,paired::Products &);
std::vector<BezierCell> ComposeEdge(const Trim &,unsigned,paired::Products &);
Target BuildTarget(const Trim &,double outputTolerance,GeneralRationalBoundaryReport &);
bool Recover(SurfaceSampleContext &,const AnalyticParameters &,const GeneralRationalBoundaryOptions &,Trim &);
std::array<double,2> Coordinates(const Trim &,V,std::array<double,2> reference);
V ChartPoint(const Trim &,std::array<double,2>);
double Fraction(const Trim &,unsigned,V,double reference);
double FiniteDistance(const Trim &,V);
GeneralRationalBoundaryReport Validate(SurfaceSampleContext &,const AnalyticParameters &,
                                      const GeneralRationalBoundaryOptions &);
}
