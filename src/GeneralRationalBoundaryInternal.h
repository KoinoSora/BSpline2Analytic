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
// Double Bernstein coefficients plus a uniform outward-rounded error radius.
// Only the current general verifier uses this filter. Inconclusive nodes are
// recomputed by CurveErrors before subdivision, with the same node budget.
class CurveEnvelopes
{
    struct Polynomial
    {
        std::vector<double> c;
        double radius=0, magnitude=0;
    };
    struct Plan
    {
        std::vector<double> weights;
        double norm=0, error=0, gamma=0;
    };
    std::map<std::pair<unsigned,unsigned>,Plan> plans;
    const Plan &Weights(unsigned,unsigned);
    Polynomial Read(const BPoly &);
    Polynomial Mul(const Polynomial &,const Polynomial &);
    Polynomial Subtract(Polynomial,Polynomial);
    BPoly Enclose(const Polynomial &);
    I Midpoint(const BPoly &);
public:
    std::optional<CurveError> Build(const BezierCell &,const BezierCell &);
    std::optional<BV> Tangent(const BezierCell &);
    std::optional<BV> TangentResidual(const BV &source,const BV &target,const BPoly &scale);
    bool PositionExceeded(const CurveError &,double tolerance);
    bool TangentExceeded(const BV &,const BV &,double sineSquaredLimitUpper);
};
BV TangentNumerator(const BezierCell &,paired::Products &);
// Sufficient same-direction angle bound from E'=T'-S'. INFINITY is inconclusive.
double ResidualTangentRatio(const CurveError &,const BezierCell &source,
                            const BV &sourceTangentNumerator,double distanceUpper);
// Positive auxiliary scale selected from point estimates, never a certificate.
// All arguments use the SAME local parameter. Rebuild after mapping refinement.
std::optional<BPoly> ChooseTangentScale(const BV &source,const BV &target);
BV BuildScaledTangentResidual(const BV &source,const BV &target,const BPoly &scale,paired::Products &);
V TangentResidualReference(const BV &source);
// ||A-lambda B||^2 / lower(||A||^2); lambda>0 and a bound <1 also prove A.B>0.
// INFINITY means inconclusive, not a witnessed angle violation.
double ScaledTangentSineSquaredUpper(const BV &source,const BV &residual,const BPoly &scale);
// Test-visible, exact tensor-to-line composition in homogeneous Bernstein form.
BezierCell ComposeLine(const paired::TargetCell &,sggk::Point2D,sggk::Point2D,
                       double,double,paired::Products &);
std::vector<BezierCell> ComposeEdge(const Trim &,unsigned,paired::Products &);
Target BuildTarget(const Trim &,double outputTolerance,GeneralRationalBoundaryReport &,
                   const ExecutionDeadline *deadline=nullptr);
// Optional diagnostic sink. Normal validation does not collect walks or change
// any recovery decision. Detailed analysis is run outside performance timers.
struct RecoveryDiagnostics
{
    std::string reason="not_started";
    int edge=-1,sample=-1;
    double parameter=NAN,value=NAN,limit=NAN;
    std::array<std::vector<std::pair<double,std::array<double,2>>>,4> walks;
};
bool Recover(SurfaceSampleContext &,const AnalyticParameters &,const GeneralRationalBoundaryOptions &,Trim &,
             RecoveryDiagnostics *diagnostics=nullptr);
std::vector<AnalyticParameters> RecoveryCharts(SurfaceSampleContext &,const AnalyticParameters &,
                                              const GeneralRationalBoundaryOptions &);
std::array<double,2> Coordinates(const Trim &,V,std::array<double,2> reference);
V ChartPoint(const Trim &,std::array<double,2>);
double Fraction(const Trim &,unsigned,V,double reference);
double FiniteDistance(const Trim &,V);
GeneralRationalBoundaryReport Validate(SurfaceSampleContext &,const AnalyticParameters &,
                                      const GeneralRationalBoundaryOptions &);
}
