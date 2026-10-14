#include "ApproximateInternal.h"
#include "SurfaceConversion.h"
#include "Dataset500.h"
#include <Foundation/init.h>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
static unsigned checks = 0;
static void Require(bool condition, const std::string &message)
{
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
// x in [-a,a], y=h*x*x/(a*a). For h=R*theta^2/8 this is a
// polynomial approximation of a circle, without using the production fitter.
static sggk::BSplineCurve3DPtr Parabola(double a, double h)
{
    return std::make_shared<sggk::BSplineCurve3D>(2,
        sggk::Point3DArray{{-a,h,0},{0,-h,0},{a,h,0}},
        sggk::RealArray{0,1},sggk::UIntArray{3,3});
}
static sggk::BSplineSurfacePtr Extrude(double h)
{
    sggk::Point3DMatrix cp{{{-1,h,0},{-1,h,1}},{{0,-h,0},{0,-h,1}},{{1,h,0},{1,h,1}}};
    return std::make_shared<sggk::BSplineSurface>(2,1,cp,sggk::RealArray{0,1},
        sggk::RealArray{0,1},sggk::UIntArray{3,3},sggk::UIntArray{2,2});
}
static void Curves()
{
    for (auto [radius, theta] : {std::pair<double,double>{1e4,.0008}, {1e8,1e-5}})
    {
        auto curve = Parabola(radius*theta*.5,radius*theta*theta/8);
        IsoCurveFitDiagnostic d;
        auto fit = FitApproximateIsoCurve(curve,BoundaryId::Umin,2e-4,20000,&d);
        std::cout << "arc R=" << radius << " theta=" << theta << " error=" << d.circleError
                  << " ratio=" << d.transverseRatio << " state=" << CircleFitStateName(d.circleState) << '\n';
        Require(!fit.hasLine && fit.hasCircle,"Visible short arc was rejected or treated as a line");
        Require(d.circleState==CircleFitState::Accepted && std::abs(d.sweepAngle)<1e-3,
                "Short-arc fixture did not exercise the removed angle gate");
        Require(std::abs(d.radius/radius-1)<1e-3,"Short-arc recovered radius is inaccurate");
        if (theta<1e-4) Require(d.transverseRatio<1e-10,"Fixture did not exercise the covariance gate");
        double dense = 0;
        for (unsigned i=0;i<=256;++i)
        {
            V p=Get(curve->CalcPoint(double(i)/256))-Get(fit.circle.circle.center);
            V n=Get(fit.circle.circle.normal);double z=Dot(p,n);
            dense=std::max(dense,std::hypot(z,Norm(p-n*z)-fit.circle.circle.radius));
        }
        Require(dense<1e-5,"Accepted curve failed held-out geometric samples");
    }
    for (double h : {0.0,1e-16})
    {
        IsoCurveFitDiagnostic d;
        auto fit=FitApproximateIsoCurve(Parabola(1,h),BoundaryId::Umin,2e-4,20000,&d);
        Require(fit.hasLine && !fit.hasCircle && d.circleState==CircleFitState::NumericallyLinear,
                "Numerically straight input produced a spurious circle");
    }
    IsoCurveFitDiagnostic d;
    auto both=FitApproximateIsoCurve(Parabola(1,1e-4),BoundaryId::Umin,2e-4,20000,&d);
    Require(both.hasLine && both.hasCircle,"Compatible line and circle hypotheses were not retained");
    auto neither=FitApproximateIsoCurve(Parabola(1,1),BoundaryId::Umin,2e-4,20000,&d);
    Require(!neither.hasLine && !neither.hasCircle && d.circleState==CircleFitState::ResidualExceeded,
            "Noncircular parabola was accepted within tolerance");
    bool nonfinite=false;
    try { FitApproximateIsoCurve(Parabola(1,0),BoundaryId::Umin,2e-4,20000,nullptr,
                                [](double){ return V{INFINITY,0,0}; }); }
    catch (const Failure &e) { nonfinite=e.status==Status::IllConditioned; }
    Require(nonfinite,"Nonfinite evaluated data was not reported as a numerical failure");
}
static bool GeneralAttempted(const SurfaceConversionResult &r)
{
    return std::any_of(r.attempts.begin(),r.attempts.end(),[](const auto &a) {
        return a.path==ConversionPath::GeneralFallback && !a.skipped;
    });
}
static void StatesAndFallback()
{
    Require(ApproximateFailureStatus(Status::InconsistentBoundary)==ConversionStatus::NoCandidate,
            "Boundary failure became invalid input");
    for (Status s : {Status::IllConditioned,Status::Ambiguous,Status::VerificationFailed})
        Require(ApproximateFailureStatus(s)==ConversionStatus::Unresolved,"Unresolved failure lost its meaning");
    auto curved=Extrude(1);
    auto specialized=ConvertApproximate(*curved);
    Require(specialized.status==ConversionStatus::NoCandidate &&
            specialized.classificationFailure==ClassificationFailureReason::NonLineOrCircleBoundary,
            "Valid unsupported boundary did not return NoCandidate");
    SurfaceConversionOptions options;
    options.general.fitting.refinement.maxRounds=0;
    options.general.fitting.validation.maxValidationPieces=64;
    options.general.fitting.validation.maxSubdivisionDepth=4;
    auto r=Convert(*curved,options);
    Require(r.status!=SurfaceConversionStatus::InvalidInput && GeneralAttempted(r),
            "Valid boundary failure did not reach general fallback");
    Require(r.attempts.size()>=3 && !r.attempts[1].isoCurveFits.empty(),
            "General fallback erased specialized fit diagnostics");
    Require(r.attempts[1].isoCurveFits.back().circleState==CircleFitState::ResidualExceeded,
            "Wrong specialized rejection reason retained");
    // A real sampling-budget failure also must reach the general branch.
    options.polynomial.sampling.maxSamples=1;
    auto budget=Convert(*Extrude(0),options);
    Require(budget.status!=SurfaceConversionStatus::InvalidInput && GeneralAttempted(budget),
            "Sampling budget exhaustion was treated as invalid input");
    options.polynomial={};options.polynomial.validation.distanceTolerance=0;
    auto invalid=Convert(*Extrude(0),options);
    Require(invalid.status==SurfaceConversionStatus::InvalidInput && !GeneralAttempted(invalid),
            "Invalid options were forwarded to general fallback");
}
static void FormerInvalidInputs(const std::filesystem::path &data)
{
    for (unsigned id : {97u,99u,105u,107u,113u,115u,121u,123u})
    {
        auto c=ReadDatasetCase(data/("extreme-"+std::to_string(id+1000).substr(1)+".b2a"));
        bool shortCircle=false;
        for (unsigned direction=0;direction<2;++direction)
            for (bool upper : {false,true})
            {
                double t=direction==0 ? (upper?c.input->MaxParamU():c.input->MinParamU()) :
                                         (upper?c.input->MaxParamV():c.input->MinParamV());
                auto curve=direction==0 ? c.input->CalcVCurve(t) : c.input->CalcUCurve(t);
                IsoCurveFitDiagnostic d;
                auto fit=FitApproximateIsoCurve(curve,static_cast<BoundaryId>(2*direction+upper),2e-4,20000,&d);
                Require(fit.hasLine || fit.hasCircle,"Former invalid input still has no fitted boundary: "+std::to_string(id));
                shortCircle |= fit.hasCircle && std::abs(d.sweepAngle)<1e-3;
            }
        Require(shortCircle,"No short circle was retained: "+std::to_string(id));
    }
}
int main(int argc,char **argv)
{
    sggk::init();int code=1;
    try
    {
        if (argc!=2) throw std::runtime_error("Expected extreme200 fixture directory");
        std::cout << std::setprecision(17);
        Curves();StatesAndFallback();FormerInvalidInputs(argv[1]);
        std::cout << "Passed " << checks << " classification gate checks\n";code=0;
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; }
    sggk::fini();return code;
}
