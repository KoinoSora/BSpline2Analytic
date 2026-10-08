#include "BSpline2Analytic.h"
#include <Foundation/init.h>
#include <Geometry/3D/Surface/TorusSurface.h>
#include <iostream>

static int Run()
{
    // Concentric parallel boundary circles must still be recognized as a torus.
    sggk::TorusSurface source(sggk::Ucs3D(), 3, 1);
    auto spline = source.ToBSpline(
        sggk::UVRange(sggk::Interval(0.2, 1.1), sggk::Interval(0, 3.141592653589793)));
    spline->InsertUKnots((spline->MinParamU() + spline->MaxParamU()) * 0.5);
    auto result = bspline2analytic::Recognize(*spline);
    std::cout << bspline2analytic::StatusName(result.status) << ": "
              << bspline2analytic::TypeName(result.parameters.type) << '\n'
              << result.message << "\nsampled distance: " << result.validation.maxSampledDistance
              << "\nBernstein mother-surface bound: " << result.validation.surfaceDistanceBound
              << '\n';
    return result.status == bspline2analytic::Status::Success ? 0 : 1;
}
int main()
{
    sggk::init();
    int result = 1;
    try
    {
        result = Run();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
    sggk::fini();
    return result;
}
