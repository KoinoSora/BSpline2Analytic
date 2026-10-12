#include "ApproximateConversion.h"
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
    sggk::Point3DMatrix points(4, std::vector<sggk::Point3D>(4));
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            points[i][j] = sggk::Point3D(i, j, (i > 0 && i < 3 && j > 0 && j < 3) ? 2e-6 : 0);
    sggk::BSplineSurface polynomial(3, 3, points, sggk::RealArray{0, 1}, sggk::RealArray{0, 1},
                                    sggk::UIntArray{4, 4}, sggk::UIntArray{4, 4});
    auto approximate = bspline2analytic::ConvertApproximate(polynomial);
    const auto &t = approximate.timing;
    std::cout << "Approximate: " << bspline2analytic::ConversionStatusName(approximate.status)
              << "\nclassification: " << t.classification.elapsedMs << " ms"
              << "\noptimization: " << t.optimization.elapsedMs << " ms (" << t.optimization.calls
              << " calls)"
              << "\nvalidation: " << t.validation.elapsedMs << " ms (" << t.validation.calls
              << " calls)"
              << "\nother: " << t.otherMs << " ms; total: " << t.totalMs << " ms\n";
    return result.status == bspline2analytic::Status::Success &&
                   approximate.status == bspline2analytic::ConversionStatus::Success
               ? 0
               : 1;
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
