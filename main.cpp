#include "SurfaceConversion.h"
#include <Foundation/init.h>
#include <Geometry/3D/Surface/TorusSurface.h>
#include <iostream>

static void Print(const bspline2analytic::SurfaceConversionResult &result)
{
    using namespace bspline2analytic;
    const auto &t = result.Timing();
    std::cout << ConversionPathName(result.path) << ": " << result.StatusText() << ": "
              << TypeName(result.Type()) << '\n'
              << result.Message() << "\nclassification: " << t.classification.elapsedMs << " ms"
              << "\noptimization: " << t.optimization.elapsedMs << " ms"
              << "\nvalidation: " << t.validation.elapsedMs << " ms"
              << "\nother: " << t.otherMs << " ms; total: " << t.totalMs << " ms\n";
    if (const auto *exact = std::get_if<RecognitionResult>(&result.details))
        std::cout << "Bernstein mother-surface bound: " << exact->validation.surfaceDistanceBound
                  << '\n';
    else
        std::cout << "Finite-patch interval bound: "
                  << std::get<ConversionResult>(result.details).validation.geometry.upperBound
                  << '\n';
}

static int Run()
{
    // Concentric parallel boundary circles must still be recognized as a torus.
    sggk::TorusSurface source(sggk::Ucs3D(), 3, 1);
    auto spline = source.ToBSpline(
        sggk::UVRange(sggk::Interval(0.2, 1.1), sggk::Interval(0, 3.141592653589793)));
    spline->InsertUKnots((spline->MinParamU() + spline->MaxParamU()) * 0.5);
    auto result = bspline2analytic::Convert(*spline);
    Print(result);
    sggk::Point3DMatrix points(4, std::vector<sggk::Point3D>(4));
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            points[i][j] = sggk::Point3D(i, j, (i > 0 && i < 3 && j > 0 && j < 3) ? 2e-6 : 0);
    sggk::BSplineSurface polynomial(3, 3, points, sggk::RealArray{0, 1}, sggk::RealArray{0, 1},
                                    sggk::UIntArray{4, 4}, sggk::UIntArray{4, 4});
    auto approximate = bspline2analytic::Convert(polynomial);
    Print(approximate);
    return result.Succeeded() && approximate.Succeeded() ? 0 : 1;
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
