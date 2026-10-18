/**
 * @file main.cpp
 * @brief SGK 初始化与精确解析识别的最小示例。
 * 构造输入后调用识别接口并展示结果；当前等参优先、失败兜底入口参见 SequentialConversion.h。
 */
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
    else if (const auto *approximate=std::get_if<ConversionResult>(&result.details))
        std::cout << "Finite-patch interval bound: "
                  << approximate->validation.geometry.upperBound
                  << '\n';
    else
        std::cout << "General finite-trim bound: "
                  << std::get<GeneralConversionResult>(result.details).validation.distanceUpperBound << '\n';
    if (result.patch)
        std::cout << "Linear PCurves: " << result.patch->PCurves().size()
                  << "; boundary Hausdorff bound: " << result.pcurveValidation.boundaryHausdorffUpperBound << '\n';
}

static int Run()
{
    // 两条边界圆即使同轴且所在平面平行，仍应识别为圆环面。
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
