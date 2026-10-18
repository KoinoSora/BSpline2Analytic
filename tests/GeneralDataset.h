/**
 * @file GeneralDataset.h
 * @brief 一般测试输入、独立真值与案例预期的共享声明。
 * V2 角点属于有理参数图；指定母面拒绝与排除所有目标母面的全局负例分别表示。
 */
#pragma once
#include "GeneralConversion.h"
#include "IntervalBounds.h"
#include <filesystem>

namespace bspline2analytic::tests
{
// 文件同时含输入几何和独立参考模型；只有输入传给 ConvertGeneral。
// 专项案例各有预期，不能统一标成无解。
enum class GeneralFixtureGroup { Feasible, Infeasible, NonlinearPCurve, Periodic,
    SpherePriority, Tolerance, InvalidOptions, ValidationBudget };
struct GeneralFixture
{
    unsigned id = 0, type = 0, bucket = 0, variant = 0, seed = 0;
    bool rational = false, extra = false;
    unsigned collapsedMask = 0;
    GeneralFixtureGroup group = GeneralFixtureGroup::Feasible;
    AnalyticParameters truth;
    std::array<sggk::Point2D, 4> corners;
    // 可选弯曲 UV 构造：v+=bend*s*(1-s)；原 500 个正例采用零弯曲。
    double bend = 0;
    sggk::BSplineSurfacePtr input;
    double storedUpper = 0, exclusionRadius = 0, storedRankMargin = 0, pairTolerance = 1e-6;
    // V2 的角点位于指定的母面有理参数图，不属于解析 UV。
    bool rationalDomain = false;
    std::array<double,4> chartDomain{};
    double storedTangentDeg = INFINITY;
};
struct FixtureCertificate
{
    bool passed = false;
    double upper = INFINITY, rankMargin = 0;
    unsigned pieces = 0, witnesses = 0;
    std::string message;
};
sggk::SurfaceType GeneralDatasetType(unsigned);
const char *FixtureGroupName(GeneralFixtureGroup);
std::string GeneralFixtureFilename(unsigned, bool extra = false);
GeneralFixture ReadGeneralFixture(const std::filesystem::path &);
void WriteGeneralFixture(const GeneralFixture &, const std::filesystem::path &);
void GenerateGeneralDataset(const std::filesystem::path &, unsigned first = 0, unsigned last = 660);
// 规范解析公式与参数图由测试侧独立实现，不使用生产投影代码。
detail::V FixtureTarget(const GeneralFixture &, double, double);
FixtureCertificate CertifyAnalyticPair(const GeneralFixture &, double tolerance,
                                      unsigned maxPieces = 20000);
// 排除所有总次数不超过 4 的非零多项式，因此覆盖五类目标母面。
FixtureCertificate CertifyNoQuartic(const sggk::BSplineSurface &, double radius);
unsigned FixtureCollapsedMask(const sggk::BSplineSurface &);
}
