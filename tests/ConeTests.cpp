/**
 * @file ConeTests.cpp
 * @brief 圆锥初值、优化停滞与候选后备的历史专项回归。
 * 涉及旧近似验证器的完整执行已暂停，不能按历史性能命令重跑。
 */
#include "ApproximateCases.h"
#include "ApproximateInternal.h"
#include "Dataset500.h"
#include <Foundation/init.h>
#include <Geometry/3D/Surface/ConeSurface.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace
{
void Require(bool value, const char *message)
{
    if (!value)
        throw std::runtime_error(message);
}
bool SameParameters(const AnalyticPatch &a, const AnalyticPatch &b)
{
    return Norm(Get(a.parameters.origin) - Get(b.parameters.origin)) < 1e-12 &&
           Angle(Get(a.parameters.axis), Get(b.parameters.axis)) < 1e-12 &&
           std::abs(a.parameters.halfAngle - b.parameters.halfAngle) < 1e-12 &&
           a.domain == b.domain;
}
void Check249(const std::filesystem::path &directory, std::ofstream &csv)
{
    auto fixture = ReadDatasetCase(directory / CaseFilename(249));
    const auto &s = *fixture.input;
    ConversionOptions o;
    ClassificationDiagnostics oldCounts, newCounts;
    auto legacy = ClassifyApproximate(s, o, &oldCounts, false);
    auto candidates = ClassifyApproximate(s, o, &newCounts);
    Require(!legacy.empty() && !candidates.empty(), "Missing cone candidates");
    Require(newCounts.sampleCacheBuilds == 1 && newCounts.cachedPointEvaluations == 254 &&
                newCounts.cachedDerivativeEvaluations == 1 &&
                newCounts.directionalConeCandidates > 0 && newCounts.patchRecoveryAttempts <= 4 &&
                newCounts.rawCandidateCount > 16,
            "Classification cache or deferred recovery regressed");
    for (const auto &p : legacy)
        Require(std::any_of(candidates.begin(), candidates.end(),
                            [&](const AnalyticPatch &q) { return SameParameters(p, q); }),
                "Directional cone discarded a legacy fallback");
    auto samples = SampleSurface(s, o.sampling);
    EvaluationWorkspace w;
    auto p = legacy.front();
    ConversionStatistics stats;
    RefineCandidateLevel(p, samples, o, stats, ParameterLevel::SizeOrOffset, 6, false, w);
    Require(stats.refinementRuns.back().bestWorst > 1, "Old cone no longer exercises escalation");
    unsigned before = stats.optimizationIterations;
    o.refinement.recordIterationTrace = true;
    RefineCandidateLevel(p, samples, o, stats, ParameterLevel::MotherSurface, 10, true, w);
    const auto &run = stats.refinementRuns.back();
    Require(stats.optimizationIterations == before + 1 && run.acceptedSteps == 1 &&
                run.stopReason == RefinementStopReason::MarginReached && run.bestWorst <= .75 &&
                ValidateConversion(s, p).Passed(),
            "Forced mother solve did not restore the margin exit after its first best update");
    std::cout << "249 legacy: size iterations=" << before
              << " mother iterations=1 worst=" << run.bestWorst << '\n';
    // 即使采样已经可行，强制求解调用也不能通过零次尝试的裕量条件退出。
    for (unsigned mode = 0; mode < 3; ++mode)
    {
        auto feasible = p;
        ConversionOptions limited;
        if (mode == 0)
            limited.refinement.maxIterations = 0;
        if (mode == 1)
            limited.refinement.minTrustRadius = .1;
        EvaluationWorkspace scratch;
        ResidualBuffer initial;
        EvaluateResiduals(feasible, samples, limited.validation, scratch, initial);
        ConversionStatistics forced;
        RefineCandidateLevel(feasible, samples, limited, forced, ParameterLevel::MotherSurface, 10,
                             true, scratch);
        const auto &d = forced.refinementRuns.back();
        if (mode < 2)
            Require(forced.optimizationIterations == 0 &&
                        d.stopReason != RefinementStopReason::MarginReached &&
                        SameParameters(p, feasible),
                    "A bounded forced call falsely reported margin progress");
        else
            Require(forced.optimizationIterations > 0 &&
                        (d.stopReason != RefinementStopReason::MarginReached ||
                         d.bestWorst < initial.worst ||
                         (d.bestWorst == initial.worst && d.bestCost < initial.cost)),
                    "Forced margin exit returned an unchanged best candidate");
    }
    // 优先级最高的坐标架失败，不能移除下一项有效边界候选。
    auto cache = BuildClassificationSamples(s);
    auto broken = legacy.front().parameters;
    broken.axis = sggk::Dir3D(0, 0, 1);
    broken.xDirection = broken.axis;
    for (bool directional : {false, true})
    {
        std::vector<ApproximateSeed> seeds{{broken, 0, true, directional},
                                           {legacy.front().parameters, 1, true, false}};
        ClassificationDiagnostics recovery;
        auto found = RecoverRankedCandidates(cache, seeds, 1e-5, &recovery);
        Require(found.size() == 1 && recovery.patchRecoveryAttempts == 2 &&
                    SameParameters(found.front(), legacy.front()),
                "Failed preferred frame suppressed a usable fallback");
    }
    // 一次预热、五次计时；生产转换器始终不接收测试真值。
    for (int repeat = -1; repeat < 5; ++repeat)
    {
        auto result = ConvertApproximate(s);
        Require(result.status == ConversionStatus::Success && result.validation.Passed() &&
                    result.statistics.optimizationIterations == 0 &&
                    result.timing.optimization.calls == 0,
                "Directional cone 249 did not pass without refinement");
        if (repeat >= 0)
        {
            const auto &t = result.timing;
            const auto &d = result.statistics.classification;
            csv << "249," << repeat << ',' << t.classification.elapsedMs << ','
                << t.optimization.elapsedMs << ',' << t.validation.elapsedMs << ',' << t.totalMs
                << ',' << t.classification.elapsedMs + t.optimization.elapsedMs << ','
                << d.cachedPointEvaluations << ',' << d.patchRecoveryAttempts << ','
                << result.validation.geometry.upperBound << ','
                << result.validation.boundaryTangent.upperBound << ','
                << result.validation.boundaryNormal.upperBound << '\n';
        }
    }
    std::cout << "249 cache: points=" << newCounts.cachedPointEvaluations
              << " derivatives=" << newCounts.cachedDerivativeEvaluations
              << " raw=" << newCounts.rawCandidateCount
              << " recoveries=" << newCounts.patchRecoveryAttempts << std::endl;
}
void CheckPerturbedCone()
{
    sggk::Ucs3D frame(sggk::Point3D(13, -21, 8), sggk::Dir3D(1, 0, 1), sggk::Dir3D(-1, 2, 1));
    auto cone =
        std::make_shared<sggk::ConeSurface>(sggk::ConeSurface::CreateFromHalfAngle(frame, .4));
    auto source = cone->ToBSpline(sggk::UVRange(sggk::Interval(.2, 1.3), sggk::Interval(2, 5)));
    auto fixture = CurvedPlateauFixture({"cone", cone, source});
    Require(ValidateConversion(*fixture.spline, fixture.feasible).Passed(),
            "Known feasible perturbed cone failed validation");
    Require(!ValidateConversion(*fixture.spline, fixture.initial).Passed(),
            "Perturbed cone no longer needs optimization");
    for (auto mode : {JacobianMode::AnalyticWithFallback, JacobianMode::FiniteDifference})
    {
        ConversionOptions o;
        o.refinement.jacobianMode = mode;
        auto result = ConvertApproximate(*fixture.spline, o);
        Require(result.status == ConversionStatus::Success && result.validation.Passed() &&
                    result.statistics.optimizationIterations > 0,
                "Perturbed cone failed its optimization regression");
        std::cout << "Perturbed cone mode=" << int(mode)
                  << " iterations=" << result.statistics.optimizationIterations
                  << " bounds=" << result.validation.geometry.upperBound << ','
                  << result.validation.boundaryTangent.upperBound << ','
                  << result.validation.boundaryNormal.upperBound << std::endl;
    }
    ValidationOptions budget;
    budget.maxValidationPieces = 1;
    Require(!ValidateConversion(*fixture.spline, fixture.feasible, budget).Passed(),
            "Insufficient verification budget falsely accepted a cone");
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init();
    int code = 1;
    try
    {
        Require(argc == 3, "Usage: ConeTests dataset_directory timing.csv");
        std::ofstream csv(argv[2]);
        Require(bool(csv), "Cannot open cone timing CSV");
        csv << std::setprecision(17)
            << "id,run,classify_ms,optimize_ms,validate_ms,total_ms,classify_optimize_ms,"
               "cached_points,patch_recoveries,geometry_bound,tangent_deg,normal_deg\n";
        Check249(argv[1], csv);
        CheckPerturbedCone();
        Require(bool(csv), "Cone timing CSV write failed");
        std::cout << "Cone-focused regressions passed" << std::endl;
        code = 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << std::endl;
    }
    sggk::fini();
    return code;
}
