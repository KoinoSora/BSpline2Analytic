/**
 * @file SurfaceSampleContext.h
 * @brief 源面惰性采样上下文接口。
 * 源面必须在调用期间保持不变；缓存避免重复求值，但不缓存候选是否通过认证。
 */
#pragma once
#include "Numerics.h"
#include "ExecutionDeadline.h"
#include <map>
#include <tuple>
#include <functional>

namespace bspline2analytic::detail
{
struct SurfaceSample
{
    double u = 0, v = 0, area = 0;
    V point, du, dv, normal;
    bool normalValid = false;
    bool derivativesReady = false;
};
// 每次转换在栈上持有独立上下文，期间源面不得改变。
// 构造时不调用 SDK 求值，各采样视图和参数点均延迟计算。
class SurfaceSampleContext
{
public:
    explicit SurfaceSampleContext(const sggk::BSplineSurface &source, unsigned budget = 20000)
        : source(source), budget(budget) {}
    const SurfaceSample &At(double u, double v,
        sggk::DerivType sideU = sggk::DerivType::NONE,
        sggk::DerivType sideV = sggk::DerivType::NONE, bool needDerivatives = true);
    const std::vector<SurfaceSample> &Greville(bool needNormals = true);
    // 对已提取等参曲线延迟写入位置值；已有缓存优先，导数仍通过普通 At 路径求值。
    const SurfaceSample &AtCurvePoint(double u, double v, const std::function<V()> &evaluate,
                                     bool &inserted);
    std::vector<std::array<double, 2>> ScoreParameters() const;
    const std::vector<size_t> &Seeds();
    const std::vector<size_t> &Scores();
    const std::vector<size_t> &CoarseScores(unsigned count);
    const std::array<std::vector<SurfaceSample>, 4> &Boundaries();
    unsigned Evaluations() const { return static_cast<unsigned>(values.size()); }
    const sggk::BSplineSurface &source;
    // 可选的共享调用时间预算；为空时保留其他转换路径原有行为。
    const ExecutionDeadline *deadline = nullptr;
    void CheckDeadline(const char *stage) const { if (deadline) deadline->Check(stage); }
private:
    unsigned budget;
    bool grevilleReady = false, grevilleNormalsReady = false, boundariesReady = false;
    std::map<std::tuple<double, double, int, int>, SurfaceSample> values;
    std::vector<SurfaceSample> greville;
    std::vector<size_t> seeds, scores;
    std::map<unsigned, std::vector<size_t>> coarseScores;
    std::vector<size_t> seedPartition;
    std::vector<double> grevilleWeights;
    std::array<std::vector<SurfaceSample>, 4> boundaries;
};
RecognitionResult RecognizeWithSamples(const sggk::BSplineSurface &, const RecognitionOptions &,
                                      SurfaceSampleContext &);
}
