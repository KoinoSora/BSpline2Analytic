/**
 * @file PairedValidationInternal.h
 * @brief 配对参数图、源缓存及候选分析上下文。
 * 源面准备可复用，候选几何和对应关系须重新认证；母面验证器使用独立源缓存。
 */
#pragma once
#include "ApproximateInternal.h"
#include "PairedConversionValidation.h"
#include <map>

namespace bspline2analytic::detail::mother { struct SourceCache; }

namespace bspline2analytic::detail::paired
{
struct TargetCell
{
    BezierCell surface;
    BV normal; // 规范法向分子已经应用 patch.normalSign。
    // 相关因子按源方向解析构造，不能从系数矩阵的近似秩推断。
    std::array<BPoly, 2> weightFactors;
    bool separableWeight = false;
};
struct Chart
{
    bool hasNormals = true;
    std::array<std::vector<double>, 2> angles;
    std::array<bool, 2> angular{};
    std::vector<TargetCell> cells;
};
// 私有系数缓存，不改变共享区间算术。
class Products
{
    std::map<std::pair<unsigned, unsigned>, std::vector<I>> factors;
    const std::vector<I> &Weights(unsigned, unsigned);

  public:
    BPoly Mul(const BPoly &, const BPoly &);
    BV Mul(const BV &, const BPoly &);
    BPoly Dot(const BV &, const BV &);
    BV Cross(const BV &, const BV &);
    BPoly Elevate(const BPoly &, unsigned, unsigned);
    BV Subtract(BV, const BV &);
    BV MulShared(const BV &, const BPoly &);
    BV MulWeight(const BV &, const TargetCell &);
};
double PositionBoxUpper(const BV &, const BPoly &);
double PositionControlUpper(const BV &, const BPoly &, Products &);
BPoly Restrict(const BPoly &, bool alongU, I lo, I hi);
BV Restrict(const BV &, bool alongU, I lo, I hi);
BezierCell Crop(const BezierCell &, double, double, double, double);
Chart BuildChart(const AnalyticPatch &, bool includeNormals = true);
void PlaceChart(Chart &, const AnalyticParameters &, V sourceOrigin);
void CheckPatch(const AnalyticPatch &, bool checkNormalOrientation = true);
double RationalCoordinate(const Chart &, unsigned dimension, double analyticParameter);
struct LocalFrame
{
    IV origin, x, y, z;
};
LocalFrame MakeLocalFrame(const AnalyticParameters &);
void Localize(BV &, const LocalFrame &);
std::vector<BezierCell> ExtractPolynomialCells(const sggk::BSplineSurface &,
                                               PairedValidationReport &);
struct BoundarySegment
{
    double begin = 0, end = 1;
    BV tangent, normal; // 压入 U 方向存储，导数仍相对于源面的实际参数。
};
using Boundaries = std::array<std::vector<BoundarySegment>, 4>;
struct PreparedSource
{
    V origin;
    std::vector<BezierCell> cells;
    Boundaries boundaries;
};
PreparedSource PrepareSource(const sggk::BSplineSurface &, PairedValidationReport &);
struct BoundaryResult
{
    bool passed = false, tangentExceeded = false, normalExceeded = false;
    double tangentRatio = 0, normalRatio = 0;
};
// 映射已定稿并包含参数图分段点，此函数不得再修改映射。
Boundaries PrepareTargetBoundaries(const Chart &, Products &);
BoundaryResult ValidateBoundaries(const Boundaries &, const Boundaries &,
                                  const PairedParameterMap &, const PairedParameterMap &,
                                  const PairedValidationOptions &, PairedValidationReport &,
                                  Products &, bool checkNormals = true);
struct Budget : std::runtime_error
{
    using std::runtime_error::runtime_error;
};
// 源面持有的不可变系数，不能仅按指针身份或有损哈希判断缓存是否有效。
struct SourceCache
{
    unsigned degreeU = 0, degreeV = 0;
    std::array<double, 4> domain{};
    sggk::RealArray knotsU, knotsV;
    sggk::UIntArray multsU, multsV;
    std::vector<std::vector<V>> controls;
    PreparedSource prepared;
    bool Matches(const sggk::BSplineSurface &) const;
    void Store(const sggk::BSplineSurface &, PreparedSource &&);
};
// 栈上上下文连接分类、恢复与配对验证；其生命周期内源面必须保持不变。
struct AnalysisContext
{
    explicit AnalysisContext(const sggk::BSplineSurface &s, unsigned budget = 50000)
        : samples(s, budget)
    {
    }
    SurfaceSampleContext samples;
    std::optional<ClassificationSamples> classification;
    SourceCache source;
    // 仅供母面距离验证器使用的独立多项式缓存。
    std::shared_ptr<mother::SourceCache> motherSource;
};
PairedValidationReport ValidatePrepared(AnalysisContext &, const RecoveredCandidate &,
                                        const PairedValidationOptions & = {});
// 共享的参数对应构造，不负责几何认证。
void CompleteChartMaps(PairedValidationReport &, const Chart &);
void RefineBoundaryMaps(AnalysisContext &, const AnalyticPatch &, const Chart &,
                        const PairedValidationOptions &, PairedValidationReport &,
                        const RecoveredCandidate *);
} // namespace bspline2analytic::detail::paired
