/**
 * @file GeneralRationalBoundaryInternal.h
 * @brief 当前有理边界恢复、输出与连续认证的内部数据契约。
 * Trim 表示有理参数图候选，Target 保存实际输出；候选恢复、数值过滤和最终证书互不混淆。
 */
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
    GeneralTrim polygon; // 仅使用角点，坐标属于有理参数图。
    std::array<std::vector<std::pair<double,std::array<double,2>>>,4> walks;
    std::array<bool,4> collapsed{};
    std::array<bool,4> refine3D{};
    unsigned refinementCalls = 0;
};
struct Target
{
    GeneralBoundaryPatch patch;
    std::array<std::vector<BezierCell>,4> edges; // 经过中心化的实际导出曲线。
    double motherBound=0;
};
struct CurveError
{
    BV numerator;
    BPoly denominator;
};
// 源/目标在同一局部参数 [0,1] 下的一维齐次误差；对应关系改变后必须重新构造。
class CurveErrors
{
    std::map<std::pair<unsigned,unsigned>,std::vector<I>> factors;
public:
    CurveError Build(const BezierCell &source,const BezierCell &target,paired::Products &,bool &fused);
};
// 双精度 Bernstein 系数配合统一向外舍入误差半径，仅供当前一般验证器过滤。
// 无结论的节点在细分前交给 CurveErrors 区间重算，仍共享同一节点预算。
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
// 由 E'=T'-S' 得到同向切向夹角的充分上界；INFINITY 表示无法判定。
double ResidualTangentRatio(const CurveError &,const BezierCell &source,
                            const BV &sourceTangentNumerator,double distanceUpper);
// 由点值估计选取正辅助比例，本身不是证书。
// 所有参数必须使用同一局部参数；映射细化后重新构造。
std::optional<BPoly> ChooseTangentScale(const BV &source,const BV &target);
BV BuildScaledTangentResidual(const BV &source,const BV &target,const BPoly &scale,paired::Products &);
V TangentResidualReference(const BV &source);
// 采用 ||A-lambda*B||² / lower(||A||²) 求界；lambda>0 且上界小于 1 还可证明 A.B>0。
// INFINITY 表示无结论，不代表已有夹角超限见证。
double ScaledTangentSineSquaredUpper(const BV &source,const BV &residual,const BPoly &scale);
// 供测试访问的齐次 Bernstein 张量曲面与参数直线的精确代数复合。
BezierCell ComposeLine(const paired::TargetCell &,sggk::Point2D,sggk::Point2D,
                       double,double,paired::Products &);
std::vector<BezierCell> ComposeEdge(const Trim &,unsigned,paired::Products &);
Target BuildTarget(const Trim &,double outputTolerance,GeneralRationalBoundaryReport &,
                   const ExecutionDeadline *deadline=nullptr);
// 可选诊断输出；常规验证不收集轨迹，也不因诊断改变恢复决策。
// 详细分析在性能计时区间之外运行。
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
V ChartPoint(const Trim &,std::array<double,2>,std::array<V,2> *derivatives=nullptr);
// 仅提出候选：沿同一有理参数域直线边界执行有界、三维距离下降的 Gauss-Newton 修正。
// 验收仍由输出误差样条决定。
double RefineFraction(const Trim &,unsigned,V,double);
bool ValidRationalPolygon(const std::array<sggk::Point2D,4> &,const std::array<bool,4> &);
// 完整笛卡尔距离界，支持将端点与微小目标段配对。
std::optional<GeneralBoundaryTangentExemption> MicroSegmentBound(const BezierCell &,
    const std::vector<BezierCell> &,double,double,double,double);
double Fraction(const Trim &,unsigned,V,double reference);
double FiniteDistance(const Trim &,V);
GeneralRationalBoundaryReport Validate(SurfaceSampleContext &,const AnalyticParameters &,
                                      const GeneralRationalBoundaryOptions &);
}
