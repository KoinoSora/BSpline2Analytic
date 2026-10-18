/**
 * @file BSpline2Analytic.h
 * @brief 精确识别公共接口及共享几何、状态和计时类型。
 * 以四边直线/圆弧关系恢复五类解析母面，采样与 Bernstein 界的保证等级分别报告。
 */
#pragma once

#include <GeomBase/Dir3D.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <array>
#include <limits>
#include <optional>
#include <string>

namespace bspline2analytic
{
enum class Status
{
    Success,
    InvalidInput,
    Degenerate,
    IllConditioned,
    InconsistentBoundary,
    Ambiguous,
    VerificationFailed,
    Unsupported
};
enum class BoundaryId
{
    Umin,
    Umax,
    Vmin,
    Vmax
};
enum class BoundaryKind
{
    Unknown,
    Line,
    Circle,
    Degenerate
};
enum class ClassificationFailureReason
{
    None,
    DegenerateBoundary,
    NonLineOrCircleBoundary,
    InconsistentGeometry,
    NoCandidate,
    IllConditionedBoundary
};
/// 验证强度；Sampled 只说明检查点通过，BernsteinBound 才表示已完成对应连续误差求界。
enum class VerificationLevel
{
    None,
    Sampled,
    BernsteinBound
};

/// 精确识别配置。距离随模型单位，角度差以弧度计；平行性也用 angularTolerance 比较单位方向叉积模。
/// 数值秩阈值与几何容差独立，不能把它们当作边界转角的度数容差。
struct RecognitionOptions
{
    // 长度单位与输入 SDK 模型一致，通常为 mm。
    double distanceTolerance = 1e-5;
    double rankTolerance = 1e-12;
    double angularTolerance = 1e-8;
    unsigned maxSubdivisionDepth = 12;
    unsigned maxValidationPieces = 20000;
    unsigned samplesPerSpan = 4;
    bool requireBound = true;
};

struct LineData
{
    sggk::Point3D origin;
    sggk::Dir3D direction{1, 0, 0};
    double start = 0, end = 0;
};

struct CircleData
{
    sggk::Point3D center;
    sggk::Dir3D normal{0, 0, 1};
    sggk::Dir3D xDirection{1, 0, 0};
    double radius = 0, startAngle = 0, sweepAngle = 0;
};

struct BoundaryInfo
{
    BoundaryId id = BoundaryId::Umin;
    BoundaryKind kind = BoundaryKind::Unknown;
    sggk::Curve3DPtr curve;
    LineData line;
    CircleData circle;
    double errorBound = std::numeric_limits<double>::infinity();
    std::string message;
};

struct BoundaryAnalysisResult
{
    Status status = Status::InvalidInput;
    std::array<BoundaryInfo, 4> boundaries;
    std::string message;
    ClassificationFailureReason failureReason = ClassificationFailureReason::None;
};

// 解析母面的规范几何坐标架，不要求复现源曲面的参数原点。
struct AnalyticParameters
{
    sggk::SurfaceType type = sggk::SurfaceType::Unknown;
    sggk::Point3D origin; // 依曲面类型分别为平面原点、轴线上参考点、锥顶、球心或圆环中心。
    sggk::Dir3D axis{0, 0, 1};
    sggk::Dir3D xDirection{1, 0, 0};
    double radius = 0;      // 球面或圆柱半径；圆环面的小半径（管半径）。
    double majorRadius = 0; // SGK Lemon 曲面采用带符号的大半径。
    double halfAngle = 0;
};

struct ParameterCorrespondence
{
    bool swappedUV = false;
    bool reversedU = false;
    bool reversedV = false;
    // 这些标志只描述方向，不假设 NURBS 参数与解析角参数之间存在线性映射。
    // 角点顺序为 (umin,vmin)、(umax,vmin)、(umin,vmax)、(umax,vmax)。
    std::array<sggk::Point2D, 4> analyticCorners;
};

struct ValidationReport
{
    VerificationLevel level = VerificationLevel::None;
    double maxSampledDistance = 0;
    double maxBoundaryBound = 0;
    double surfaceDistanceBound = std::numeric_limits<double>::infinity();
    unsigned piecesVisited = 0;
    // BernsteinBound 包含浮点误差保护，但不等同于采用向外舍入的区间算术证书。
    std::string message;
};

struct StageTiming
{
    double elapsedMs = 0;
    unsigned calls = 0;
};
struct ConversionTiming
{
    // 互不重叠的阶段总计。
    // 专用路径中 classification 包含初值验证，optimization 包含优化后的验证；validation 记录阶段外检查（如 Convert 最终的 PCurve
    // 检查）。
    // 一般兜底沿用其阶段划分。
    StageTiming classification, optimization, validation;
    double otherMs = 0, totalMs = 0;
    // 专用路径的细分计时属于上述阶段总计的子项，不可再次累加；Convert 也会传递这些明细。
    StageTiming classificationCompute, classificationValidation;
    StageTiming optimizationCompute, optimizationValidation;
};

struct RecognitionResult
{
    Status status = Status::InvalidInput;
    AnalyticParameters parameters;
    sggk::SurfacePtr surface;
    std::array<BoundaryInfo, 4> boundaries;
    std::optional<sggk::UVRange> domain;
    ParameterCorrespondence mapping;
    ValidationReport validation;
    ConversionTiming timing;
    std::string message;
    ClassificationFailureReason classificationFailure = ClassificationFailureReason::None;
};

/** @brief 提取并分类源参数矩形的四边。
 * @param surface 调用期间保持不变的源 NURBS 曲面。
 * @param options 几何容差、采样和求界预算。
 * @return 按 Umin、Umax、Vmin、Vmax 排列的边界及诊断；边界通过不等于曲面转换成功。
 */
BoundaryAnalysisResult AnalyzeBoundaries(const sggk::BSplineSurface &surface,
                                         const RecognitionOptions &options = {});
/** @brief 从满足等参边界结构的源面识别解析母面与有限参数范围。
 * @return 仅 status==Success 时可消费解析几何；验证强度与误差界见 validation。
 * 多个不等价候选同时通过、预算或误差检查失败均保留原因，不推断所有解析解都不存在。
 */
RecognitionResult Recognize(const sggk::BSplineSurface &surface,
                            const RecognitionOptions &options = {});
const char *StatusName(Status status);
const char *TypeName(sggk::SurfaceType type);
} // namespace bspline2analytic
