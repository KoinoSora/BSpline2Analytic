/**
 * @file MotherDistanceValidation.h
 * @brief 非有理源面的母面距离与边界 U/V 联合切向验证接口。
 * 全域母面距离、对应边界位置和联合切向分别报告，不将联合充分条件混称为实际法向误差。
 */
#pragma once
#include "PairedConversionValidation.h"

namespace bspline2analytic
{
// 面向多项式源面的独立母面距离验证器。
// 全域距离采用精确识别配置与带浮点保护的 Bernstein 界；四边对应位置和两个参数切向采用区间认证。
// 不提供内部逐参数配对距离或有限裁剪面片的 Hausdorff 证书。
struct MotherDistanceValidationOptions
{
    RecognitionOptions position;
    // position.distanceTolerance 同时控制边界位置与参数对应。
    // boundary.validation 中旧的位置、切向和法向容差不用于本验证器验收，但其中的细分与映射预算仍然有效。
    PairedValidationOptions boundary;
    // 每个边界子区间均要求 sin²(angle(Su,Qu))+sin²(angle(Sv,Qv)) < sin²(tolerance)，且对应方向的点积为正。
    // 目标参数方向在解析构造上正交。
    // 不构造法向场，失败时也不回退至直接法向验证。
    double jointTangentToleranceDeg = 0.1;
};
struct MotherDistanceValidationReport
{
    PairedValidationStatus status = PairedValidationStatus::NotCertified;
    ConstraintReport motherDistance, boundaryPosition;
    ConstraintReport boundaryTangentU, boundaryTangentV, boundaryJointTangent;
    // 联合上界取自同一叶子区间的联合预算最大值，不能将 U/V 各自独立最大值相加替代。
    // boundaryJointTangent.upperBound 以度为单位，是 asin(sqrt(联合正弦平方预算)) 的上界，不是实测法向转角。
    double jointTangentSineSquaredUpperBound = INFINITY;
    // 通过或尚未检查时为空；稳定的失败代码区分有见证的联合条件超限、无法认证及预算耗尽。
    std::string jointTangentFailure;
    CheckState correspondence = CheckState::NotChecked;
    ValidationTiming timing;
    std::string message;
    PairedParameterMap uMap, vMap;
    unsigned surfaceNodes = 0, surfaceDepth = 0, boundaryPositionNodes = 0;
    unsigned boundaryNodes = 0, mappingAttempts = 0, sampleCount = 0;
    unsigned sourceExtractions = 0, sourceCacheHits = 0;
    // 正比例切向残差拟合次数（每个联合节点两次）与联合细分次数。
    unsigned tangentResidualFits = 0, tangentResidualSplits = 0;
    double sourceExtractionMs = 0, sourceBoundaryMs = 0, targetPreparationMs = 0;
    double mappingMs = 0, sampleMs = 0, residualBuildMs = 0, residualBoundMs = 0;
    double boundaryPositionMs = 0, boundaryAngleMs = 0; // 包含两个参数切向及联合预算的检查耗时。
    double sourceEnclosureRadius = 0;
    // 流式提取诊断；所有源面准备耗时均为 preparation 的子项。
    unsigned sourceCells = 0, extractedCells = 0, cachedCellHits = 0;
    unsigned extractionPlanBuilds = 0, tightenedCells = 0, boundaryFallbacks = 0;
    double extractionPlanMs = 0, sourceContractionMs = 0, sourceEnvelopeMs = 0;
    double maxExtractionError = 0;
    bool Passed() const noexcept
    {
        return status == PairedValidationStatus::Passed &&
               motherDistance.state == CheckState::Passed &&
               boundaryPosition.state == CheckState::Passed &&
               boundaryTangentU.state == CheckState::Passed &&
               boundaryTangentV.state == CheckState::Passed &&
               boundaryJointTangent.state == CheckState::Passed &&
               correspondence == CheckState::Passed;
    }
};
/** @brief 验证给定候选的全域母面距离、四边位置与 U/V 联合切向条件。
 * 输入为非有理、受支持次数的源面及有限候选；函数不优化几何。
 * @return Passed() 要求全部子检查通过；未认证不等价于不存在其他候选或对应关系。
 * @note 不检查相邻面片 G1，也不输出直接法向转角证书。
 */
MotherDistanceValidationReport ValidateMotherDistanceConversion(
    const sggk::BSplineSurface &, const AnalyticPatch &,
    const MotherDistanceValidationOptions & = {});
} // namespace bspline2analytic
