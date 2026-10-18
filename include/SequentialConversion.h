/**
 * @file SequentialConversion.h
 * @brief 当前等参优先、失败兜底的统一限时转换接口。
 * 先执行有理精确识别或非有理母面认证，再按失败状态进入一般兜底；只以 Passed() 判定最终成功。
 */
#pragma once
#include "MotherDistanceValidation.h"
#include "GeneralRationalBoundaryValidation.h"

namespace bspline2analytic
{
/// 串联调用配置：总预算约束最终成功，等参子预算用于给一般兜底保留执行机会。
struct SequentialConversionOptions
{
    // 必须为正且不超过 50 ms；两条分支共享这一墙钟总预算。
    double maxTotalMs = 50;
    // 等参搜索未得结论时为兜底预留时间；子预算与总预算从同一时刻开始。
    // 子预算超限只代表等参分支失败。
    double maxIsoMs = 20;
    RecognitionOptions exact;
    ConversionOptions iso;
    MotherDistanceValidationOptions mother;
    GeneralRationalBoundaryOptions fallback;
};
/// 分支诊断与最终结果。即使某子报告曾通过，仍须以本层 Passed() 检查共享预算。
struct SequentialConversionResult
{
    bool accepted = false, isoAttempted = false, fallbackAttempted = false;
    bool timeBudgetExceeded = false;
    bool isoTimeBudgetExceeded = false;
    // path 为最后尝试的路径：iso_exact、iso_mother 或 general_rational，不单独表示成功。
    std::string path, message, isoMessage, timeoutStage;
    // 单位均为 ms；totalMs 包含失败尝试、调度与异常清理，不包括调用方日志/文件输出。
    double totalMs = 0, isoMs = 0, fallbackMs = 0;
    unsigned isoCandidates = 0, fallbackCandidates = 0, optimizationIterations = 0;
    RecognitionResult exact;
    std::optional<AnalyticPatch> isoPatch;
    MotherDistanceValidationReport mother;
    GeneralRationalBoundaryReport fallback;
    /// 最终成功判据；子预算耗尽后，兜底仍可能在总预算内成功。
    bool Passed() const noexcept { return accepted && !timeBudgetExceeded; }
};
// 当前验收合同：有理输入走精确识别；非有理等参输入采用母面距离与边界 U/V 联合切向认证；兜底采用内部采样与连续有理边界位置/切向认证。
// 不调用已废弃验证器。
// 兜底只获得输入几何，不注入测试真值或预设参数图。
/**
 * @brief 先等参识别/认证，失败后尝试一般有理边界兜底。
 * 输入曲面在调用期间须保持不变，长度单位跟随模型；不接收测试真值。
 * @return 包含实际耗时、失败阶段和对应分支报告的结果；仅 Passed() 为真时可消费几何。
 * @note 50 ms 是成功结果的验收上限，协作式取消不保证函数在该时刻被强制抢占返回。
 */
SequentialConversionResult ConvertIsoThenFallback(const sggk::BSplineSurface &,
                                                  const SequentialConversionOptions & = {});
}
