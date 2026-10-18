/**
 * @file GeneralExtremeDataset.h
 * @brief 一般极端数据的来源几何合同声明。
 * 控制凸包、世界范围及参数域退化分别检查，不将空间极点等同于参数折叠。
 */
#pragma once
#include "GeneralDataset.h"
namespace bspline2analytic::tests
{
inline constexpr double GeneralExtremeWorldLimit=1e6;
struct GeneralExtremeSpec
{
    GeneralFixture fixture;
    std::string group, variant;
    double length=0, width=0, sag=0, translation=0;
    unsigned orientation=0;
};
std::vector<GeneralExtremeSpec> GeneralExtremeSpecs();
std::string GeneralExtremeFilename(unsigned);
void GenerateGeneralExtreme(const std::filesystem::path &,unsigned first=0,unsigned last=320);
// 采用实际存储控制点，并显式计入中心化区间余量；三角形可坍缩一条 UV 边，内部不可坍缩。
FixtureCertificate CertifyGeneralExtremePair(const GeneralFixture &,double,unsigned maxPieces=10000);
}
