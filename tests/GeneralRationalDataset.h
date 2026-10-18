/**
 * @file GeneralRationalDataset.h
 * @brief V2 测试参考参数图与来源证书接口。
 * 测试参考构造独立于生产图恢复，普通资格认证覆盖全域，不只检查边界。
 */
#pragma once
#include "GeneralDataset.h"

namespace bspline2analytic::tests
{
struct RationalSourceCertificate : FixtureCertificate
{
    double tangentDeg=0;
    double boundaryUpper=0;
};
struct RationalReferencePiece
{
    detail::BezierCell target;
    // 当前参数图公式适用的源参数包围矩形。
    std::array<double,4> domain;
};
// 测试侧独立参考构造，不调用生产参数图或恢复实现。
struct RationalReference
{
    GeneralFixture spec;
    std::vector<detail::BezierCell> chart;
    explicit RationalReference(const GeneralFixture &);
    detail::V Value(double,double) const;
    std::vector<detail::BezierCell> On(double,double,double,double,detail::V anchor) const;
    std::vector<RationalReferencePiece> Cover(const detail::BezierCell &,detail::V anchor) const;
    std::array<std::vector<double>,2> TensorBreaks() const;
};
// boundariesOnly 仅用于构造内部故意不可行的输入；普通数据资格认证必须使用默认 false。
RationalSourceCertificate CertifyRationalDomain(const GeneralFixture &,unsigned maxPieces=40000,bool boundariesOnly=false);
GeneralFixture RebuildRationalFixture(const GeneralFixture &,bool pressure=false);
std::string RationalFixtureName(const std::string &suite,unsigned id);
void GenerateRationalDataset(const std::filesystem::path &legacyMain,
    const std::filesystem::path &legacyExtreme,const std::filesystem::path &output,
    const std::string &suite,unsigned first,unsigned last);
void AuditRationalDataset(const std::filesystem::path &,const std::string &,unsigned,unsigned);
}
