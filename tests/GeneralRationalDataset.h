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
    // Enclosing source rectangle where this chart formula applies.
    std::array<double,4> domain;
};
// Test-owned reference construction. Does not call production chart/recovery code.
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
// boundariesOnly is for constructing an intentionally infeasible interior;
// ordinary qualification always uses its default false.
RationalSourceCertificate CertifyRationalDomain(const GeneralFixture &,unsigned maxPieces=40000,bool boundariesOnly=false);
GeneralFixture RebuildRationalFixture(const GeneralFixture &,bool pressure=false);
std::string RationalFixtureName(const std::string &suite,unsigned id);
void GenerateRationalDataset(const std::filesystem::path &legacyMain,
    const std::filesystem::path &legacyExtreme,const std::filesystem::path &output,
    const std::string &suite,unsigned first,unsigned last);
void AuditRationalDataset(const std::filesystem::path &,const std::string &,unsigned,unsigned);
}
