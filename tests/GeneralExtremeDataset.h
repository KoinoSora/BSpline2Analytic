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
// Uses actual stored control points, with an explicit interval allowance for
// recentering. A triangle may collapse one UV side, never its interior.
FixtureCertificate CertifyGeneralExtremePair(const GeneralFixture &,double,unsigned maxPieces=10000);
}
