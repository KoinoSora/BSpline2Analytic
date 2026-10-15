#include "ApproximateCases.h"
#include "Dataset500.h"
#include "PairedValidationInternal.h"
#include <Foundation/init.h>
#include <cfenv>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace bspline2analytic;
using namespace bspline2analytic::detail;
using namespace bspline2analytic::tests;
namespace
{
unsigned checks = 0;
void Require(bool b, const std::string &s)
{
    ++checks;
    if (!b)
        throw std::runtime_error(s);
}
std::string Quote(std::string s)
{
    std::string q = "\"";
    for (char c : s)
    {
        if (c == '\"')
            q += '\"';
        q += c;
    }
    return q + '\"';
}
void Timing(const PairedValidationReport &r)
{
    const auto &t = r.certificate.timing;
    double sum = t.preparation.elapsedMs + t.surfaceBounds.elapsedMs + t.subdivision.elapsedMs +
                 t.boundaryChecks.elapsedMs + t.reporting.elapsedMs;
    Require(std::isfinite(t.totalMs) && t.totalMs >= 0 &&
                std::abs(sum - t.totalMs) < 1e-7 * std::max(1.0, t.totalMs),
            "Exclusive paired timing mismatch");
    Require(r.nurbsMs + r.mappingMs + r.sourceExtractionMs + r.sourceCacheMs + r.localTransformMs +
                    r.sourceCenterMs + r.sourceBoundaryMs + r.targetTransformMs +
                    r.targetBoundaryMs <=
                t.preparation.elapsedMs + 1e-5,
            "Nested preparation timings exceed preparation");
    Require(r.surfaceComposeMs + r.surfaceProductMs + r.surfaceDifferenceMs + r.surfaceFastBoundMs +
                    r.surfaceTightBoundMs <=
                t.surfaceBounds.elapsedMs + 1e-5,
            "Nested position timings exceed surface bounds");
    Require(r.edgeComposeMs + r.edgeAlgebraMs + r.edgeBoundsMs + r.edgeSubdivisionMs <=
                t.boundaryChecks.elapsedMs + 1e-5,
            "Nested boundary timings exceed boundary checks");
    Require(r.fastBoundPasses + r.tightBoundChecks == r.surfaceNodes,
            "Position bound counters do not cover visited nodes");
    if (r.Passed())
        Require(r.certificate.Passed() && std::isfinite(r.certificate.geometry.upperBound),
                "Incomplete certificate accepted");
}
void CheckExtraction(const sggk::BSplineSurface &s)
{
    CheckApproximateInput(s, {});
    auto reference = ExtractCells(s);
    PairedValidationReport counts;
    auto cells = paired::ExtractPolynomialCells(s, counts);
    Require(cells.size() == reference.size() &&
                cells.size() == size_t(counts.uExtractionMatrices) * counts.vExtractionMatrices,
            "Extraction coverage or one-matrix-per-span count changed");
    PairedValidationReport preparation;
    auto prepared = paired::PrepareSource(s, preparation);
    Require(preparation.sourceBoundaryBuilds == 1, "Source boundaries not prepared exactly once");
    for (unsigned edge = 0; edge < 4; ++edge)
        for (const auto &segment : prepared.boundaries[edge])
        {
            auto cell = std::find_if(reference.begin(), reference.end(), [&](const BezierCell &c) {
                bool outer = edge == 0   ? c.u0 == reference.front().u0
                             : edge == 1 ? c.u1 == reference.back().u1
                             : edge == 2 ? c.v0 == reference.front().v0
                                         : c.v1 == reference.back().v1;
                return outer && (edge < 2 ? c.v0 : c.u0) == segment.begin &&
                       (edge < 2 ? c.v1 : c.u1) == segment.end;
            });
            Require(cell != reference.end(), "Direct boundary partition differs from 2D reference");
            for (double t : {0., .29, 1.})
            {
                double value = edge == 0   ? s.MinParamU()
                               : edge == 1 ? s.MaxParamU()
                               : edge == 2 ? s.MinParamV()
                                           : s.MaxParamV();
                double low = edge < 2 ? cell->u0 : cell->v0, high = edge < 2 ? cell->u1 : cell->v1;
                I fixed = value == low    ? I(0)
                          : value == high ? I(1)
                                          : (I(value) - I(low)) / (I(high) - I(low));
                auto eval = [&](bool along) {
                    BV p = paired::Restrict(Derivative(cell->p, along), edge < 2, fixed, fixed);
                    return Range(Evaluate(p, edge < 2 ? 0 : t, edge < 2 ? t : 0));
                };
                IV du = eval(true) * (I(1) / (I(cell->u1) - I(cell->u0)));
                IV dv = eval(false) * (I(1) / (I(cell->v1) - I(cell->v0)));
                for (auto pair :
                     {std::make_pair(Range(Evaluate(segment.tangent, t, 0)), edge < 2 ? dv : du),
                      std::make_pair(Range(Evaluate(segment.normal, t, 0)), Cross(du, dv))})
                    for (auto coordinate : {std::make_pair(pair.first.x, pair.second.x),
                                            std::make_pair(pair.first.y, pair.second.y),
                                            std::make_pair(pair.first.z, pair.second.z)})
                        Require(coordinate.first.finite() && coordinate.second.finite() &&
                                    coordinate.first.lo <= coordinate.second.hi &&
                                    coordinate.second.lo <= coordinate.first.hi,
                                "Direct boundary derivative/normal misses independent enclosure");
            }
        }
    for (size_t i = 0; i < cells.size(); ++i)
    {
        const auto &a = cells[i], &b = reference[i];
        Require(a.u0 == b.u0 && a.u1 == b.u1 && a.v0 == b.v0 && a.v1 == b.v1,
                "Source extraction domain changed");
        for (auto pair : {std::make_pair(&a.p.x, &b.p.x), std::make_pair(&a.p.y, &b.p.y),
                          std::make_pair(&a.p.z, &b.p.z)})
        {
            Require(pair.first->c.size() == pair.second->c.size(), "Extraction degree changed");
            for (size_t k = 0; k < pair.first->c.size(); ++k)
            {
                I x = pair.first->c[k], y = pair.second->c[k];
                Require(x.finite() && y.finite() && x.lo <= y.hi && y.lo <= x.hi,
                        "Independent extraction coefficient enclosures do not overlap");
            }
        }
        for (double u : {.0, .31, 1.0})
            for (double v : {.0, .73, 1.0})
            {
                IV value = Range(Evaluate(a.p, u, v));
                V point = Get(s.CalcPoint(a.u0 + (a.u1 - a.u0) * u, a.v0 + (a.v1 - a.v0) * v));
                for (auto pair :
                     {std::make_pair(value.x, point.x), std::make_pair(value.y, point.y),
                      std::make_pair(value.z, point.z)})
                    Require(pair.first.lo - 1e-10 <= pair.second &&
                                pair.second <= pair.first.hi + 1e-10,
                            "SDK evaluation outside extraction enclosure");
            }
    }
}
sggk::BSplineSurfacePtr Plane(unsigned mode)
{
    sggk::Point3DMatrix cp(4, std::vector<sggk::Point3D>(4));
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
        {
            // mode 4: x=(u+u^2)/2, represented exactly in degree three.
            double x =
                mode == 4 ? (double(i) / 6 + double(i * (i ? i - 1 : 0)) / 12) : double(i) / 3;
            if (mode == 5)
                x = double(i) / 3 +
                    .1 * (double(i) / 3 - double(i * (i ? i - 1 : 0)) / 6) * double(j) / 3;
            cp[i][j] = sggk::Point3D(x * (mode >= 4 ? 1 : .001), double(j) / 3, 0);
        }
    if (mode == 1)
        for (unsigned j = 1; j < 3; ++j)
            cp[1][j] = sggk::Point3D(.001 / 3, double(j) / 3, 8e-6);
    if (mode == 2)
        for (unsigned i = 1; i < 3; ++i)
            cp[i][0] = sggk::Point3D(double(i) * .001 / 3, 8e-6, 0);
    if (mode == 3)
        for (unsigned i = 1; i < 3; ++i)
            for (unsigned j = 1; j < 3; ++j)
                cp[i][j] = sggk::Point3D(double(i) * .001 / 3, double(j) / 3, 1e-3);
    return std::make_shared<sggk::BSplineSurface>(3, 3, cp, sggk::RealArray{0, 1},
                                                  sggk::RealArray{0, 1}, sggk::UIntArray{4, 4},
                                                  sggk::UIntArray{4, 4});
}
void PositionArithmetic()
{
    paired::Products products;
    auto sameEnclosure = [&](const BV &a, const BV &b) {
        for (auto pair :
             {std::make_pair(&a.x, &b.x), std::make_pair(&a.y, &b.y), std::make_pair(&a.z, &b.z)})
        {
            Require(pair.first->u == pair.second->u && pair.first->v == pair.second->v,
                    "Weighted product degrees differ");
            for (size_t i = 0; i < pair.first->c.size(); ++i)
                Require(pair.first->c[i].finite() && pair.second->c[i].finite() &&
                            pair.first->c[i].lo <= pair.second->c[i].hi &&
                            pair.second->c[i].lo <= pair.first->c[i].hi,
                        "Weighted product enclosures do not overlap");
        }
    };
    // Independent reference: the shared implementation uses successive elevation.
    for (auto degree : {std::make_pair(0u, 0u), {0u, 3u}, {2u, 0u}, {2u, 3u}, {7u, 7u}, {32u, 1u}})
    {
        BPoly p(degree.first, degree.second);
        for (size_t i = 0; i < p.c.size(); ++i)
            p.c[i] = I(std::sin(double(i)), std::sin(double(i)) + 1e-10);
        unsigned u = p.u + 3, v = p.v + 2;
        auto direct = products.Elevate(p, u, v), reference = p + BPoly(u, v);
        Require(direct.u == u && direct.v == v, "Wrong direct elevation degrees");
        for (size_t i = 0; i < direct.c.size(); ++i)
            Require(direct.c[i].finite() && direct.c[i].lo <= reference.c[i].hi &&
                        reference.c[i].lo <= direct.c[i].hi,
                    "Direct/reference elevation enclosures are disjoint");
        BV a{p, p, p}, b{BPoly(1, 1, I(.3)), BPoly(0, 0, I(-.2)), BPoly(2, 1, I(.7))};
        auto difference = products.Subtract(a, b), old = a - b;
        for (auto pair :
             {std::make_pair(&difference.x, &old.x), std::make_pair(&difference.y, &old.y),
              std::make_pair(&difference.z, &old.z)})
            for (size_t i = 0; i < pair.first->c.size(); ++i)
                Require(pair.first->c[i].lo <= pair.second->c[i].hi &&
                            pair.second->c[i].lo <= pair.first->c[i].hi,
                        "Direct subtraction differs from reference enclosure");
    }
    BPoly quadratic(2, 0);
    quadratic.c = {I(0), I(0), I(1)};
    auto elevated = products.Elevate(quadratic, 4, 0);
    // Exact Bernstein coefficients of t^2 at degree four.
    const double expected[] = {0, 0, 1.0 / 6, .5, 1};
    for (unsigned i = 0; i < 5; ++i)
        Require(elevated.c[i].lo <= expected[i] && elevated.c[i].hi >= expected[i],
                "Known polynomial elevation failed");
    BPoly w(1, 0);
    w.c = {I(1), I(4)};
    BV e{w * I(.5), BPoly(), BPoly()};
    Require(paired::PositionBoxUpper(e, w) > 1 &&
                paired::PositionControlUpper(e, w, products) < .5000000001,
            "Tight rational bound did not resolve a loose box bound");
    Require(!std::isfinite(paired::PositionBoxUpper(e, BPoly(0, 0, I(0)))) &&
                !std::isfinite(paired::PositionControlUpper(e, BPoly(0, 0, I(-1)), products)),
            "Nonpositive weights were accepted");
    // Known chart factorization must survive orientation, chart seams and
    // restriction. Compare against the independent, full tensor product.
    using T = sggk::SurfaceType;
    for (auto type : {T::Plane, T::Cylinder, T::Cone, T::Sphere, T::Torus})
        for (unsigned flags = 0; flags < 8; ++flags)
        {
            AnalyticPatch p;
            p.parameters.type = type;
            p.parameters.radius = .8;
            p.parameters.majorRadius = 2;
            p.parameters.halfAngle = .3;
            p.domain = {-.4, 3.1, .3, 2.7};
            p.mapping.reversedU = (flags & 1) != 0;
            p.mapping.reversedV = (flags & 2) != 0;
            p.mapping.swappedUV = (flags & 4) != 0;
            auto chart = paired::BuildChart(p);
            for (auto cell : chart.cells)
            {
                Require(cell.separableWeight && !cell.weightFactors[0].v &&
                            !cell.weightFactors[1].u,
                        "Wrong oriented factor axes");
                cell.surface.w = paired::Restrict(
                    paired::Restrict(cell.surface.w, true, I(.17), I(.83)), false, I(.11), I(.71));
                cell.weightFactors[0] =
                    paired::Restrict(cell.weightFactors[0], true, I(.17), I(.83));
                cell.weightFactors[1] =
                    paired::Restrict(cell.weightFactors[1], false, I(.11), I(.71));
                BV a{BPoly(2, 3), BPoly(2, 3), BPoly(2, 3)};
                for (size_t i = 0; i < a.x.c.size(); ++i)
                {
                    a.x.c[i] = I(std::sin(double(i)));
                    a.y.c[i] = I(double(i) - 2);
                    a.z.c[i] = I(.1);
                }
                auto reference = products.Mul(a, cell.surface.w);
                sameEnclosure(products.MulWeight(a, cell), reference);
                cell.separableWeight = false;
                sameEnclosure(products.MulWeight(a, cell), reference);
            }
        }
}
void Unit()
{
    PositionArithmetic();
    AnalyticPatch patch;
    patch.parameters.type = sggk::SurfaceType::Plane;
    patch.domain = {0, .001, 0, 1};
    MaterializePatch(patch);
    auto exact = ValidatePairedConversion(*Plane(0), patch);
    Timing(exact);
    Require(exact.Passed(), "Exact plane: " + exact.certificate.message);
    Require(exact.fastBoundPasses == exact.surfaceNodes && !exact.tightBoundChecks,
            "Exact plane did not use the fast position bound");
    {
        auto p = patch;
        p.domain = {0, 1, 0, 1};
        auto s = Plane(4);
        // Error control vectors point in different directions. The box bound is
        // above tolerance, while every control vector is strictly below it.
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                s->SetControlPoint(sggk::Point3D(double(i) / 3 + 8e-6 * (1 - double(i) / 3),
                                                 double(j) / 3 + 8e-6 * double(i) / 3, 0),
                                   i, j);
        PairedValidationOptions o;
        o.enablePiecewiseMapping = false;
        auto fallback = ValidatePairedConversion(*s, p, o);
        Timing(fallback);
        Require(fallback.Passed() && fallback.tightBoundChecks > 0 && !fallback.fastBoundPasses,
                "Loose fast bound was not resolved by the tight fallback");
        const double bulge[] = {0, 2e-5, -2e-5, 0};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                s->SetControlPoint(sggk::Point3D(double(i) / 3, double(j) / 3, bulge[i]), i, j);
        auto subdivided = ValidatePairedConversion(*s, p, o);
        Timing(subdivided);
        Require(subdivided.Passed() && subdivided.surfaceNodes > subdivided.commonCells &&
                    subdivided.tightBoundChecks > 0 &&
                    subdivided.certificate.timing.subdivision.elapsedMs > 0,
                "Loose bounds incorrectly rejected a surface requiring subdivision");
        for (double offset : {9.9e-6, 1e-5, 1.01e-5})
        {
            auto shiftedSource = Plane(0);
            for (unsigned i = 0; i < 4; ++i)
                for (unsigned j = 0; j < 4; ++j)
                {
                    const auto &point = shiftedSource->ControlPoints()[i][j];
                    shiftedSource->SetControlPoint(sggk::Point3D(point.X(), point.Y(), offset), i,
                                                   j);
                }
            auto r = ValidatePairedConversion(*shiftedSource, patch, o);
            Timing(r);
            Require(r.Passed() == (offset < 1e-5), "Strict position tolerance boundary failed");
        }
    }
    for (unsigned i = 1; i < 4; ++i)
    {
        PairedValidationOptions o;
        o.enablePiecewiseMapping = false;
        auto r = ValidatePairedConversion(*Plane(i), patch, o);
        Timing(r);
        Require(!r.Passed(), "Missed distance/tangent/normal violation " + std::to_string(i));
        Require(i != 1 || r.certificate.boundaryNormal.state == CheckState::Exceeded,
                "Normal-only violation not detected");
        Require(i != 2 || r.certificate.boundaryTangent.state == CheckState::Exceeded,
                "Tangent-only violation not detected");
        Require(i != 3 || r.certificate.geometry.state == CheckState::Exceeded,
                "Interior bulge not detected");
    }
    auto warped = Plane(4);
    patch.domain = {0, 1, 0, 1};
    MaterializePatch(patch);
    PairedValidationOptions affine;
    affine.enablePiecewiseMapping = false;
    auto noMap = ValidatePairedConversion(*warped, patch, affine);
    Require(!noMap.Passed(), "Affine correspondence should fail on warped plane");
    auto mapped = ValidatePairedConversion(*warped, patch);
    Timing(mapped);
    std::cout << "warped " << PairedValidationStatusName(mapped.status)
              << " segments=" << mapped.uMap.source.size() - 1 << " " << mapped.certificate.message
              << std::endl;
    Require(mapped.Passed() && mapped.mappingAttempts == 2 && mapped.uMap.source.size() > 2,
            "Piecewise mapping failed on warped plane");
    {
        paired::AnalysisContext shared(*warped);
        shared.classification = BuildClassificationSamples(*warped, nullptr, &shared.samples);
        RecoveredCandidate candidate;
        candidate.patch =
            RecoverApproxPatch(*shared.classification, patch.parameters, &candidate.trace);
        candidate.traceParameters = candidate.patch.parameters;
        candidate.hasTrace = true;
        candidate.sourceStamp = shared.classification->stamp;
        auto a = paired::ValidatePrepared(shared, candidate);
        auto b = ValidatePairedConversion(*warped, candidate.patch);
        Require(a.Passed() && b.Passed() && a.recoveryTraceHits == 1 && !a.recoveryTraceBuilds &&
                    a.sourceExtractions == 1 && a.sourceBoundaryBuilds == 1 &&
                    a.uMap.source == b.uMap.source && a.uMap.target == b.uMap.target &&
                    a.certificate.geometry.upperBound == b.certificate.geometry.upperBound &&
                    a.certificate.boundaryNormal.upperBound ==
                        b.certificate.boundaryNormal.upperBound,
                "Recovery trace reuse changed the certificate");
        auto warm = paired::ValidatePrepared(shared, candidate);
        Require(warm.Passed() && warm.sourceCacheHits == 1 && !warm.sourceBoundaryBuilds &&
                    !warm.sourceExtractions && warm.sourceBoundaryMs == 0 &&
                    warm.sourceCenterMs == 0,
                "Shared source was rebuilt on revalidation");
        auto wrongStamp = candidate;
        wrongStamp.sourceStamp = std::make_shared<const int>(0);
        auto rebuilt = paired::ValidatePrepared(shared, wrongStamp);
        Require(rebuilt.Passed() && rebuilt.recoveryTraceBuilds == 1 && !rebuilt.recoveryTraceHits,
                "Foreign source trace was reused");
        candidate.patch.parameters.origin = sggk::Point3D(0, 0, .01);
        MaterializePatch(candidate.patch);
        auto changed = paired::ValidatePrepared(shared, candidate);
        Require(!changed.Passed() && changed.recoveryTraceBuilds == 1 &&
                    !changed.recoveryTraceHits && changed.sourceCacheHits == 1,
                "Stale mother-parameter trace was reused");
        Timing(a);
        Timing(warm);
        Timing(rebuilt);
        Timing(changed);
    }
    {
        // Identical varying directions require subdivision of their dot-product
        // hull; endpoint/midpoint samples alone must not be the certificate.
        paired::Boundaries edges;
        BV tangent{BPoly(1, 0), BPoly(0, 0, I(.01)), BPoly(0, 0)};
        tangent.x.c = {I(-.9), I(1.1)};
        BV normal{BPoly(0, 0), BPoly(0, 0), BPoly(0, 0, I(1))};
        for (auto &edge : edges)
            edge.push_back({0, 1, tangent, normal});
        PairedParameterMap map{{0, 1}, {0, 1}};
        PairedValidationReport report;
        paired::Products products;
        auto r = paired::ValidateBoundaries(edges, edges, map, map, {}, report, products);
        Require(r.passed && report.boundarySegments == 4 && report.boundaryNodes > 4 &&
                    report.edgeSubdivisionMs > 0,
                "One-dimensional angle subdivision failed");
        auto wrong = edges;
        wrong[3][0].normal = normal * I(-1);
        report = {};
        Require(!paired::ValidateBoundaries(edges, wrong, map, map, {}, report, products).passed,
                "Reversed boundary normal was accepted");
        auto constant = edges;
        BV direction{BPoly(0, 0, I(1)), BPoly(0, 0), BPoly(0, 0)};
        for (auto &edge : constant)
            edge[0].tangent = direction;
        for (double degrees : {.099, .101})
        {
            auto tilted = constant;
            double angle = degrees * Pi / 180;
            tilted[2][0].tangent = {BPoly(0, 0, I(std::cos(angle))),
                                    BPoly(0, 0, I(std::sin(angle))), BPoly(0, 0)};
            report = {};
            Require(paired::ValidateBoundaries(constant, tilted, map, map, {}, report, products)
                            .passed == (degrees < .1),
                    "Boundary angle threshold was weakened");
        }
        PairedValidationOptions limit;
        limit.validation.maxValidationPieces = 3;
        bool exhausted = false;
        report = {};
        try
        {
            (void)paired::ValidateBoundaries(constant, constant, map, map, limit, report, products);
        }
        catch (const paired::Budget &)
        {
            exhausted = true;
        }
        Require(exhausted, "Boundary checks did not share the node budget");
    }
    for (unsigned flags = 1; flags < 8; ++flags)
    {
        auto oriented = Plane(4);
        if (flags & 1)
            oriented->UReverse();
        if (flags & 2)
            oriented->VReverse();
        if (flags & 4)
            oriented->ExchangeUV();
        oriented = NormalizeDatasetSurface(*oriented);
        auto p = RecoverApproxPatch(*oriented, patch.parameters);
        auto r = ValidatePairedConversion(*oriented, p);
        Timing(r);
        Require(r.Passed(), "Reversed/swapped piecewise map failed");
    }
    auto differentEdges = Plane(5);
    Require(ValidateConversion(*differentEdges, patch).Passed(),
            "Different-speed plane is not a valid geometric fixture");
    auto conservative = ValidatePairedConversion(*differentEdges, patch);
    Timing(conservative);
    Require(!conservative.Passed() && conservative.status == PairedValidationStatus::NotCertified,
            "Incompatible separable map was incorrectly certified");
    // A second call must use changed geometry/domain, not a stale successful map.
    auto shifted = patch;
    shifted.parameters.origin = sggk::Point3D(0, 0, 1e-3);
    MaterializePatch(shifted);
    Require(!ValidatePairedConversion(*warped, shifted, affine).Passed(),
            "Stale candidate geometry accepted");
    auto smaller = patch;
    smaller.domain[1] = .9;
    MaterializePatch(smaller);
    Require(!ValidatePairedConversion(*warped, smaller, affine).Passed(),
            "Stale candidate domain accepted");
    auto flat = Plane(0);
    auto optimized = RecoverApproxPatch(*flat, patch.parameters);
    optimized.parameters.origin = sggk::Point3D(0, 0, 1e-3);
    MaterializePatch(optimized);
    auto samples = SampleSurface(*flat, {});
    ConversionStatistics optimization;
    PairedValidationContext optimizationContext;
    Require(!ValidatePairedConversion(*flat, optimized, optimizationContext, affine).Passed(),
            "Displaced candidate unexpectedly passed before refinement");
    RefineCandidate(*flat, optimized, samples, {}, optimization);
    auto afterOptimization = ValidatePairedConversion(*flat, optimized, optimizationContext);
    Timing(afterOptimization);
    Require(optimization.optimizationIterations > 0 && afterOptimization.Passed() &&
                afterOptimization.sourceCacheHits == 1 && !afterOptimization.sourceExtractions,
            "Independent verification of an actually optimized candidate failed");
    PairedValidationContext context;
    auto cachedSource = Plane(0);
    auto cachedPatch = RecoverApproxPatch(*cachedSource, patch.parameters);
    auto cold = ValidatePairedConversion(*cachedSource, cachedPatch, context);
    auto warm = ValidatePairedConversion(*cachedSource, cachedPatch, context);
    Timing(cold);
    Timing(warm);
    Require(cold.Passed() && warm.Passed() && cold.sourceExtractions == 1 &&
                !cold.sourceCacheHits && warm.sourceCacheHits == 1 && !warm.sourceExtractions &&
                !warm.uExtractionMatrices && !warm.vExtractionMatrices &&
                warm.certificate.geometry.upperBound == cold.certificate.geometry.upperBound,
            "Prepared source was not reused exactly");
    auto displaced = cachedPatch;
    displaced.parameters.origin = sggk::Point3D(0, 0, .001);
    auto changedCandidate = ValidatePairedConversion(*cachedSource, displaced, context, affine);
    Require(!changedCandidate.Passed() && changedCandidate.sourceCacheHits == 1,
            "Source cache retained a stale candidate certificate");
    auto narrow = cachedPatch;
    narrow.domain[1] *= .8;
    Require(!ValidatePairedConversion(*cachedSource, narrow, context, affine).Passed(),
            "Cached validation retained a stale candidate domain");
    const auto oldPoint = cachedSource->ControlPoints()[1][1];
    cachedSource->SetControlPoint(sggk::Point3D(oldPoint.X(), oldPoint.Y(), .001), 1, 1);
    auto changedSource = ValidatePairedConversion(*cachedSource, cachedPatch, context, affine);
    Require(!changedSource.Passed() && changedSource.sourceExtractions == 1 &&
                !changedSource.sourceCacheHits,
            "In-place source mutation reused stale cells");
    cachedSource->SetControlPoint(oldPoint, 1, 1);
    cachedSource->InsertUKnots(.31);
    auto inserted = ValidatePairedConversion(*cachedSource, cachedPatch, context);
    Require(inserted.Passed() && inserted.sourceExtractions == 1 && !inserted.sourceCacheHits &&
                inserted.uExtractionMatrices == 2 && inserted.vExtractionMatrices == 1,
            "Changed knots did not invalidate source preparation");
    PairedValidationOptions smallBudget;
    smallBudget.maxCommonCells = 1;
    Require(ValidatePairedConversion(*cachedSource, cachedPatch, context, smallBudget).status ==
                PairedValidationStatus::BudgetExceeded,
            "Cached source bypassed a reduced cell budget");
    auto movedContext = std::move(context);
    Require(ValidatePairedConversion(*cachedSource, cachedPatch, movedContext).sourceCacheHits == 1,
            "Moved context lost valid source preparation");
    Require(ValidatePairedConversion(*cachedSource, cachedPatch, context).sourceExtractions == 1,
            "Moved-from context was not reusable");
    CheckExtraction(*cachedSource);
    sggk::Point3DMatrix net(6, std::vector<sggk::Point3D>(7));
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 7; ++j)
            net[i][j] = sggk::Point3D(.3 * i, .2 * j, std::sin(double(i + 2 * j)));
    auto asymmetric = std::make_shared<sggk::BSplineSurface>(
        2, 3, net, sggk::RealArray{-.3, .1, .4, 1.7}, sggk::RealArray{2, 2.2, 3.3, 5},
        sggk::UIntArray{3, 2, 1, 3}, sggk::UIntArray{4, 2, 1, 4});
    CheckExtraction(*asymmetric);
    auto unclamped = std::make_shared<sggk::BSplineSurface>(
        2, 3, net, sggk::RealArray{0, 1, 2, 3, 4, 5, 6, 7, 8},
        sggk::RealArray{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, sggk::UIntArray(9, 1),
        sggk::UIntArray(11, 1));
    CheckExtraction(*unclamped);
    asymmetric->ExchangeUV();
    asymmetric->UReverse();
    asymmetric->VReverse();
    CheckExtraction(*asymmetric);
    PairedValidationOptions budget;
    budget.validation.maxValidationPieces = 1;
    auto limited = ValidatePairedConversion(
        *Plane(0), RecoverApproxPatch(*Plane(0), patch.parameters), budget);
    Timing(limited);
    Require(limited.status == PairedValidationStatus::BudgetExceeded && !limited.Passed(),
            "Node budget not enforced");
    budget = {};
    budget.maxMappingSegments = 2;
    Require(ValidatePairedConversion(*warped, patch, budget).status ==
                PairedValidationStatus::BudgetExceeded,
            "Mapping budget not enforced");
    auto wrong = patch;
    wrong.normalSign = -1;
    Require(!ValidatePairedConversion(*warped, wrong, affine).Passed(),
            "Opposite normals accepted");
    PairedValidationOptions bad;
    bad.validation.distanceTolerance = 0;
    Require(ValidatePairedConversion(*warped, patch, bad).status ==
                PairedValidationStatus::InvalidInput,
            "Invalid tolerance accepted");
    int rounding = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    auto roundingResult = ValidatePairedConversion(*warped, patch);
    std::fesetround(rounding);
    Require(!roundingResult.Passed(), "Unsupported rounding mode accepted");
    // Verify each analytic rational chart against direct analytic evaluation,
    // including multi-span arcs, swaps/reversals and both angular directions.
    for (auto type : {sggk::SurfaceType::Plane, sggk::SurfaceType::Cylinder,
                      sggk::SurfaceType::Cone, sggk::SurfaceType::Sphere, sggk::SurfaceType::Torus})
        for (unsigned flags = 0; flags < 8; ++flags)
        {
            AnalyticPatch p;
            p.parameters.type = type;
            p.parameters.radius = 2;
            p.parameters.majorRadius = 4;
            p.parameters.halfAngle = .4;
            p.parameters.origin = sggk::Point3D(1e6, -2e6, 3e6);
            p.parameters.axis = Direction({.3, .7, .9});
            p.parameters.xDirection = Direction({.8, -.4, .2});
            p.domain = {6.1, 8.3, .4, 2.1};
            p.mapping.swappedUV = bool(flags & 1);
            p.mapping.reversedU = bool(flags & 2);
            p.mapping.reversedV = bool(flags & 4);
            auto chart = paired::BuildChart(p);
            auto placed = chart;
            V sourceOrigin{1e6 + 4, -2e6 + 3, 3e6 - 2};
            paired::PlaceChart(placed, p.parameters, sourceOrigin);
            auto frame = FrameOf(p.parameters);
            for (size_t cellIndex = 0; cellIndex < chart.cells.size(); ++cellIndex)
            {
                const auto &c = chart.cells[cellIndex];
                double u = .31, v = .67, ou = c.surface.u0 + (c.surface.u1 - c.surface.u0) * u,
                       ov = c.surface.v0 + (c.surface.v1 - c.surface.v0) * v;
                double coordinates[2] = {p.mapping.reversedU ? 1 - ou : ou,
                                         p.mapping.reversedV ? 1 - ov : ov};
                if (p.mapping.swappedUV)
                    std::swap(coordinates[0], coordinates[1]);
                double angles[2]{};
                for (unsigned d = 0; d < 2; ++d)
                {
                    double scaled = coordinates[d] * (chart.angles[d].size() - 1);
                    unsigned k = static_cast<unsigned>(scaled);
                    double t = scaled - k, a = chart.angles[d][k], b = chart.angles[d][k + 1];
                    angles[d] = chart.angular[d]
                                    ? a + 2 * std::atan2(t * std::sin((b - a) * .5),
                                                         1 - t + t * std::cos((b - a) * .5))
                                    : a + t * (b - a);
                    Require(std::abs(paired::RationalCoordinate(chart, d, angles[d]) -
                                     coordinates[d]) < 1e-13,
                            "Angular inverse mismatch");
                }
                IV num = Range(Evaluate(c.surface.p, u, v));
                I w = Evaluate(c.surface.w, u, v).range();
                V expected = CanonicalPoint(p.parameters, angles[0], angles[1]);
                for (auto pair :
                     {std::make_pair(num.x / w, expected.x), std::make_pair(num.y / w, expected.y),
                      std::make_pair(num.z / w, expected.z)})
                    Require(pair.first.lo - 1e-12 <= pair.second &&
                                pair.second <= pair.first.hi + 1e-12,
                            "Rational chart evaluation mismatch");
                IV transformed = Range(Evaluate(placed.cells[cellIndex].surface.p, u, v));
                IV normal = Range(Evaluate(placed.cells[cellIndex].normal, u, v));
                V point = frame.global(expected) - sourceOrigin;
                V expectedNormal =
                    frame.vector(CanonicalNormal(p.parameters, angles[0], angles[1]));
                for (auto pair : {std::make_pair(transformed.x / w, point.x),
                                  std::make_pair(transformed.y / w, point.y),
                                  std::make_pair(transformed.z / w, point.z),
                                  std::make_pair(normal.x / w, expectedNormal.x),
                                  std::make_pair(normal.y / w, expectedNormal.y),
                                  std::make_pair(normal.z / w, expectedNormal.z)})
                    Require(pair.first.lo - 2e-9 <= pair.second &&
                                pair.second <= pair.first.hi + 2e-9,
                            "Fixed source frame changes analytic position/normal");
            }
        }
    {
        auto translatedPlane = Plane(0);
        V delta{1e6, -2e6, 3e6};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                translatedPlane->SetControlPoint(Point(Get(translatedPlane->ControlPoints()[i][j]) + delta), i, j);
        AnalyticPatch p;
        p.parameters.type = sggk::SurfaceType::Plane;
        p.parameters.origin = Point(delta);
        p.domain = {0, .001, 0, 1};
        MaterializePatch(p);
        auto result = ValidatePairedConversion(*translatedPlane, p);
        Require(result.Passed(), "Large translated source lost its certificate");
        Timing(result);
    }
    // Exercise the entire common-partition certificate across rational chart
    // seams, including a periodic seam and a swapped/reversed source axis.
    AnalyticPatch cylinder;
    cylinder.parameters.type = sggk::SurfaceType::Cylinder;
    cylinder.parameters.radius = 2;
    cylinder.domain = {6.1, 8.3, -.4, 1.2};
    MaterializePatch(cylinder);
    auto rational = cylinder.surface->ToBSpline(
        sggk::UVRange(sggk::Interval(6.1, 8.3), sggk::Interval(-.4, 1.2)));
    Require(bool(rational), "Multi-span reference creation failed");
    for (unsigned flags : {0u, 7u})
    {
        auto polynomial = RefitPolynomial(*rational, 12);
        if (flags & 1)
            polynomial->UReverse();
        if (flags & 2)
            polynomial->VReverse();
        if (flags & 4)
            polynomial->ExchangeUV();
        polynomial = NormalizeDatasetSurface(*polynomial);
        CheckExtraction(*polynomial);
        auto recovered = RecoverApproxPatch(*polynomial, cylinder.parameters);
        auto result = ValidatePairedConversion(*polynomial, recovered);
        Timing(result);
        Require(result.Passed() && result.targetPatches > 1,
                "Multi-span cylinder certificate failed: " + result.certificate.message);
    }
    std::cout << "Paired unit assertions: " << checks << std::endl;
}
void SpotCheck(const sggk::BSplineSurface &s, const AnalyticPatch &p,
               const PairedValidationReport &r)
{
    auto chart = paired::BuildChart(p);
    auto frame = FrameOf(p.parameters);
    auto map = [](const PairedParameterMap &m, double t) {
        if (t == m.source.back())
            return m.target.back();
        size_t i = static_cast<size_t>(std::upper_bound(m.source.begin(), m.source.end(), t) -
                                       m.source.begin() - 1);
        return m.target[i] + (m.target[i + 1] - m.target[i]) * (t - m.source[i]) /
                                 (m.source[i + 1] - m.source[i]);
    };
    // Independent double evaluation checks the certificate, never establishes it.
    SurfaceSampleContext cache(s);
    for (unsigned i = 0; i <= 8; ++i)
        for (unsigned j = 0; j <= 8; ++j)
        {
            double u =
                i == 8 ? s.MaxParamU() : s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 8;
            double v =
                j == 8 ? s.MaxParamV() : s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * j / 8;
            double q[2] = {map(r.uMap, u), map(r.vMap, v)};
            if (p.mapping.reversedU)
                q[0] = 1 - q[0];
            if (p.mapping.reversedV)
                q[1] = 1 - q[1];
            if (p.mapping.swappedUV)
                std::swap(q[0], q[1]);
            double angle[2];
            for (unsigned d = 0; d < 2; ++d)
            {
                unsigned n = static_cast<unsigned>(chart.angles[d].size() - 1);
                unsigned k = std::min(n - 1, static_cast<unsigned>(std::max(0.0, q[d]) * n));
                double t = q[d] * n - k, a = chart.angles[d][k], b = chart.angles[d][k + 1];
                angle[d] = chart.angular[d] ? a + 2 * std::atan2(t * std::sin((b - a) / 2),
                                                                 1 - t + t * std::cos((b - a) / 2))
                                            : a + t * (b - a);
            }
            V target = frame.global(CanonicalPoint(p.parameters, angle[0], angle[1]));
            const auto &sample = cache.At(u, v);
            Require(Norm(sample.point - target) <= r.certificate.geometry.upperBound + 1e-10,
                    "Sample exceeds paired position certificate");
            if (i == 0 || i == 8 || j == 0 || j == 8)
            {
                V n = frame.vector(CanonicalNormal(p.parameters, angle[0], angle[1])) *
                      double(p.normalSign);
                Require(Angle(sample.normal, n) * 180 / Pi <=
                            r.certificate.boundaryNormal.upperBound + 1e-7,
                        "Sample exceeds paired normal certificate");
            }
        }
}
int Dataset(const std::filesystem::path &dir, const std::filesystem::path &file, bool sphere,
            bool compareLegacy = true)
{
    std::ofstream csv(file);
    Require(bool(csv), "Cannot open paired CSV");
    csv << std::setprecision(17);
    csv << "id,type,rational,collapsed,controls,classify_ms,candidates,old_status,new_status,old_"
           "ms,new_ms,speedup,geometry_bound,tangent_deg,normal_deg,attempts,u_segments,v_segments,"
           "common_cells,surface_nodes,boundary_nodes,prepare_ms,surface_ms,subdivide_ms,boundary_"
           "ms,report_ms,nurbs_ms,mapping_ms,extraction_ms,source_cache_ms,local_transform_ms,"
           "source_extractions,cache_hits,u_matrices,v_matrices,compose_ms,product_ms,difference_"
           "ms,"
           "fast_bound_ms,tight_bound_ms,fast_passes,tight_checks,source_center_ms,source_boundary_"
           "ms,"
           "target_transform_ms,target_boundary_ms,edge_compose_ms,edge_algebra_ms,edge_bounds_ms,"
           "edge_subdivision_ms,source_boundary_builds,boundary_segments,trace_hits,trace_builds,"
           "timing_schema,message\n";
    unsigned passed = 0, notCertified = 0, skipped = 0;
    bool warmed[5]{};
    for (unsigned id = sphere ? 301u : 0u; id < (sphere ? 400u : 500u); id += sphere ? 2u : 1u)
    {
        auto c = ReadDatasetCase(dir / CaseFilename(id));
        CheckDataset500Contract(c, id);
        Require(c.id == id, "Fixture identity mismatch");
        auto inputBound = BoundPair(*c.source, *c.input, c.refitted ? 1e-6 : 1e-9);
        Require(inputBound.passed, "Fixture source certificate failed");
        unsigned controls = static_cast<unsigned>(c.input->ControlPoints().size() *
                                                  c.input->ControlPoints()[0].size());
        std::string oldStatus = "NotRun", newStatus = "NotRun", message;
        ConversionValidationReport old;
        PairedValidationReport newer;
        double classifyMs = 0, oldMs = 0, newMs = 0;
        size_t candidates = 0;
        if (c.input->IsRational())
        {
            // New module intentionally supports polynomial input only.
            AnalyticPatch p;
            p.parameters = c.truth;
            p.domain = c.analyticDomain;
            auto r = ValidatePairedConversion(*c.input, p);
            Require(r.status == PairedValidationStatus::Unsupported,
                    "Rational input was unexpectedly accepted");
            newStatus = "Unsupported";
            message = r.certificate.message;
            ++skipped;
        }
        else
        {
            try
            {
                auto start = std::chrono::steady_clock::now();
                paired::AnalysisContext analysis(*c.input);
                std::vector<RecoveredCandidate> recovered;
                analysis.classification.emplace();
                auto list =
                    ClassifyApproximate(*c.input, ConversionOptions{}, nullptr, true,
                                        &analysis.samples, &recovered, &*analysis.classification);
                classifyMs = std::chrono::duration<double, std::milli>(
                                 std::chrono::steady_clock::now() - start)
                                 .count();
                candidates = list.size();
                if (list.empty())
                {
                    oldStatus = newStatus = "NoCandidate";
                    ++skipped;
                }
                else
                {
                    // First actual recovered candidate; never substitute fixture truth.
                    const auto &p = list.front();
                    if (sphere && !c.collapsedMask)
                        Require(p.parameters.type == sggk::SurfaceType::Sphere,
                                "Regular sphere must precede near-sphere torus candidates");
                    SurfaceSampleContext referenceSamples(*c.input, 50000);
                    auto ordinary = ClassifyApproximate(*c.input, ConversionOptions{}, nullptr,
                                                        true, &referenceSamples);
                    Require(ordinary.size() == list.size() && recovered.size() == list.size(),
                            "Prepared classification changed candidate count");
                    for (size_t i = 0; i < list.size(); ++i)
                        Require(SameMotherParameters(ordinary[i].parameters, list[i].parameters) &&
                                    ordinary[i].domain == list[i].domain &&
                                    ordinary[i].mapping.swappedUV == list[i].mapping.swappedUV &&
                                    ordinary[i].mapping.reversedU == list[i].mapping.reversedU &&
                                    ordinary[i].mapping.reversedV == list[i].mapping.reversedV,
                                "Prepared classification changed candidate order/geometry");
                    if (!warmed[c.type] && !c.collapsedMask)
                    {
                        if (compareLegacy) (void)ValidateConversion(*c.input, p);
                        (void)ValidatePairedConversion(*c.input, p);
                        warmed[c.type] = true;
                    }
                    auto oldRun = [&] {
                        if (!compareLegacy) return;
                        auto t = std::chrono::steady_clock::now();
                        old = ValidateConversion(*c.input, p);
                        oldMs = std::chrono::duration<double, std::milli>(
                                    std::chrono::steady_clock::now() - t)
                                    .count();
                    };
                    auto newRun = [&] {
                        auto t = std::chrono::steady_clock::now();
                        newer = paired::ValidatePrepared(analysis, recovered.front());
                        newMs = std::chrono::duration<double, std::milli>(
                                    std::chrono::steady_clock::now() - t)
                                    .count();
                    };
                    if (c.pair % 2)
                    {
                        newRun();
                        oldRun();
                    }
                    else
                    {
                        oldRun();
                        newRun();
                    }
                    Timing(newer);
                    oldStatus = compareLegacy ? (old.Passed() ? "Passed" : "NotCertified") : "NotRun";
                    newStatus = PairedValidationStatusName(newer.status);
                    message = newer.certificate.message;
                    if (newer.Passed())
                    {
                        Require(!c.collapsedMask, "Degenerate input unexpectedly certified");
                        Require(newer.certificate.geometry.upperBound < 1e-5 &&
                                    newer.certificate.boundaryTangent.upperBound < .1 &&
                                    newer.certificate.boundaryNormal.upperBound < .1,
                                "Tolerance violated");
                        SpotCheck(*c.input, p, newer);
                        ++passed;
                    }
                    else
                        ++notCertified;
                }
            }
            catch (const Failure &e)
            {
                oldStatus = newStatus = "ClassificationFailed";
                message = e.what();
                ++skipped;
            }
        }
        const auto &r = newer.certificate;
        const auto &t = r.timing;
        csv << id << ',' << TypeName(c.truth.type) << ',' << c.input->IsRational() << ','
            << c.collapsedMask << ',' << controls << ',' << classifyMs << ',' << candidates << ','
            << oldStatus << ',' << newStatus << ',' << oldMs << ',' << newMs << ','
            << (newMs ? oldMs / newMs : 0) << ',' << r.geometry.upperBound << ','
            << r.boundaryTangent.upperBound << ',' << r.boundaryNormal.upperBound << ','
            << newer.mappingAttempts << ','
            << (newer.uMap.source.empty() ? 0 : newer.uMap.source.size() - 1) << ','
            << (newer.vMap.source.empty() ? 0 : newer.vMap.source.size() - 1) << ','
            << newer.commonCells << ',' << newer.surfaceNodes << ',' << newer.boundaryNodes << ','
            << t.preparation.elapsedMs << ',' << t.surfaceBounds.elapsedMs << ','
            << t.subdivision.elapsedMs << ',' << t.boundaryChecks.elapsedMs << ','
            << t.reporting.elapsedMs << ',' << newer.nurbsMs << ',' << newer.mappingMs << ','
            << newer.sourceExtractionMs << ',' << newer.sourceCacheMs << ','
            << newer.localTransformMs << ',' << newer.sourceExtractions << ','
            << newer.sourceCacheHits << ',' << newer.uExtractionMatrices << ','
            << newer.vExtractionMatrices << ',' << newer.surfaceComposeMs << ','
            << newer.surfaceProductMs << ',' << newer.surfaceDifferenceMs << ','
            << newer.surfaceFastBoundMs << ',' << newer.surfaceTightBoundMs << ','
            << newer.fastBoundPasses << ',' << newer.tightBoundChecks << ',' << newer.sourceCenterMs
            << ',' << newer.sourceBoundaryMs << ',' << newer.targetTransformMs << ','
            << newer.targetBoundaryMs << ',' << newer.edgeComposeMs << ',' << newer.edgeAlgebraMs
            << ',' << newer.edgeBoundsMs << ',' << newer.edgeSubdivisionMs << ','
            << newer.sourceBoundaryBuilds << ',' << newer.boundarySegments << ','
            << newer.recoveryTraceHits << ',' << newer.recoveryTraceBuilds << ','
            << newer.preparationSchema << ',' << Quote(message) << '\n';
        csv.flush();
        std::cout << id << " old=" << oldStatus << " new=" << newStatus << " ms=" << newMs
                  << std::endl;
    }
    std::cout << "Paired comparison passed=" << passed << " not_certified=" << notCertified
              << " skipped=" << skipped << std::endl;
    // Every regular polynomial fixture has an independently certified source.
    // Guard coverage as well as the correctness of individual successes.
    Require(passed == (sphere ? 40u : 280u), "Regular polynomial dataset coverage regressed");
    return 0;
}
} // namespace
int main(int argc, char **argv)
{
    sggk::init();
    int code = 0;
    try
    {
        Require(argc == 2 || argc == 4, "Expected --unit or --sphere/--dataset[-new-only] directory csv");
        if (std::string(argv[1]) == "--unit")
        {
            Require(argc == 2, "Unexpected unit arguments");
            Unit();
        }
        else
        {
            Require(argc == 4 &&
                        (std::string(argv[1]) == "--sphere" || std::string(argv[1]) == "--dataset" ||
                         std::string(argv[1]) == "--sphere-new-only" ||
                         std::string(argv[1]) == "--dataset-new-only"),
                    "Invalid comparison arguments");
            const std::string mode = argv[1];
            code = Dataset(argv[2], argv[3], mode.rfind("--sphere", 0) == 0,
                           mode.find("-new-only") == std::string::npos);
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << std::endl;
        code = 1;
    }
    sggk::fini();
    return code;
}
