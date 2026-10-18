/**
 * @file ExtremeDataset.cpp
 * @brief 等参极端固定输入的构造与测量。
 * 覆盖尺度、短边和退化边，独立计算参考几何；测量弧长与严格来源证书分别记录。
 */
#include "ExtremeDataset.h"
#include "ApproximateInternal.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
ExtremeSpec Base(unsigned type, double rho, double theta)
{
    ExtremeSpec s;
    s.type = type;
    s.rho = rho;
    s.theta = theta;
    s.truth.type = DatasetType(type);
    s.truth.radius = rho;
    s.truth.halfAngle = Pi / 6;
    s.domain = {-.5 * theta, .5 * theta, 0, 1};
    if (type == 2)
    {
        double z = rho / std::tan(s.truth.halfAngle);
        s.domain[2] = z;
        s.domain[3] = z + 1;
    }
    if (type == 3)
    {
        s.truth.radius = rho / std::cos(.01);
        s.domain[2] = Pi / 2 - .01;
        s.domain[3] = Pi / 2 + .01;
    }
    if (type == 4)
    {
        s.truth.radius = .1 * rho;
        s.truth.majorRadius = rho - s.truth.radius * std::cos(.01);
        s.domain[2] = -.01;
        s.domain[3] = .01;
    }
    return s;
}
std::array<double, 4> EdgeLengths(const sggk::BSplineSurface &s, unsigned n)
{
    // 逐非零结点区间采用复合 Simpson 积分；这是弧长测量，不是严格弧长包络，同时报告 n 与 2n 两种分辨率。
    std::array<double, 4> out{};
    for (unsigned e = 0; e < 4; ++e)
    {
        const auto &knots = e < 2 ? s.KnotsV() : s.KnotsU();
        for (size_t k = 1; k < knots.size(); ++k)
        {
            double a = knots[k - 1], b = knots[k], sum = 0;
            for (unsigned j = 0; j <= n; ++j)
            {
                double t = a + (b - a) * j / n, u = e < 2 ? (e == 0 ? 0 : 1) : t,
                       v = e < 2 ? t : (e == 2 ? 0 : 1);
                auto side = j == n ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
                auto fixed = e % 2 ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
                sggk::Point3D p;
                auto d = s.CalcPartialD1(u, v, e < 2 ? fixed : side, e < 2 ? side : fixed, p);
                sum += (j == 0 || j == n ? 1 : (j % 2 ? 4 : 2)) * Norm(Get(d[e < 2 ? 1 : 0]));
            }
            out[e] += sum * (b - a) / (3 * n);
        }
    }
    return out;
}
void OrientSurface(sggk::BSplineSurface &s, unsigned flags)
{
    if (flags & 1)
        s.UReverse();
    if (flags & 2)
        s.VReverse();
    if (flags & 4)
        s.ExchangeUV();
}
std::string Quote(std::string s)
{
    std::string r = "\"";
    for (char c : s)
    {
        if (c == '\"')
            r += '\"';
        r += (c == '\n' || c == '\r') ? ' ' : c;
    }
    return r + '\"';
}
} // namespace
std::vector<ExtremeSpec> ExtremeSpecs()
{
    std::vector<ExtremeSpec> out;
    auto add = [&](ExtremeSpec s, const std::string &group, const std::string &label) {
        s.pair = static_cast<unsigned>(out.size());
        s.group = group;
        s.label = label;
        out.push_back(s);
    };
    for (unsigned t = 1; t < 5; ++t)
        for (double rho : {1e2, 1e4, 1e6})
            for (int level = 0; level < 4; ++level)
            {
                double chi = level == 1 ? .1 : (level == 2 ? 1 : 10);
                double theta =
                    level == 0 ? .0051 / rho : 4 * std::asin(std::sqrt(chi * 1e-5 / (2 * rho)));
                add(Base(t, rho, theta), "radius_arc",
                    level == 0 ? "fixed_arc_0.0051" : "sagitta_ratio_" + std::to_string(chi));
            }
    for (unsigned t = 1; t < 5; ++t)
        for (double factor : {.8, .99, 1.01, 1.2})
            add(Base(t, 1e4, factor * 1e-3), "angle_gate",
                "sweep_" + std::to_string(factor) + "e-3");
    for (unsigned t : {3u, 4u})
        for (double chi : {.1, 1., 10.})
        {
            auto s = Base(t, 1e4, 4 * std::asin(std::sqrt(chi * 1e-5 / (2 * 1e4))));
            double dv = 4 * std::asin(std::sqrt(chi * 1e-5 / (2 * s.truth.radius)));
            double mid = t == 3 ? Pi / 2 : 0;
            s.domain[2] = mid - dv / 2;
            s.domain[3] = mid + dv / 2;
            add(s, "weak_both", "two_weak_curvatures_" + std::to_string(chi));
        }
    for (double beta : {1e-2, 1e-4, 1e-6})
    {
        auto s = Base(2, 100, .01);
        s.truth.halfAngle = beta;
        s.domain[2] = 100 / std::tan(beta);
        s.domain[3] = s.domain[2] + 100;
        add(s, "shallow_cone", "half_angle_" + std::to_string(beta));
    }
    for (double ratio : {1e6, 1.001})
    {
        auto s = Base(4, 100, .01);
        s.truth.radius = 1;
        s.truth.majorRadius = ratio;
        if (ratio > 2)
        {
            s.domain = {-.005, .005, -.1, .1};
            s.rho = ratio + std::cos(.1);
            s.theta = .01;
        }
        else
        {
            s.domain = {-.5, .5, Pi - .2, Pi + .2};
            s.rho = ratio + std::cos(Pi - .2);
            s.theta = 1;
        }
        add(s, "torus_ratio", "major_minor_" + std::to_string(ratio));
    }
    for (unsigned t : {2u, 3u})
        for (double radius : {.006, .06})
        {
            auto s = Base(t, 1, 1);
            if (t == 2)
            {
                s.domain[2] = radius / std::tan(s.truth.halfAngle);
                s.domain[3] = s.domain[2] + .1;
            }
            else
            {
                s.truth.radius = 100;
                s.domain[2] = std::asin(radius / 100);
                s.domain[3] = s.domain[2] + .001;
            }
            s.rho = radius;
            add(s, "near_singular", "small_boundary_radius_" + std::to_string(radius));
        }
    for (unsigned t : {1u, 3u, 4u})
        for (double offset : {0., 1e6, 1e9})
        {
            auto s = Base(t, 100, .01);
            s.truth.origin = Point({offset, offset, offset});
            add(s, "translation", "offset_" + std::to_string(offset));
        }
    for (unsigned t : {0u, 1u})
        for (double ratio : {1e6, 1e8})
        {
            auto s = Base(t, 100, .01);
            if (t == 0)
                s.domain = {0, 1, 0, 1};
            s.knotRatio = ratio;
            add(s, "knots", "span_ratio_" + std::to_string(ratio));
        }
    for (unsigned t : {1u, 3u, 4u})
    {
        auto s = Base(t, 100, .01);
        s.orientation = 7;
        add(s, "orientation", "reverse_both_exchange_uv");
    }
    for (unsigned t : {1u, 3u, 4u})
    {
        auto s = Base(t, 1e4, .0012);
        s.domain[0] += 2 * Pi;
        s.domain[1] += 2 * Pi;
        if (t == 4)
        {
            s.domain[2] += 2 * Pi;
            s.domain[3] += 2 * Pi;
        }
        add(s, "seam", "cross_periodic_seam");
    }
    for (unsigned t : {0u, 1u})
    {
        auto s = Base(t, 1e6, .0051 / 1e6);
        s.domain[2] = 0;
        s.domain[3] = 49000;
        if (t == 0)
            s.domain[0] = 0, s.domain[1] = .0051;
        add(s, "aspect", "edge_0.0051_by_49000");
    }
    if (out.size() != 100)
        throw std::runtime_error("Extreme matrix must contain 100 pairs");
    return out;
}
std::string ExtremeFilename(unsigned id)
{
    std::ostringstream s;
    s << "extreme-" << std::setfill('0') << std::setw(3) << id << ".b2a";
    return s.str();
}
double ExtremeMotherDistance(const AnalyticParameters &p, V world)
{
    // 轴线与原点运算显式实现，不调用生产投影代码。
    V z = Unit(Get(p.axis)), d = world - Get(p.origin);
    double height = Dot(d, z), rho = Norm(d - z * height);
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        return std::abs(height);
    case sggk::SurfaceType::Cylinder:
        return std::abs(rho - p.radius);
    case sggk::SurfaceType::Sphere:
        return std::abs(Norm(d) - p.radius);
    case sggk::SurfaceType::Cone:
        return rho * std::sin(p.halfAngle) + height * std::cos(p.halfAngle) < 0
                   ? Norm(d)
                   : std::abs(rho * std::cos(p.halfAngle) - height * std::sin(p.halfAngle));
    case sggk::SurfaceType::Torus:
        return std::abs(std::hypot(rho - p.majorRadius, height) - p.radius);
    default:
        return INFINITY;
    }
}
V ExtremeMotherNormal(const AnalyticParameters &p, V world)
{
    V z = Unit(Get(p.axis)), d = world - Get(p.origin);
    double height = Dot(d, z);
    V radial = d - z * height;
    switch (p.type)
    {
    case sggk::SurfaceType::Plane:
        return z;
    case sggk::SurfaceType::Cylinder:
        return Unit(radial);
    case sggk::SurfaceType::Sphere:
        return Unit(d);
    case sggk::SurfaceType::Cone:
        return Unit(radial) * std::cos(p.halfAngle) - z * std::sin(p.halfAngle);
    case sggk::SurfaceType::Torus:
        return Unit(radial - Unit(radial) * p.majorRadius + z * height);
    default:
        throw std::runtime_error("Unknown mother type");
    }
}
ExtremeQualification QualifyExtreme(const DatasetCase &c, bool business)
{
    ExtremeQualification q;
    const auto &s = *c.input;
    const auto &cp = s.ControlPoints();
    q.controls = static_cast<unsigned>(cp.size() * cp[0].size());
    q.pair = BoundPair(*c.source, s, 1e-6);
    q.lengths = EdgeLengths(s, 32);
    auto coarse = EdgeLengths(s, 16);
    for (unsigned e = 0; e < 4; ++e)
        q.lengthDisagreement = std::max(q.lengthDisagreement, std::abs(coarse[e] - q.lengths[e]));
    if (!s.Weights().empty())
    {
        q.minWeight = INFINITY;
        q.maxWeight = 0;
        for (const auto &row : s.Weights())
            for (double w : row)
            {
                q.minWeight = std::min(q.minWeight, w);
                q.maxWeight = std::max(q.maxWeight, w);
            }
    }
    for (const auto &row : cp)
        for (const auto &p : row)
            for (double x : {p.X(), p.Y(), p.Z()})
                q.ulp = std::max(q.ulp, std::nextafter(std::abs(x), INFINITY) - std::abs(x));
    bool lengthOK = true;
    for (unsigned e = 0; e < 4; ++e)
    {
        double length = q.lengths[e];
        if (business && (c.collapsedMask & (1u << e)))
        {
            // 要求控制点精确重合，不能把数值很短的边当作点退化边。
            auto first = cp[e == 1 ? cp.size()-1 : 0][e == 3 ? cp[0].size()-1 : 0];
            for (size_t j = 0; j < (e < 2 ? cp[0].size() : cp.size()); ++j)
            {
                auto p = cp[e < 2 ? (e == 0 ? 0 : cp.size()-1) : j]
                           [e < 2 ? j : (e == 2 ? 0 : cp[0].size()-1)];
                lengthOK &= p.X() == first.X() && p.Y() == first.Y() && p.Z() == first.Z();
            }
        }
        else lengthOK &= std::isfinite(length) &&
            length >= (business ? .01*(1-1e-10) : .005) &&
            length <= (business ? 1000*(1+1e-10) : 5e4);
    }
    q.dimensionsPassed = !business || BusinessDimensions(c);
    bool anglesOK = true;
    if (business && !c.collapsedMask)
    {
        q.sourceAngles = BusinessBoundaryAngles(*c.source,s);
        anglesOK = q.sourceAngles[0]<.05 && q.sourceAngles[1]<.05;
    }
    else if (c.collapsedMask) q.sourceAngles.fill(std::numeric_limits<double>::quiet_NaN());
    q.passed = q.controls >= 50 && q.controls <= 200 && q.pair.passed && lengthOK && q.dimensionsPassed && anglesOK &&
               (!c.refitted || !s.IsRational());
    q.message = q.passed ? "qualified_pair_bound_and_measured_lengths"
                         : (!anglesOK ? "source_angle_sample_exceeded" : !q.dimensionsPassed ? "business_dimensions_out_of_range" : !q.pair.passed ? "pair_certificate_unresolved"
                                           : (!lengthOK ? "measured_edge_length_out_of_range"
                                                        : "representation_or_control_budget"));
    return q;
}
void GenerateExtremeDataset(const std::filesystem::path &dir, bool business)
{
    if (std::filesystem::exists(dir) && !std::filesystem::is_empty(dir))
        throw std::runtime_error(
            "Generate into an empty directory; existing fixtures are never overwritten");
    std::filesystem::create_directories(dir);
    if (business) { std::ofstream profile(dir / "profile.txt"); profile << "business-extreme200-v1\n"; }
    unsigned failed = 0;
    std::ofstream csv(dir / "manifest.csv");
    if (!csv)
        throw std::runtime_error("Cannot write manifest");
    csv << std::setprecision(17)
        << "id,pair,group,label,type,source_kind,rho,theta,sagitta_ratio,controls,rational,min_"
           "weight,max_weight,pair_bound,pair_pieces,edge0,edge1,edge2,edge3,length_disagreement,"
           "coordinate_ulp,qualified,message,profile,dimensions_passed,collapsed_mask,control_bucket,source_tangent_sampled_deg,source_normal_sampled_deg\n";
    for (const auto &spec : (business ? BusinessExtremeSpecs() : ExtremeSpecs()))
        for (unsigned mode = 0; mode < 2; ++mode)
        {
            unsigned id = 2 * spec.pair + mode;
            DatasetCase c;
            c.id = id;
            c.pair = spec.pair;
            c.type = spec.type;
            c.variant = spec.orientation;
            c.bucket = spec.pair % 5;
            c.truth = spec.truth;
            c.analyticDomain = spec.domain;
            c.refitted = mode != 0;
            c.seed = 0xe2002026u + id;
            ExtremeQualification q;
            try
            {
                AnalyticPatch patch;
                patch.parameters = c.truth;
                patch.domain = c.analyticDomain;
                MakeNative(patch);
                c.source = NormalizeDatasetSurface(*patch.surface->ToBSpline(*patch.sdkDomain));
                auto snap = [&](sggk::BSplineSurface &surface) {
                    if (!spec.collapsed) return;
                    auto point = c.source->ControlPoints()[0][0];
                    for (unsigned i = 0; i < surface.ControlPoints().size(); ++i)
                        surface.SetControlPoint(point,i,0);
                };
                snap(*c.source);
                if (mode)
                {
                    for (unsigned degree : {3u, 5u, 7u, 9u, 11u, 13u})
                    {
                        unsigned qdegree = c.type >= 3 ? degree : 3;
                        c.input = FitDatasetPolynomial(*c.source, degree, qdegree);
                        snap(*c.input);
                        // 优先选择能通过连续上界认证的最低次数。
                        if (BoundPair(*c.source, *c.input, 7e-7).passed)
                        {
                            auto angles = business && !spec.collapsed ?
                                BusinessBoundaryAngles(*c.source,*c.input) : std::array<double,2>{};
                            if (angles[0]<.025 && angles[1]<.025) break;
                        }
                    }
                }
                else
                    c.input = std::make_shared<sggk::BSplineSurface>(*c.source);
                if (spec.knotRatio > 1)
                {
                    c.input->InsertUKnots(1 / spec.knotRatio, 1);
                    c.input->InsertVKnots(1 / spec.knotRatio, 1);
                }
                if (business) PadBusinessExtreme(*c.input,c.bucket);
                else {
                unsigned nu = std::max(8u, static_cast<unsigned>(c.input->ControlPoints().size()));
                unsigned nv =
                    std::max(8u, static_cast<unsigned>(c.input->ControlPoints()[0].size()));
                if (nu * nv > 200)
                    throw std::runtime_error("Control count exceeds 200 before padding");
                PadDatasetSurface(*c.input, nu, nv, 1);
                }
                snap(*c.input);
                OrientSurface(*c.source, spec.orientation);
                OrientSurface(*c.input, spec.orientation);
                c.source = NormalizeDatasetSurface(*c.source);
                c.input = NormalizeDatasetSurface(*c.input);
                c.collapsedMask = spec.collapsed ? CollapsedMask(*c.input) : 0;
                if (spec.collapsed && !c.collapsedMask) throw std::runtime_error("Collapsed edge lost");
                q = QualifyExtreme(c,business);
                c.storedBound = std::isfinite(q.pair.upper) ? q.pair.upper : 0;
                SaveDatasetCase(c, dir / ExtremeFilename(id));
            }
            catch (const std::exception &e)
            {
                q.passed = false;
                q.message = e.what();
            }
            csv << id << ',' << c.pair << ',' << spec.group << ',' << spec.label << ','
                << TypeName(c.truth.type) << ',' << (mode ? "polynomial" : "analytic_nurbs") << ','
                << spec.rho << ',' << spec.theta << ','
                << 2 * spec.rho * std::pow(std::sin(spec.theta / 4), 2) / 1e-5 << ',' << q.controls
                << ',' << (c.input && c.input->IsRational()) << ',' << q.minWeight << ','
                << q.maxWeight << ',' << q.pair.upper << ',' << q.pair.pieces;
            for (double x : q.lengths)
                csv << ',' << x;
            csv << ',' << q.lengthDisagreement << ',' << q.ulp << ',' << q.passed << ','
                << Quote(q.message) << ',' << (business ? "business-extreme200-v1" : "legacy-extreme200")
                << ',' << q.dimensionsPassed << ',' << c.collapsedMask << ',' << c.bucket
                << ',' << q.sourceAngles[0] << ',' << q.sourceAngles[1] << '\n';
            failed += !q.passed;
            csv.flush();
            std::cout << "generated " << id << " " << q.message << std::endl;
        }
    if (!csv)
        throw std::runtime_error("Manifest write failed");
    if (failed) throw std::runtime_error(std::to_string(failed)+" fixture qualifications failed; records preserved");
}
} // namespace bspline2analytic::tests
