#include "ExtremeDataset.h"
#include "ApproximateInternal.h"
#include <fstream>
namespace bspline2analytic::tests
{
using namespace detail;
namespace
{
ExtremeSpec Base(unsigned type, double radius, double sweep)
{
    ExtremeSpec s;
    s.type = type; s.truth.type = DatasetType(type);
    s.truth.radius = radius; s.truth.halfAngle = Pi / 4;
    s.rho = radius; s.theta = sweep; s.domain = {-sweep/2, sweep/2, 0, 1};
    if (type == 2)
    {
        s.rho = std::min(radius, 999.);
        s.domain[2] = s.rho; s.domain[3] = s.rho + 1;
    }
    if (type == 3)
    {
        s.domain[2] = Pi/2 - .01; s.domain[3] = Pi/2 + .01;
        s.rho = radius * std::cos(.01);
    }
    if (type == 4)
    {
        s.truth.radius = .1 * radius;
        s.truth.majorRadius = .9 * radius;
        s.domain[2] = -.01; s.domain[3] = .01;
        s.rho = s.truth.majorRadius + s.truth.radius * std::cos(.01);
    }
    return s;
}
void Sweep(ExtremeSpec &s, double theta)
{ s.theta = theta; s.domain[0] = -theta/2; s.domain[1] = theta/2; }
bool InRange(double x)
{
    // Only rounding allowance for generated length measurements, not algorithm tolerance.
    return std::isfinite(x) && x >= .01 * (1-1e-10) && x <= 1000 * (1+1e-10);
}
}
std::vector<ExtremeSpec> BusinessExtremeSpecs()
{
    std::vector<ExtremeSpec> out;
    auto add = [&](ExtremeSpec s, const std::string &group, const std::string &label) {
        s.pair = static_cast<unsigned>(out.size()); s.group = group; s.label = label;
        out.push_back(s);
    };
    // 30 pairs: six length scales, five families. Keep all noncollapsed edges in range.
    for (unsigned type = 0; type < 5; ++type)
        for (double size : {.01, .1, 1., 10., 100., 1000.})
        {
            auto s = Base(type, size, size == .01 ? 1.1 : .6);
            double length = size == .01 ? .011 : .6 * size;
            if (type == 0) s.domain = {0, length, 0, length};
            if (type == 1) s.domain[3] = length;
            if (type == 2)
            {
                s.rho = size == .01 ? size : .6 * size;
                s.domain[2] = s.rho; s.domain[3] = s.rho + (size == .01 ? .011 : .3*size);
            }
            if (type == 3)
            {
                double dv = size == .01 ? 1.1 : .6;
                Sweep(s, size == .01 ? 1.3 : .6);
                s.domain[2] = Pi/2-dv/2; s.domain[3] = Pi/2+dv/2;
                s.rho = size * std::cos(dv/2);
            }
            if (type == 4)
            {
                s.truth.radius = size == 1000 ? 200 : size;
                s.truth.majorRadius = size == 1000 ? 1000 : 3*size;
                double dv = size == .01 ? 1.2 : .6;
                Sweep(s, .6); s.domain[2] = -dv/2; s.domain[3] = dv/2;
                s.rho = s.truth.majorRadius+s.truth.radius*std::cos(dv/2);
            }
            add(s, "scale_regular", "scale_" + std::to_string(size));
        }
    // 24 pairs: difficulty determined by sagitta / final distance tolerance.
    for (unsigned type = 1; type < 5; ++type)
        for (double radius : {100., 1000.})
            for (double chi : {.1, 1., 10.})
            {
                auto s = Base(type, radius, .01);
                Sweep(s, 4*std::asin(std::sqrt(chi*1e-5/(2*s.rho))));
                add(s, "short_arc", "radius_"+std::to_string(radius)+"_chi_"+std::to_string(chi));
            }
    // Eight pairs bracket the former fixed angular gate at a permitted radius.
    for (unsigned type = 1; type < 5; ++type)
        for (double factor : {.99, 1.01})
            add(Base(type, 1000, factor*1e-3), "angle_gate", "sweep_"+std::to_string(factor)+"e-3");
    for (unsigned type : {3u, 4u})
        for (double chi : {.1, 1., 10.})
        {
            auto s = Base(type, 1000, .01);
            double dv = 4*std::asin(std::sqrt(chi*1e-5/(2*s.truth.radius)));
            double mid = type == 3 ? Pi/2 : 0;
            s.domain[2] = mid-dv/2; s.domain[3] = mid+dv/2;
            s.rho = type == 3 ? s.truth.radius*std::cos(dv/2) :
                s.truth.majorRadius+s.truth.radius*std::cos(dv/2);
            Sweep(s, 4*std::asin(std::sqrt(chi*1e-5/(2*s.rho))));
            add(s, "weak_both", "chi_"+std::to_string(chi));
        }
    for (double beta : {.01, .001})
    {
        auto s = Base(2, 1, .2); s.truth.halfAngle = beta;
        s.domain[2] = 100; s.domain[3] = 200; s.rho = 100*std::tan(beta);
        add(s, "shallow_cone", "half_angle_"+std::to_string(beta));
    }
    for (unsigned mode = 0; mode < 2; ++mode)
    {
        auto s = Base(4, 100, .01);
        s.truth.radius = mode ? 1 : .01; s.truth.majorRadius = mode ? 1.01 : 1000;
        s.domain = mode ? std::array<double,4>{-.5,.5,Pi-.1,Pi+.1} :
                          std::array<double,4>{-.0005,.0005,-.6,.6};
        s.theta = s.domain[1]-s.domain[0];
        s.rho = s.truth.majorRadius+s.truth.radius*std::cos(s.domain[2]);
        add(s, "torus_ratio", mode ? "inner_gap_0.01" : "major_minor_100000");
    }
    for (unsigned type : {2u, 3u})
        for (double rho : {.011, .06})
        {
            auto s = Base(type, 100, 1);
            s.domain[2] = type == 2 ? rho : std::asin(rho/100);
            s.domain[3] = s.domain[2]+(type == 2 ? .1 : .001); s.rho = rho;
            add(s, "near_singular", "boundary_radius_"+std::to_string(rho));
        }
    for (unsigned type : {2u, 3u})
        for (double size : {.1, 100.})
        {
            auto s = Base(type, size, 1); s.domain[2] = 0;
            s.domain[3] = type == 2 ? size : .2; s.collapsed = true; s.rho = 0;
            add(s, "collapsed", "scale_"+std::to_string(size));
        }
    for (unsigned type : {1u, 3u, 4u})
        for (double offset : {1e6, 1e9})
        {
            auto s = Base(type, 100, .01);
            s.truth.origin = Point({offset,offset,offset});
            add(s, "translation", "offset_"+std::to_string(offset));
        }
    for (unsigned type : {0u, 1u})
        for (double ratio : {1e6, 1e8})
        {
            auto s = Base(type, 100, .01); if (type == 0) s.domain = {0,1,0,1};
            s.knotRatio = ratio; add(s,"knots","span_ratio_"+std::to_string(ratio));
        }
    for (unsigned type : {1u, 3u, 4u})
    {
        auto s = Base(type,100,.01); s.orientation = 7;
        // Exercise spatial rotation as well as UV reversal/exchange.
        s.truth.axis = Direction(Unit(V{1,2,3}));
        s.truth.xDirection = Direction(Perpendicular(Unit(V{1,2,3})));
        add(s,"orientation","rotated_reverse_exchange");
    }
    for (unsigned type : {1u, 3u, 4u})
    {
        auto s = Base(type,100,.02); s.domain[0] += 2*Pi; s.domain[1] += 2*Pi;
        if (type == 4) { s.domain[2] += 2*Pi; s.domain[3] += 2*Pi; }
        add(s,"seam","cross_periodic_seam");
    }
    for (unsigned type : {0u, 1u})
        for (double length : {100., 1000.})
        {
            auto s = Base(type,1000,.011/1000); s.domain[2] = 0; s.domain[3] = length;
            if (type == 0) s.domain[0] = 0, s.domain[1] = .011;
            add(s,"aspect","short_0.011_long_"+std::to_string(length));
        }
    if (out.size() != 100) throw std::runtime_error("Business matrix must contain 100 pairs");
    return out;
}
bool IsBusinessExtreme(const std::filesystem::path &dir)
{
    if (!std::filesystem::exists(dir/"profile.txt")) return false;
    std::ifstream file(dir/"profile.txt"); std::string tag; std::getline(file,tag);
    if (tag != "business-extreme200-v1") throw std::runtime_error("Unknown extreme matrix profile");
    return true;
}
bool BusinessDimensions(const DatasetCase &c)
{
    const auto &a = c.truth; const auto &d = c.analyticDomain;
    if (c.type == 0) return InRange(d[1]-d[0]) && InRange(d[3]-d[2]);
    if (c.type == 1) return InRange(a.radius) && InRange(d[3]-d[2]);
    if (c.type == 2)
    {
        double lo = d[2]*std::tan(a.halfAngle), hi = d[3]*std::tan(a.halfAngle);
        return (c.collapsedMask ? lo == 0 : InRange(lo)) && InRange(hi) &&
            InRange(d[3]-d[2]) && d[3] <= 1000*(1+1e-10);
    }
    if (c.type == 3) return InRange(a.radius);
    return InRange(a.radius) && InRange(a.majorRadius) && a.majorRadius > a.radius;
}
std::array<double,2> BusinessBoundaryAngles(const sggk::BSplineSurface &source,
                                                 const sggk::BSplineSurface &input)
{
    std::array<double,2> worst{};
    auto angle = [](V a, V b) {
        if (!(Norm(a)>0 && Norm(b)>0)) return std::numeric_limits<double>::infinity();
        a = Unit(a); b = Unit(b);
        return std::atan2(Norm(Cross(a,b)),std::abs(Dot(a,b)))*180/Pi;
    };
    for (unsigned e=0;e<4;++e)
    {
        const auto &knots=e<2?input.KnotsV():input.KnotsU();
        for (size_t k=1;k<knots.size();++k)
            for (unsigned j=0;j<=32;++j)
            {
                double t=knots[k-1]+(knots[k]-knots[k-1])*j/32;
                double u=e<2?(e==0?0.:1.):t, v=e<2?t:(e==2?0.:1.);
                auto su=u==1?sggk::DerivType::LEFT:sggk::DerivType::RIGHT;
                auto sv=v==1?sggk::DerivType::LEFT:sggk::DerivType::RIGHT;
                sggk::Point3D point;
                auto a=source.CalcPartialD1(u,v,su,sv,point);
                auto b=input.CalcPartialD1(u,v,su,sv,point);
                worst[0]=std::max(worst[0],angle(Get(a[e<2?1:0]),Get(b[e<2?1:0])));
                worst[1]=std::max(worst[1],angle(Cross(Get(a[0]),Get(a[1])),Cross(Get(b[0]),Get(b[1]))));
            }
    }
    return worst;
}
void PadBusinessExtreme(sggk::BSplineSurface &s, unsigned bucket)
{
    const unsigned lows[]{50,80,110,140,170}, highs[]{79,109,139,169,200};
    unsigned u0 = static_cast<unsigned>(s.ControlPoints().size());
    unsigned v0 = static_cast<unsigned>(s.ControlPoints()[0].size());
    unsigned bestU = 0, bestV = 0, best = bucket == 4 ? 0 : 201;
    // Prefer the assigned band; if degree requires more controls, record a larger band.
    for (unsigned b = bucket; b < 5 && !bestU; ++b)
        for (unsigned u = u0; u <= 25; ++u)
            for (unsigned v = v0; v <= 25; ++v)
                if (u*v >= lows[b] && u*v <= highs[b] && (bucket == 4 ? u*v > best : u*v < best))
                { bestU = u; bestV = v; best = u*v; }
    if (!bestU) throw std::runtime_error("No control grid within business budget");
    PadDatasetSurface(s,bestU,bestV,1);
}
}
