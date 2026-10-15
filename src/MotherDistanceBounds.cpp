#include "MotherDistanceInternal.h"
#include <functional>

namespace bspline2analytic::detail::mother
{
namespace
{
double Guard(double x) { return 512 * std::numeric_limits<double>::epsilon() * x; }
double Choose(unsigned n, unsigned k)
{
    double x = 1;
    for (unsigned i = 1; i <= k; ++i) x *= double(n + 1 - i) / i;
    return x;
}
// Candidate-local cached convolution factors; constant weights never become
// dense polynomials. Subdivision operates on residuals already constructed.
class Products
{
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> cache;
    const std::vector<double> &Weights(unsigned p, unsigned q)
    {
        auto key = std::make_pair(p, q);
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
        std::vector<double> w((p + 1) * (q + 1));
        for (unsigned i = 0; i <= p; ++i)
            for (unsigned j = 0; j <= q; ++j)
                w[i * (q + 1) + j] = Choose(p, i) * Choose(q, j) / Choose(p + q, i + j);
        return cache.emplace(key, std::move(w)).first->second;
    }
  public:
    Poly Mul(const Poly &a, const Poly &b)
    {
        Poly c(a.u + b.u, a.v + b.v);
        const auto &u = Weights(a.u, b.u), &v = Weights(a.v, b.v);
        for (unsigned i = 0; i <= a.u; ++i)
            for (unsigned j = 0; j <= a.v; ++j)
                for (unsigned k = 0; k <= b.u; ++k)
                    for (unsigned l = 0; l <= b.v; ++l)
                        c.at(i + k, j + l) += a.at(i, j) * b.at(k, l) *
                            u[i * (b.u + 1) + k] * v[j * (b.v + 1) + l];
        return c;
    }
    Poly Elevate(Poly a, unsigned u, unsigned v)
    {
        // Repeated one-degree elevation is linear in coefficient count and
        // avoids multiplication by a dense polynomial representing one.
        while (a.u < u || a.v < v)
        {
            bool along = a.u < u;
            Poly b(a.u + unsigned(along), a.v + unsigned(!along));
            for (unsigned i = 0; i <= b.u; ++i)
                for (unsigned j = 0; j <= b.v; ++j)
                {
                    unsigned k = along ? i : j, n = along ? b.u : b.v;
                    double t = double(k) / n;
                    if (k) b.at(i, j) += t * a.at(i - unsigned(along), j - unsigned(!along));
                    if (k < n) b.at(i, j) += (1 - t) * a.at(i, j);
                }
            a = std::move(b);
        }
        return a;
    }
};
struct Residual
{
    Poly f, radial, z;
    double error = 0, radialError = 0, zError = 0, enclosure = 0;
};
Residual Build(const BezierCell &cell, V sourceOrigin, const paired::LocalFrame &frame,
               const AnalyticParameters &p, double scale, Products &products)
{
    Poly x(cell.p.x.u, cell.p.x.v), y = x, z = x;
    double enclosure = 0;
    for (size_t k = 0; k < x.c.size(); ++k)
    {
        IV v = IV{cell.p.x.c[k], cell.p.y.c[k], cell.p.z.c[k]} + Exact(sourceOrigin) - frame.origin;
        std::array<I, 3> c{Dot(v, frame.x) / I(scale), Dot(v, frame.y) / I(scale),
                           Dot(v, frame.z) / I(scale)};
        std::array<double, 3> mid{}, radius{};
        for (unsigned d = 0; d < 3; ++d)
        {
            if (!c[d].finite()) throw Failure(Status::VerificationFailed, "Nonfinite local enclosure");
            mid[d] = c[d].lo * .5 + c[d].hi * .5;
            radius[d] = (c[d] - I(mid[d])).absMax();
        }
        x.c[k] = mid[0]; y.c[k] = mid[1]; z.c[k] = mid[2];
        enclosure = std::max(enclosure,
            Sqrt(Sqr(I(radius[0])) + Sqr(I(radius[1])) + Sqr(I(radius[2]))).hi);
    }
    Residual out;
    out.z = z; out.enclosure = enclosure;
    out.zError = Guard(z.maxAbs());
    if (p.type == sggk::SurfaceType::Plane)
    {
        out.f = z;
        out.error = Guard(x.maxAbs() + y.maxAbs() + z.maxAbs());
        return out;
    }
    out.radial = products.Mul(x, x) + products.Mul(y, y);
    out.radialError = Guard(x.maxAbs() * x.maxAbs() + y.maxAbs() * y.maxAbs());
    Poly zz = products.Mul(z, z), norm = out.radial + zz;
    double r = p.radius / scale, R = p.majorRadius / scale;
    if (p.type == sggk::SurfaceType::Sphere || p.type == sggk::SurfaceType::Cylinder)
    {
        out.f = p.type == sggk::SurfaceType::Sphere ? norm : out.radial;
        out.error = Guard(out.f.maxAbs() + r * r);
        for (auto &c : out.f.c) c -= r * r;
    }
    else if (p.type == sggk::SurfaceType::Cone)
    {
        double ca = std::cos(p.halfAngle), sa = std::sin(p.halfAngle);
        Poly a = out.radial * (ca * ca), b = zz * (sa * sa);
        out.f = a - b;
        out.error = Guard(a.maxAbs() + b.maxAbs());
    }
    else
    {
        Poly sum = norm;
        for (auto &c : sum.c) c += R * R - r * r;
        Poly a = products.Mul(sum, sum);
        Poly b = products.Elevate(out.radial, a.u, a.v) * (4 * R * R);
        out.f = a - b;
        out.error = Guard(a.maxAbs() + b.maxAbs() +
            std::pow(norm.maxAbs() + R * R + r * r, 2));
    }
    if (!std::isfinite(out.f.maxAbs() + out.error))
        throw Failure(Status::VerificationFailed, "Nonfinite mother residual");
    return out;
}
double Bound(const Residual &f, const AnalyticParameters &p, double scale)
{
    double denominator = 1, r = p.radius / scale;
    if (p.type == sggk::SurfaceType::Sphere || p.type == sggk::SurfaceType::Cylinder)
        denominator = r;
    else if (p.type == sggk::SurfaceType::Cone || p.type == sggk::SurfaceType::Torus)
    {
        double rho = std::sqrt(std::max(0., f.radial.min() - f.radialError));
        double zlo = f.z.min() - f.zError, zhi = f.z.max() + f.zError;
        if (p.type == sggk::SurfaceType::Cone)
        {
            if (zlo < 0) return INFINITY;
            denominator = rho * std::cos(p.halfAngle) + zlo * std::sin(p.halfAngle);
        }
        else
        {
            double R = p.majorRadius / scale;
            double zmin = zlo <= 0 && zhi >= 0 ? 0 : std::min(zlo * zlo, zhi * zhi);
            // For supported ring tori R>r the other factor is positive.
            double positive = (rho + R) * (rho + R) + zmin;
            denominator = (positive - r * r - Guard(positive + r * r)) * r;
        }
    }
    denominator -= Guard(std::abs(denominator));
    if (!(denominator > 0)) return INFINITY;
    // Distance to a closed set is 1-Lipschitz. Keep the source coefficient
    // enclosure lost by midpoint conversion; convex hull bounds it everywhere.
    double value = scale * ((f.f.maxAbs() + f.error) / denominator + f.enclosure);
    return value + Guard(std::abs(value));
}
std::pair<Residual, Residual> Subdivide(const Residual &f, bool along)
{
    auto r = Split(f.f, along), a = Split(f.radial, along), z = Split(f.z, along);
    Residual l{std::move(r.first), std::move(a.first), std::move(z.first),
        f.error + Guard(f.f.maxAbs()), f.radialError + Guard(f.radial.maxAbs()),
        f.zError + Guard(f.z.maxAbs()), f.enclosure};
    Residual h{std::move(r.second), std::move(a.second), std::move(z.second),
        l.error, l.radialError, l.zError, l.enclosure};
    return {std::move(l), std::move(h)};
}
} // namespace
bool ValidatePosition(const paired::PreparedSource &source, const AnalyticParameters &p,
                      const RecognitionOptions &options, MotherDistanceValidationReport &out,
                      ValidationTimer &timer)
{
    ValidationTimer::Scope stage(timer, out.timing.surfaceBounds);
    Products products;
    const auto frame = paired::MakeLocalFrame(p);
    double scale = std::max(p.radius, std::abs(p.majorRadius));
    for (const auto &cell : source.cells)
        for (size_t k = 0; k < cell.p.x.c.size(); ++k)
        {
            IV v = IV{cell.p.x.c[k], cell.p.y.c[k], cell.p.z.c[k]} + Exact(source.origin) - frame.origin;
            scale = std::max(scale, Sqrt(Sqr(v.x) + Sqr(v.y) + Sqr(v.z)).hi);
        }
    if (!(scale > 0 && std::isfinite(scale)))
        throw Failure(Status::VerificationFailed, "Invalid mother normalization scale");
    double maximum = 0;
    std::function<bool(const Residual &, unsigned)> visit;
    visit = [&](const Residual &f, unsigned depth) {
        if (++out.surfaceNodes > options.maxValidationPieces)
            throw paired::Budget("Mother distance node budget exhausted");
        out.surfaceDepth = std::max(out.surfaceDepth, depth);
        double upper;
        {
            Measure time(out.residualBoundMs);
            upper = Bound(f, p, scale);
        }
        if (upper < options.distanceTolerance)
        {
            maximum = std::max(maximum, upper);
            return true;
        }
        if (depth >= options.maxSubdivisionDepth) return false;
        std::pair<Residual, Residual> children;
        {
            ValidationTimer::Scope split(timer, out.timing.subdivision);
            children = Subdivide(f, depth % 2 == 0);
        }
        return visit(children.first, depth + 1) && visit(children.second, depth + 1);
    };
    for (const auto &cell : source.cells)
    {
        Residual f;
        {
            Measure time(out.residualBuildMs);
            f = Build(cell, source.origin, frame, p, scale, products);
        }
        out.sourceEnclosureRadius = std::max(out.sourceEnclosureRadius, f.enclosure * scale);
        if (!visit(f, 0))
        {
            out.motherDistance.state = CheckState::Unresolved;
            return false;
        }
    }
    out.motherDistance.state = CheckState::Passed;
    out.motherDistance.upperBound = maximum;
    return true;
}
} // namespace bspline2analytic::detail::mother
