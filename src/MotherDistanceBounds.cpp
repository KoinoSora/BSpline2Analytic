#include "MotherDistanceInternal.h"
#include <functional>

namespace bspline2analytic::detail::mother
{
namespace
{
double Guard(double x)
{
    return 512 * std::numeric_limits<double>::epsilon() * x;
}
double Choose(unsigned n, unsigned k)
{
    double x = 1;
    for (unsigned i = 1; i <= k; ++i)
        x *= double(n + 1 - i) / i;
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
        if (it != cache.end())
            return it->second;
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
                        c.at(i + k, j + l) +=
                            a.at(i, j) * b.at(k, l) * u[i * (b.u + 1) + k] * v[j * (b.v + 1) + l];
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
                    if (k)
                        b.at(i, j) += t * a.at(i - unsigned(along), j - unsigned(!along));
                    if (k < n)
                        b.at(i, j) += (1 - t) * a.at(i, j);
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
struct Projection
{
    V offset, offsetError;
    std::array<std::array<double, 3>, 3> axes{};
    std::array<double, 3> rowNorm{}, rowError{};
    Projection(V sourceOrigin, const AnalyticParameters &p)
    {
        offset = sourceOrigin - Get(p.origin);
        offsetError = {RoundedError(offset.x), RoundedError(offset.y), RoundedError(offset.z)};
        if (p.type == sggk::SurfaceType::Sphere)
        {
            for (unsigned d = 0; d < 3; ++d)
            {
                axes[d][d] = 1;
                rowNorm[d] = 1;
            }
            return;
        }
        // Parameter-only interval work, once per candidate. Never rotate a
        // source coordinate interval. Two transverse projections give rho^2
        // with two squares, cheaper than three Cartesian perpendicular squares.
        auto f = paired::MakeLocalFrame(p);
        std::array<IV, 3> basis{f.x, f.y, f.z};
        for (unsigned d = 0; d < 3; ++d)
        {
            std::array<I, 3> row{basis[d].x, basis[d].y, basis[d].z};
            for (unsigned k = 0; k < 3; ++k)
            {
                axes[d][k] = row[k].lo * .5 + row[k].hi * .5;
                rowNorm[d] = ErrorAdd(rowNorm[d], std::abs(axes[d][k]));
                rowError[d] = ErrorAdd(rowError[d], (row[k] - I(axes[d][k])).absMax());
            }
        }
    }
};
Residual Build(const SourceCell &cell, const Projection &frame, const AnalyticParameters &p,
               double scale, Products &products)
{
    Poly x(cell.x.u, cell.x.v), y = x, z = x;
    double enclosure = 0;
    V maxInput{};
    for (size_t k = 0; k < x.c.size(); ++k)
    {
        V v = V{cell.x.c[k], cell.y.c[k], cell.z.c[k]} + frame.offset;
        maxInput = {std::max(maxInput.x, std::abs(v.x)), std::max(maxInput.y, std::abs(v.y)),
                    std::max(maxInput.z, std::abs(v.z))};
        std::array<double, 3> a{v.x, v.y, v.z}, mid{};
        for (unsigned d = 0; d < 3; ++d)
        {
            mid[d] = p.type == sggk::SurfaceType::Sphere
                         ? a[d] / scale
                         : (a[0] * frame.axes[d][0] + a[1] * frame.axes[d][1] +
                            a[2] * frame.axes[d][2]) /
                               scale;
            if (!std::isfinite(mid[d]))
                throw Failure(Status::VerificationFailed, "Nonfinite local coefficient");
        }
        x.c[k] = mid[0];
        y.c[k] = mid[1];
        z.c[k] = mid[2];
    }
    // Uniform coefficient envelopes also enclose every point by partition of
    // unity. Scalar upward rounding is only needed here, not in the 3D loops.
    double M = std::max({maxInput.x, maxInput.y, maxInput.z});
    double r = std::max({ErrorAdd(cell.error.x, frame.offsetError.x),
                         ErrorAdd(cell.error.y, frame.offsetError.y),
                         ErrorAdd(cell.error.z, frame.offsetError.z)});
    r = ErrorAdd(r, RoundedError(M));
    std::array<double, 3> error{};
    for (unsigned d = 0; d < 3; ++d)
    {
        double A = frame.rowNorm[d], e = frame.rowError[d];
        double projectionError = ErrorAdd(ErrorAdd(ErrorMul(ErrorAdd(A, e), r), ErrorMul(e, M)),
                                          ErrorMul(ErrorMul(Gamma(6), A), M));
        projectionError = ErrorAdd(projectionError, 64 * std::numeric_limits<double>::denorm_min());
        error[d] = ErrorAdd(std::nextafter(projectionError / scale, INFINITY),
                            RoundedError(d == 0   ? x.maxAbs()
                                         : d == 1 ? y.maxAbs()
                                                  : z.maxAbs()));
    }
    enclosure = ErrorNorm({error[0], error[1], error[2]});
    Residual out;
    out.z = z;
    out.enclosure = enclosure;
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
    r = p.radius / scale;
    double R = p.majorRadius / scale;
    if (p.type == sggk::SurfaceType::Sphere || p.type == sggk::SurfaceType::Cylinder)
    {
        out.f = p.type == sggk::SurfaceType::Sphere ? norm : out.radial;
        out.error = Guard(out.f.maxAbs() + r * r);
        for (auto &c : out.f.c)
            c -= r * r;
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
        for (auto &c : sum.c)
            c += R * R - r * r;
        Poly a = products.Mul(sum, sum);
        Poly b = products.Elevate(out.radial, a.u, a.v) * (4 * R * R);
        out.f = a - b;
        out.error = Guard(a.maxAbs() + b.maxAbs() + std::pow(norm.maxAbs() + R * R + r * r, 2));
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
            if (zlo < 0)
                return INFINITY;
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
    if (!(denominator > 0))
        return INFINITY;
    // Distance to a closed set is 1-Lipschitz. Keep the source coefficient
    // enclosure lost by midpoint conversion; convex hull bounds it everywhere.
    double value = scale * ((f.f.maxAbs() + f.error) / denominator + f.enclosure);
    return value + Guard(std::abs(value));
}
std::pair<Residual, Residual> Subdivide(const Residual &f, bool along)
{
    auto r = Split(f.f, along), a = Split(f.radial, along), z = Split(f.z, along);
    Residual l{std::move(r.first),
               std::move(a.first),
               std::move(z.first),
               f.error + Guard(f.f.maxAbs()),
               f.radialError + Guard(f.radial.maxAbs()),
               f.zError + Guard(f.z.maxAbs()),
               f.enclosure};
    Residual h{std::move(r.second), std::move(a.second),
               std::move(z.second), l.error,
               l.radialError,       l.zError,
               l.enclosure};
    return {std::move(l), std::move(h)};
}
} // namespace
bool ValidatePosition(SourceCache &source, const AnalyticParameters &p,
                      const RecognitionOptions &options, MotherDistanceValidationReport &out,
                      ValidationTimer &timer)
{
    ValidationTimer::Scope stage(timer, out.timing.surfaceBounds);
    Products products;
    const Projection frame(source.origin, p);
    // A scale need not be a tight radius. The input control box is available
    // before extraction and avoids visiting all cell coefficients twice.
    double scale = std::max({std::abs(p.radius), std::abs(p.majorRadius),
                             source.magnitude.x + std::abs(frame.offset.x),
                             source.magnitude.y + std::abs(frame.offset.y),
                             source.magnitude.z + std::abs(frame.offset.z)});
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
        if (depth >= options.maxSubdivisionDepth)
            return false;
        std::pair<Residual, Residual> children;
        {
            ValidationTimer::Scope split(timer, out.timing.subdivision);
            children = Subdivide(f, depth % 2 == 0);
        }
        return visit(children.first, depth + 1) && visit(children.second, depth + 1);
    };
    // All active cells cover the closed source domain, including its four edges.
    // Streaming permits early failure, but Passed requires every cell to pass.
    for (size_t id = 0; id < source.cells.size(); ++id)
    {
        const SourceCell *cell;
        {
            ValidationTimer::Scope preparation(timer, out.timing.preparation);
            Measure time(out.sourceExtractionMs);
            cell = &source.Get(id, out);
            if (ErrorNorm(cell->error) > options.distanceTolerance * .125)
                cell = &source.Tighten(id, out);
        }
        Residual f;
        {
            Measure time(out.residualBuildMs);
            f = Build(*cell, frame, p, scale, products);
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
