#include "ApproximateCases.h"
#include "ApproximateInternal.h"
namespace bspline2analytic::tests
{
sggk::BSplineSurfacePtr PlateauFixture(const sggk::BSplineSurface &source)
{
    auto result = RefitPolynomial(source, 12);
    auto exact = Recognize(source);
    auto z = detail::Get(exact.parameters.axis);
    const auto cp = result->ControlPoints();
    for (unsigned i = 4; i + 4 < cp.size(); ++i)
        for (unsigned j = 4; j + 4 < cp[i].size(); ++j)
            result->SetControlPoint(detail::Point(detail::Get(cp[i][j]) + z * 1.6e-5), i, j);
    return result;
}
PerturbationCase CurvedPlateauFixture(const ReferenceCase &reference)
{
    using namespace detail;
    auto exact = Recognize(*reference.spline);
    if (exact.status != Status::Success)
        throw std::runtime_error("Perturbation reference recognition failed");
    PerturbationCase fixture;
    fixture.name = reference.name;
    fixture.spline = exact.parameters.type == sggk::SurfaceType::Plane
                         ? PlateauFixture(*reference.spline)
                         : RefitPolynomial(*reference.spline, 16);
    fixture.initial = RecoverApproxPatch(*fixture.spline, exact.parameters);
    auto f = FrameOf(exact.parameters);
    if (exact.parameters.type != sggk::SurfaceType::Plane)
    {
        const auto cp = fixture.spline->ControlPoints();
        // Positive B-spline basis coefficients make a broad plateau with a
        // transition zone. Four untouched outer rows preserve boundary jets.
        for (unsigned i = 4; i + 4 < cp.size(); ++i)
            for (unsigned j = 4; j + 4 < cp[i].size(); ++j)
            {
                V a = f.local(Get(cp[i][j]));
                auto uv = MotherParameters(exact.parameters, a);
                V displacement;
                if (exact.parameters.type == sggk::SurfaceType::Cone)
                    displacement = V{std::cos(uv[0]), std::sin(uv[0]), 0} *
                                   (1.6e-5 * a.z / fixture.initial.domain[3]);
                else
                    displacement = CanonicalNormal(exact.parameters, uv[0], uv[1]) * 1.6e-5;
                fixture.spline->SetControlPoint(Point(Get(cp[i][j]) + f.vector(displacement)), i,
                                                j);
            }
    }
    fixture.feasible = fixture.initial;
    auto &a = fixture.feasible.parameters;
    if (a.type == sggk::SurfaceType::Plane)
        a.origin = Point(Get(a.origin) + f.z * 8e-6);
    else if (a.type == sggk::SurfaceType::Cone)
        a.halfAngle = std::atan(std::tan(a.halfAngle) + 8e-6 / fixture.initial.domain[3]);
    else
        a.radius += 8e-6;
    MaterializePatch(fixture.feasible);
    // Feasibility is intentionally asserted by the independent validator in
    // tests, not assumed from the control-point displacement magnitude.
    return fixture;
}
sggk::BSplineSurfacePtr RefitPolynomial(const sggk::BSplineSurface &source, unsigned count)
{
    using namespace detail;
    auto grid = [&](const sggk::RealArray &knots, double lo, double hi) {
        std::vector<double> r{lo};
        double previous = lo;
        for (double k : knots)
        {
            if (k <= lo || k > hi)
                continue;
            for (unsigned j = 1; j <= count; ++j)
                r.push_back(previous + (k - previous) * j / count);
            previous = k;
        }
        if (previous < hi)
            for (unsigned j = 1; j <= count; ++j)
                r.push_back(previous + (hi - previous) * j / count);
        return r;
    };
    auto us = grid(source.KnotsU(), source.MinParamU(), source.MaxParamU()),
         vs = grid(source.KnotsV(), source.MinParamV(), source.MaxParamV());
    // Tensor bicubic interpolation of actual rational source points.
    // Using values avoids importing the wrong one-sided derivative at a
    // rational span junction. No source weights are copied or deleted.
    auto coefficients = [](std::array<V, 4> f) {
        return std::array<V, 4>{f[0], f[1] * 3 - f[2] * 1.5 - f[0] * (5.0 / 6) + f[3] / 3,
                                f[2] * 3 - f[1] * 1.5 + f[0] / 3 - f[3] * (5.0 / 6), f[3]};
    };
    sggk::Point3DMatrix cp(3 * (us.size() - 1) + 1,
                           std::vector<sggk::Point3D>(3 * (vs.size() - 1) + 1));
    for (size_t i = 0; i + 1 < us.size(); ++i)
        for (size_t j = 0; j + 1 < vs.size(); ++j)
        {
            std::array<std::array<V, 4>, 4> values;
            for (unsigned a = 0; a < 4; ++a)
            {
                for (unsigned b = 0; b < 4; ++b)
                {
                    double u = a == 3 ? us[i + 1] : us[i] + (us[i + 1] - us[i]) * a / 3;
                    double v = b == 3 ? vs[j + 1] : vs[j] + (vs[j + 1] - vs[j]) * b / 3;
                    values[a][b] = Get(source.CalcPoint(u, v));
                }
                values[a] = coefficients(values[a]);
            }
            for (unsigned b = 0; b < 4; ++b)
            {
                auto c = coefficients({values[0][b], values[1][b], values[2][b], values[3][b]});
                for (unsigned a = 0; a < 4; ++a)
                    cp[3 * i + a][3 * j + b] = Point(c[a]);
            }
        }
    sggk::UIntArray mu(us.size(), 3), mv(vs.size(), 3);
    mu.front() = mu.back() = mv.front() = mv.back() = 4;
    return std::make_shared<sggk::BSplineSurface>(3, 3, cp, us, vs, mu, mv);
}
} // namespace bspline2analytic::tests
