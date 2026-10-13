// First-order forward chain rule, seeded in Increment's local coordinates.
// No perturbation size or numerical subtraction is used on regular branches.
#include "ApproximateInternal.h"
#include <chrono>
namespace bspline2analytic::detail
{
namespace
{
struct BranchFailure
{
};
template <size_t N> struct Jet
{
    double v = 0;
    std::array<double, N> d{};
    Jet() = default;
    Jet(double value) : v(value) {}
    Jet operator-() const
    {
        Jet r(-v);
        for (size_t j = 0; j < N; ++j)
            r.d[j] = -d[j];
        return r;
    }
    Jet operator+(const Jet &b) const
    {
        Jet r(v + b.v);
        for (size_t j = 0; j < N; ++j)
            r.d[j] = d[j] + b.d[j];
        return r;
    }
    Jet operator-(const Jet &b) const
    {
        return *this + (-b);
    }
    Jet operator*(const Jet &b) const
    {
        Jet r(v * b.v);
        for (size_t j = 0; j < N; ++j)
            r.d[j] = d[j] * b.v + v * b.d[j];
        return r;
    }
    Jet operator/(const Jet &b) const
    {
        if (std::abs(b.v) < 1e-30)
            throw BranchFailure{};
        Jet r(v / b.v);
        for (size_t j = 0; j < N; ++j)
            r.d[j] = (d[j] - r.v * b.d[j]) / b.v;
        return r;
    }
};
template <size_t N> Jet<N> Sin(Jet<N> a)
{
    Jet<N> r(std::sin(a.v));
    for (size_t j = 0; j < N; ++j)
        r.d[j] = std::cos(a.v) * a.d[j];
    return r;
}
template <size_t N> Jet<N> Cos(Jet<N> a)
{
    Jet<N> r(std::cos(a.v));
    for (size_t j = 0; j < N; ++j)
        r.d[j] = -std::sin(a.v) * a.d[j];
    return r;
}
template <size_t N> Jet<N> Sqrt(Jet<N> a)
{
    if (a.v < 1e-28)
        throw BranchFailure{};
    Jet<N> r(std::sqrt(a.v));
    for (size_t j = 0; j < N; ++j)
        r.d[j] = a.d[j] / (2 * r.v);
    return r;
}
template <size_t N> Jet<N> Atan2(Jet<N> y, Jet<N> x)
{
    double q = x.v * x.v + y.v * y.v;
    if (q < 1e-28)
        throw BranchFailure{};
    Jet<N> r(std::atan2(y.v, x.v));
    for (size_t j = 0; j < N; ++j)
        r.d[j] = (x.v * y.d[j] - y.v * x.d[j]) / q;
    return r;
}
template <size_t N> Jet<N> Unwrapped(Jet<N> a, double ref)
{
    a.v = Unwrap(a.v, ref);
    if (std::abs(std::abs(a.v - ref) - Pi) < 1e-9)
        throw BranchFailure{};
    return a;
}
template <size_t N> struct Vec
{
    Jet<N> x, y, z;
    Vec() = default;
    Vec(Jet<N> a, Jet<N> b, Jet<N> c) : x(a), y(b), z(c) {}
    explicit Vec(V a) : x(a.x), y(a.y), z(a.z) {}
    Vec operator+(const Vec &b) const
    {
        return {x + b.x, y + b.y, z + b.z};
    }
    Vec operator-(const Vec &b) const
    {
        return {x - b.x, y - b.y, z - b.z};
    }
    Vec operator*(Jet<N> b) const
    {
        return {x * b, y * b, z * b};
    }
    Vec operator/(Jet<N> b) const
    {
        return {x / b, y / b, z / b};
    }
    V value() const
    {
        return {x.v, y.v, z.v};
    }
    V derivative(size_t j) const
    {
        return {x.d[j], y.d[j], z.d[j]};
    }
};
template <size_t N> Jet<N> DotJ(Vec<N> a, Vec<N> b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <size_t N> Vec<N> CrossJ(Vec<N> a, Vec<N> b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <size_t N> Vec<N> UnitJ(Vec<N> a)
{
    return a / Sqrt(DotJ(a, a));
}

template <size_t N> struct Model
{
    using S = Jet<N>;
    using W = Vec<N>;
    sggk::SurfaceType type;
    W o, x, y, z;
    S radius, major, slope;
    std::array<S, 4> domain;
    struct Edge
    {
        S start, end, g0, g1, radius;
        W o, d, ex, ey;
        bool circle = false;
    };
    std::array<Edge, 4> edges;
    Model(const AnalyticPatch &p, const ActiveParameters &active, double scale)
        : type(p.parameters.type)
    {
        auto f = FrameOf(p.parameters);
        auto seed = [&](double value, unsigned index, double derivative) {
            S a(value);
            for (size_t j = 0; j < active.indices.size(); ++j)
                if (active.indices[j] == index)
                    a.d[j] = derivative;
            return a;
        };
        unsigned i = 0;
        W shift;
        if (type == sggk::SurfaceType::Plane)
            shift.z = seed(0, i++, scale);
        else if (type == sggk::SurfaceType::Cylinder)
        {
            shift.x = seed(0, i++, scale);
            shift.y = seed(0, i++, scale);
        }
        else
        {
            shift.x = seed(0, i++, scale);
            shift.y = seed(0, i++, scale);
            shift.z = seed(0, i++, scale);
        }
        W omega;
        omega.x = seed(0, i++, 1);
        omega.y = seed(0, i++, 1);
        if (type == sggk::SurfaceType::Plane)
            omega.z = seed(0, i++, 1);
        // Exact first derivative of Rodrigues at the zero local increment.
        W globalOmega = W(f.x) * omega.x + W(f.y) * omega.y + W(f.z) * omega.z;
        o = W(f.o) + W(f.x) * shift.x + W(f.y) * shift.y + W(f.z) * shift.z;
        x = W(f.x) + CrossJ(globalOmega, W(f.x));
        y = W(f.y) + CrossJ(globalOmega, W(f.y));
        z = W(f.z) + CrossJ(globalOmega, W(f.z));
        if (type == sggk::SurfaceType::Cone)
        {
            double t = std::tan(p.parameters.halfAngle);
            slope = seed(t, i++, t);
        }
        else if (type != sggk::SurfaceType::Plane)
            radius = seed(p.parameters.radius, i++, p.parameters.radius);
        if (type == sggk::SurfaceType::Torus)
            major = seed(p.parameters.majorRadius, i++, scale);
        for (unsigned j = 0; j < 4; ++j)
        {
            bool length = type == sggk::SurfaceType::Plane ||
                          (j >= 2 && (type == sggk::SurfaceType::Cylinder ||
                                      type == sggk::SurfaceType::Cone));
            domain[j] = seed(p.domain[j], i++, length ? scale : 1);
        }
        for (unsigned k = 0; k < 4; ++k)
        {
            auto &e = edges[k];
            unsigned dim = p.mapping.swappedUV ? 1 : 0;
            unsigned fixed = k < 2 ? dim : 1 - dim;
            bool rev = k < 2 ? p.mapping.reversedU : p.mapping.reversedV;
            S value = domain[2 * fixed + ((k % 2 == 1) != rev ? 1 : 0)];
            unsigned varying = 1 - fixed;
            bool revV = k < 2 ? p.mapping.reversedV : p.mapping.reversedU;
            e.start = domain[2 * varying + (revV ? 1 : 0)];
            e.end = domain[2 * varying + (revV ? 0 : 1)];
            S c = Cos(value), s = Sin(value);
            if (type == sggk::SurfaceType::Plane)
            {
                e.o = fixed == 0 ? W(value, 0, 0) : W(0, value, 0);
                e.d = fixed == 0 ? W(0, 1, 0) : W(1, 0, 0);
            }
            else if (fixed == 0 &&
                     (type == sggk::SurfaceType::Cylinder || type == sggk::SurfaceType::Cone))
            {
                if (type == sggk::SurfaceType::Cylinder)
                {
                    e.o = {radius * c, radius * s, 0};
                    e.d = {0, 0, 1};
                }
                else
                    e.d = {slope * c, slope * s, 1};
            }
            else
            {
                e.circle = true;
                if (fixed == 1)
                {
                    e.ex = {1, 0, 0};
                    e.ey = {0, 1, 0};
                    if (type == sggk::SurfaceType::Cylinder)
                    {
                        e.radius = radius;
                        e.o.z = value;
                    }
                    else if (type == sggk::SurfaceType::Cone)
                    {
                        e.radius = value * slope;
                        e.o.z = value;
                    }
                    else if (type == sggk::SurfaceType::Sphere)
                    {
                        e.radius = radius * s;
                        e.o.z = -radius * c;
                    }
                    else
                    {
                        e.radius = major + radius * c;
                        e.o.z = radius * s;
                    }
                }
                else if (type == sggk::SurfaceType::Sphere)
                {
                    e.radius = radius;
                    e.ex = {0, 0, -1};
                    e.ey = {c, s, 0};
                }
                else
                {
                    e.radius = radius;
                    e.o = {major * c, major * s, 0};
                    e.ex = {c, s, 0};
                    e.ey = {0, 0, 1};
                }
            }
        }
    }
    W local(V a) const
    {
        W d = W(a) - o;
        return {DotJ(d, x), DotJ(d, y), DotJ(d, z)};
    }
    W vector(W a) const
    {
        return x * a.x + y * a.y + z * a.z;
    }
    W global(W a) const
    {
        return o + vector(a);
    }
    S coordinate(const Edge &e, W a) const
    {
        a = a - e.o;
        return e.circle ? Atan2(DotJ(a, e.ey), DotJ(a, e.ex)) : DotJ(a, e.d) / DotJ(e.d, e.d);
    }
    W edgePoint(const Edge &e, S t) const
    {
        return e.circle ? e.o + (e.ex * Cos(t) + e.ey * Sin(t)) * e.radius : e.o + e.d * t;
    }
    W canonical(S u, S v) const
    {
        S c = Cos(u), s = Sin(u);
        switch (type)
        {
        case sggk::SurfaceType::Plane:
            return {u, v, 0};
        case sggk::SurfaceType::Cylinder:
            return {radius * c, radius * s, v};
        case sggk::SurfaceType::Cone:
            return {v * slope * c, v * slope * s, v};
        case sggk::SurfaceType::Sphere:
            return {radius * Sin(v) * c, radius * Sin(v) * s, -radius * Cos(v)};
        default:
            return {(major + radius * Cos(v)) * c, (major + radius * Cos(v)) * s, radius * Sin(v)};
        }
    }
    W project(W a) const
    {
        S u = a.x, v = a.y;
        if (type != sggk::SurfaceType::Plane)
        {
            u = Atan2(a.y, a.x);
            S rho = Sqrt(a.x * a.x + a.y * a.y);
            if (type == sggk::SurfaceType::Cylinder)
                v = a.z;
            else if (type == sggk::SurfaceType::Cone)
            {
                v = (rho * slope + a.z) / (S(1) + slope * slope);
                if (std::abs(v.v) < 1e-10)
                    throw BranchFailure{};
                if (v.v < 0)
                    v = 0;
            }
            else if (type == sggk::SurfaceType::Sphere)
                v = Atan2(rho, -a.z);
            else
                v = Atan2(a.z, rho - major);
            u = Unwrapped(u, (domain[0].v + domain[1].v) * .5);
            if (type == sggk::SurfaceType::Torus)
                v = Unwrapped(v, (domain[2].v + domain[3].v) * .5);
        }
        // Projection changes branch on the finite patch boundary.
        for (unsigned j = 0; j < 4; ++j)
            if (std::abs((j < 2 ? u.v : v.v) - domain[j].v) < 1e-10)
                throw BranchFailure{};
        if (u.v >= domain[0].v && u.v <= domain[1].v && v.v >= domain[2].v && v.v <= domain[3].v)
            return canonical(u, v);
        struct Candidate
        {
            W point;
            double distance;
            bool unstable;
        };
        std::vector<Candidate> candidates;
        candidates.reserve(12);
        for (const auto &e : edges)
        {
            S t = coordinate(e, a);
            bool unstable = false;
            if (e.circle)
            {
                double ref = (e.start.v + e.end.v) * .5;
                t.v = Unwrap(t.v, ref);
                unstable = std::abs(std::abs(t.v - ref) - Pi) < 1e-9;
            }
            S lo = e.start.v < e.end.v ? e.start : e.end,
              hi = e.start.v < e.end.v ? e.end : e.start;
            unstable |= std::abs(t.v - lo.v) < 1e-10 || std::abs(t.v - hi.v) < 1e-10;
            if (t.v < lo.v)
                t = lo;
            else if (t.v > hi.v)
                t = hi;
            auto add = [&](S value, bool bad) {
                W q = edgePoint(e, value);
                V d = q.value() - a.value();
                candidates.push_back({q, Dot(d, d), bad});
            };
            add(t, unstable);
            add(e.start, false);
            add(e.end, false);
        }
        const auto best =
            std::min_element(candidates.begin(), candidates.end(),
                             [](const auto &a, const auto &b) { return a.distance < b.distance; });
        if (best->unstable)
            throw BranchFailure{};
        for (const auto &q : candidates)
            if (std::abs(q.distance - best->distance) < 1e-12 * std::max(1.0, best->distance) &&
                Norm(q.point.value() - best->point.value()) > 1e-9)
                throw BranchFailure{};
        return best->point;
    }
};

// A block contains one distance row and, on a boundary, one tangent row.
static std::vector<size_t> BlockRows(const SampleCache &samples)
{
    std::vector<size_t> rows(samples.size());
    size_t row = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        rows[i] = row;
        row += samples[i].edge < 0 ? 1 : 2;
    }
    return rows;
}
template <size_t N>
std::vector<bool> Analytic(const AnalyticPatch &p, const SampleCache &samples,
                           const ConversionOptions &options, const ActiveParameters &active,
                           double scale, EvaluationWorkspace &w, RefinementDiagnostics &diag)
{
    Model<N> m(p, active, scale);
    std::array<bool, 4> badEdge{};
    for (unsigned k = 0; k < 4; ++k)
    {
        auto &e = m.edges[k];
        const auto &order = samples.boundaryOrder[k];
        try
        {
            if (order.empty())
                throw BranchFailure{};
            double last = e.start.v;
            for (size_t j = 0; j < order.size(); ++j)
            {
                auto g = m.coordinate(e, m.local(samples[order[j]].point));
                if (e.circle)
                    g = Unwrapped(g, last);
                if (!j)
                    e.g0 = g;
                last = g.v;
                e.g1 = g;
            }
            if (std::abs(e.g1.v - e.g0.v) < 1e-10)
                throw BranchFailure{};
        }
        catch (const BranchFailure &)
        {
            badEdge[k] = true;
        }
    }
    const double eps = options.validation.distanceTolerance;
    const double angle = options.validation.boundaryTangentToleranceDeg * Pi / 180,
                 chord = 2 * std::sin(angle / 2);
    auto rows = BlockRows(samples);
    std::vector<bool> fallback(samples.size(), false);
    for (size_t i = 0; i < samples.size(); ++i)
    {
        const auto &s = samples[i];
        size_t row = rows[i];
        try
        {
            auto a = m.local(s.point);
            Vec<N> q, t;
            if (s.edge < 0)
                q = m.global(m.project(a));
            else
            {
                if (badEdge[s.edge])
                    throw BranchFailure{};
                const auto &e = m.edges[s.edge];
                auto g = m.coordinate(e, a);
                if (e.circle)
                    g = Unwrapped(g, e.g0.v + (e.g1.v - e.g0.v) * s.fraction);
                auto phi = e.start + (e.end - e.start) * (g - e.g0) / (e.g1 - e.g0);
                q = m.global(m.edgePoint(e, phi));
                t = UnitJ(m.vector(e.circle ? e.ex * (-Sin(phi)) + e.ey * Cos(phi) : e.d)) *
                    Jet<N>(e.end.v > e.start.v ? 1 : -1);
            }
            V r{w.current.leastSquares[3 * row], w.current.leastSquares[3 * row + 1],
                w.current.leastSquares[3 * row + 2]};
            double distance = Norm(r);
            bool zeroDistance = distance < 1e-8;
            if (zeroDistance)
                ++diag.nonsmoothConstraintRows;
            V cross{};
            double sine = 0, cosine = 0;
            if (s.edge >= 0)
            {
                cross = Cross(s.tangent, t.value());
                sine = Norm(cross);
                cosine = Dot(s.tangent, t.value());
                if (sine < 1e-12)
                {
                    if (cosine < 0)
                        throw BranchFailure{};
                    ++diag.nonsmoothConstraintRows;
                }
            }
            for (size_t j = 0; j < active.indices.size(); ++j)
            {
                V dr = q.derivative(j) * (-1 / eps);
                w.jacobian[3 * row][j] = dr.x;
                w.jacobian[3 * row + 1][j] = dr.y;
                w.jacobian[3 * row + 2][j] = dr.z;
                w.constraintJacobian[row][j] = zeroDistance ? 0 : Dot(r, dr) / distance;
                if (s.edge >= 0)
                {
                    V dt = t.derivative(j), rt = dt * (-1 / chord);
                    w.jacobian[3 * row + 3][j] = rt.x;
                    w.jacobian[3 * row + 4][j] = rt.y;
                    w.jacobian[3 * row + 5][j] = rt.z;
                    double ds = sine < 1e-12 ? 0 : Dot(cross, Cross(s.tangent, dt)) / sine;
                    w.constraintJacobian[row + 1][j] =
                        sine < 1e-12 ? 0
                                     : (cosine * ds - sine * Dot(s.tangent, dt)) /
                                           ((sine * sine + cosine * cosine) * angle);
                }
            }
            ++diag.analyticSampleBlocks;
        }
        catch (const BranchFailure &)
        {
            fallback[i] = true;
            ++diag.fallbackSampleBlocks;
        }
    }
    return fallback;
}
// For these size-only levels the chosen finite-patch point scales or
// translates exactly, and every boundary tangent is unchanged. Recover Q
// from the current distance residual; trial geometry is never consulted.
static bool SimpleSizeJacobian(const AnalyticPatch &p, const SampleCache &samples,
                               const ConversionOptions &o, const ActiveParameters &active,
                               double scale, EvaluationWorkspace &w, RefinementDiagnostics &diag)
{
    using T = sggk::SurfaceType;
    auto type = p.parameters.type;
    unsigned index = type == T::Plane ? 0 : type == T::Cylinder ? 4 : 5;
    if (active.indices.size() != 1 || active.indices[0] != index ||
        (type != T::Plane && type != T::Cylinder && type != T::Sphere))
        return false;
    auto f = FrameOf(p.parameters);
    double eps = o.validation.distanceTolerance;
    size_t row = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        V r{w.current.leastSquares[3 * row], w.current.leastSquares[3 * row + 1],
            w.current.leastSquares[3 * row + 2]};
        V q = samples[i].point - r * eps;
        V dq = type == T::Plane ? f.z * scale : q - f.o;
        if (type == T::Cylinder)
            dq = dq - f.z * Dot(dq, f.z);
        V dr = dq * (-1 / eps);
        w.jacobian[3 * row][0] = dr.x;
        w.jacobian[3 * row + 1][0] = dr.y;
        w.jacobian[3 * row + 2][0] = dr.z;
        double distance = Norm(r);
        if (distance < 1e-8)
            ++diag.nonsmoothConstraintRows;
        else
            w.constraintJacobian[row][0] = Dot(r, dr) / distance;
        ++row;
        if (samples[i].edge >= 0)
        {
            if (w.current.constraints[row] < 1e-9)
                ++diag.nonsmoothConstraintRows;
            ++row;
        }
        ++diag.analyticSampleBlocks;
    }
    return true;
}
static bool Differences(const AnalyticPatch &p, const SampleCache &samples,
                        const ConversionOptions &o, const ActiveParameters &active, double scale,
                        EvaluationWorkspace &w, const std::vector<bool> &fallback,
                        RefinementDiagnostics &diag)
{
    bool full = o.refinement.jacobianMode == JacobianMode::FiniteDifference;
    SampleCache subset;
    std::vector<size_t> original;
    // Keep the complete ordered boundary walks: endpoints alone lose winding.
    for (size_t i = 0; i < samples.size(); ++i)
        if (!full && (fallback[i] || samples[i].edge >= 0))
        {
            subset.points.push_back(samples[i]);
            original.push_back(i);
        }
    if (!full)
        subset.RebuildIndices();
    const auto &input = full ? samples : subset;
    auto sourceRows = BlockRows(samples), localRows = BlockRows(input);
    EvaluationWorkspace scratch;
    auto &evaluation = full ? w : scratch;
    ResidualBuffer base;
    if (!full)
    {
        EvaluateResiduals(p, input, o.validation, scratch, base);
        ++diag.fallbackResidualEvaluations;
    }
    const auto &current = full ? w.current : base;
    std::vector<double> delta(
        ParametersFor(p.parameters.type, ParameterLevel::FullPatch).indices.size(), 0);
    for (size_t j = 0; j < active.indices.size(); ++j)
    {
        constexpr double h = 2e-6;
        bool plus = false, minus = false;
        const unsigned before = evaluation.evaluations;
        delta[active.indices[j]] = h;
        try
        {
            EvaluateResiduals(Increment(p, delta, scale), input, o.validation, evaluation,
                              evaluation.plus);
            plus = true;
        }
        catch (const Failure &)
        {
        }
        delta[active.indices[j]] = -h;
        try
        {
            EvaluateResiduals(Increment(p, delta, scale), input, o.validation, evaluation,
                              evaluation.minus);
            minus = true;
        }
        catch (const Failure &)
        {
        }
        delta[active.indices[j]] = 0;
        if (!full)
            diag.fallbackResidualEvaluations += evaluation.evaluations - before;
        if (!plus && !minus)
        {
            w.evaluations += scratch.evaluations;
            return false;
        }
        const auto &a = plus ? evaluation.plus : current;
        const auto &b = minus ? evaluation.minus : current;
        double denominator = plus && minus ? 2 * h : h;
        for (size_t i = 0; i < input.size(); ++i)
        {
            size_t orig = full ? i : original[i];
            if (!full && !fallback[orig])
                continue;
            size_t sr = sourceRows[orig], lr = localRows[i], count = input[i].edge < 0 ? 1 : 2;
            for (size_t k = 0; k < count; ++k)
            {
                w.constraintJacobian[sr + k][j] =
                    (a.constraints[lr + k] - b.constraints[lr + k]) / denominator;
                for (unsigned c = 0; c < 3; ++c)
                    w.jacobian[3 * (sr + k) + c][j] =
                        (a.leastSquares[3 * (lr + k) + c] - b.leastSquares[3 * (lr + k) + c]) /
                        denominator;
            }
        }
    }
    w.evaluations += scratch.evaluations;
    return true;
}
} // namespace
bool BuildJacobian(const AnalyticPatch &p, const SampleCache &samples, const ConversionOptions &o,
                   const ActiveParameters &active, double scale, EvaluationWorkspace &w,
                   RefinementDiagnostics &diag)
{
    const auto begin = std::chrono::steady_clock::now();
    auto resize = [&](auto &matrix, size_t rows) {
        bool changed = matrix.size() != rows;
        matrix.resize(rows);
        for (auto &row : matrix)
        {
            changed |= row.size() != active.indices.size();
            row.assign(active.indices.size(), 0);
        }
        if (changed)
            ++w.resizes;
    };
    resize(w.jacobian, w.current.leastSquares.size());
    resize(w.constraintJacobian, w.current.constraints.size());
    std::vector<bool> fallback(samples.size(), false);
    bool ok = true;
    if (o.refinement.jacobianMode == JacobianMode::AnalyticWithFallback &&
        !SimpleSizeJacobian(p, samples, o, active, scale, w, diag))
    {
        // Specialize the frequent one/two-variable levels to avoid unused lanes.
        if (active.indices.size() == 1)
            fallback = Analytic<1>(p, samples, o, active, scale, w, diag);
        else if (active.indices.size() == 2)
            fallback = Analytic<2>(p, samples, o, active, scale, w, diag);
        else
            fallback = Analytic<11>(p, samples, o, active, scale, w, diag);
    }
    if (o.refinement.jacobianMode == JacobianMode::FiniteDifference ||
        std::any_of(fallback.begin(), fallback.end(), [](bool b) { return b; }))
        ok = Differences(p, samples, o, active, scale, w, fallback, diag);
    for (const auto *matrix : {&w.jacobian, &w.constraintJacobian})
        for (const auto &row : *matrix)
            for (double v : row)
                ok &= std::isfinite(v);
    diag.jacobianMs +=
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - begin).count();
    return ok;
}
} // namespace bspline2analytic::detail
