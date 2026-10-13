#include "ApproximateInternal.h"
namespace bspline2analytic::detail
{
void EvaluateResiduals(const AnalyticPatch &p, const SampleCache &samples,
                       const ValidationOptions &options, EvaluationWorkspace &w, ResidualBuffer &r)
{
    ++w.evaluations;
    auto resize = [&](auto &buffer, size_t size) {
        if (buffer.size() != size)
        {
            buffer.resize(size);
            ++w.resizes;
        }
    };
    size_t count = samples.size(), boundaryCount = 0;
    for (const auto &order : samples.boundaryOrder)
        boundaryCount += order.size();
    resize(w.local, count);
    resize(w.phase, count);
    resize(w.mappedPoints, count);
    resize(w.mappedTangents, count);
    resize(r.leastSquares, 3 * (count + boundaryCount));
    resize(r.constraints, count + boundaryCount);
    auto &context = w.context;
    context.frame = FrameOf(p.parameters);
    context.coneTangent =
        p.parameters.type == sggk::SurfaceType::Cone ? std::tan(p.parameters.halfAngle) : 0;
    for (unsigned k = 0; k < 4; ++k)
        context.edges[k] = MakeEdge(p, k);
    for (size_t i = 0; i < count; ++i)
    {
        const auto &s = samples[i];
        w.local[i] = context.frame.local(s.point);
        if (s.edge >= 0)
            w.phase[i] = EdgeCoordinateRaw(context.edges[s.edge], w.local[i]);
    }
    for (unsigned k = 0; k < 4; ++k)
    {
        auto &e = context.edges[k];
        const auto &order = samples.boundaryOrder[k];
        if (order.empty())
            throw Failure(Status::Degenerate, "Missing boundary samples");
        double last = e.start;
        for (size_t j = 0; j < order.size(); ++j)
        {
            double raw = w.phase[order[j]];
            last = e.circle ? Unwrap(raw, last) : raw;
            if (!j)
                e.g0 = last;
        }
        e.g1 = last;
        if (std::abs(e.g1 - e.g0) < 1e-12)
            throw Failure(Status::Degenerate, "Collapsed optimization correspondence");
    }
    r.cost = r.worst = 0;
    size_t row = 0;
    auto add = [&](V a, double metric) {
        r.leastSquares[3 * row] = a.x;
        r.leastSquares[3 * row + 1] = a.y;
        r.leastSquares[3 * row + 2] = a.z;
        r.constraints[row++] = metric;
        r.cost += Dot(a, a);
        r.worst = std::max(r.worst, metric);
    };
    double angleTolerance = options.boundaryTangentToleranceDeg * Pi / 180;
    double tangentScale = 2 * std::sin(angleTolerance / 2);
    for (size_t i = 0; i < count; ++i)
    {
        const auto &s = samples[i];
        if (s.edge < 0)
        {
            w.mappedPoints[i] = context.frame.global(
                ProjectInteriorPointLocal(p, w.local[i], context.edges, context.coneTangent));
            w.mappedTangents[i] = {};
        }
        else
        {
            const auto &e = context.edges[s.edge];
            // Preserve the original fraction-based unwrap branch. The sequential
            // boundary walk above supplies endpoints, not a new correspondence.
            double g =
                e.circle ? Unwrap(w.phase[i], e.g0 + (e.g1 - e.g0) * s.fraction) : w.phase[i];
            auto q = MapBoundaryPointAndTangent(context.frame, e, g);
            w.mappedPoints[i] = q.point;
            w.mappedTangents[i] = q.tangent;
        }
        V d = (s.point - w.mappedPoints[i]) / options.distanceTolerance;
        add(d, Norm(d));
        if (s.edge >= 0)
            add((s.tangent - w.mappedTangents[i]) / tangentScale,
                Angle(s.tangent, w.mappedTangents[i]) / angleTolerance);
    }
    if (!std::isfinite(r.cost) || !std::isfinite(r.worst))
        throw Failure(Status::VerificationFailed, "Nonfinite trial residual");
}
} // namespace bspline2analytic::detail
