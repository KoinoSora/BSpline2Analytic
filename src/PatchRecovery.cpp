#include "PatchRecovery.h"
#include "Numerics.h"

namespace bspline2analytic::detail
{
void RecoverPatch(const sggk::BSplineSurface &input, const Candidate &c,
                  const std::array<BoundaryInfo, 4> &edges, const RecognitionOptions &o,
                  sggk::UVRange &range, ParameterCorrespondence &mapping)
{
    using UV = std::array<double, 2>;
    std::array<std::vector<UV>, 4> paths;
    bool periodU = c.parameters.type != sggk::SurfaceType::Plane;
    bool periodV = c.parameters.type == sggk::SurfaceType::Torus &&
                   c.parameters.majorRadius > c.parameters.radius;
    auto eval = [&](const sggk::Point3D &p) {
        auto q = c.surface->CalcParam(p);
        UV r{q.X(), q.Y()};
        if (!std::isfinite(r[0]) || !std::isfinite(r[1]))
            throw Failure(Status::VerificationFailed, "Nonfinite analytic parameter");
        return r;
    };
    auto trace = [&](unsigned index, UV start) {
        auto bc = std::dynamic_pointer_cast<sggk::BSplineCurve3D>(edges[index].curve);
        UV last = start;
        for (double t : SampleParameters(*bc, std::max(8u, o.samplesPerSpan)))
        {
            UV q = eval(bc->CalcPoint(t));
            if (periodU)
                q[0] = Unwrap(q[0], last[0]);
            if (periodV)
                q[1] = Unwrap(q[1], last[1]);
            paths[index].push_back(q);
            last = q;
        }
    };
    UV corner = eval(input.CalcPoint(input.MinParamU(), input.MinParamV()));
    trace(2, corner);
    trace(0, corner);
    trace(1, paths[2].back());
    trace(3, paths[0].back());
    double angularScale =
        std::max({c.parameters.radius, std::abs(c.parameters.majorRadius), 1e-12});
    std::array<double, 2> metric{periodU ? angularScale : 1.0,
                                 c.parameters.type == sggk::SurfaceType::Sphere ||
                                         c.parameters.type == sggk::SurfaceType::Torus
                                     ? std::max(c.parameters.radius, 1e-12)
                                     : 1.0};
    std::array<unsigned, 4> fixed{};
    UV lo{INFINITY, INFINITY}, hi{-INFINITY, -INFINITY};
    for (unsigned i = 0; i < 4; ++i)
    {
        UV low = paths[i].front(), high = low;
        for (const auto &q : paths[i])
            for (unsigned j = 0; j < 2; ++j)
            {
                low[j] = std::min(low[j], q[j]);
                high[j] = std::max(high[j], q[j]);
                lo[j] = std::min(lo[j], q[j]);
                hi[j] = std::max(hi[j], q[j]);
            }
        fixed[i] = (high[0] - low[0]) * metric[0] < (high[1] - low[1]) * metric[1] ? 0 : 1;
        if ((high[fixed[i]] - low[fixed[i]]) * metric[fixed[i]] > 8 * o.distanceTolerance)
            throw Failure(Status::InconsistentBoundary,
                          "A boundary is not an isoparametric curve in the recovered frame");
        unsigned varying = 1 - fixed[i];
        double direction = 0;
        for (std::size_t j = 1; j < paths[i].size(); ++j)
        {
            double delta = paths[i][j][varying] - paths[i][j - 1][varying];
            if (std::abs(delta) * metric[varying] > o.distanceTolerance)
            {
                if (direction * delta < 0)
                    throw Failure(Status::InconsistentBoundary,
                                  "Analytic boundary parameter retraces");
                direction = std::copysign(1.0, delta);
            }
        }
    }
    if (fixed[0] != fixed[1] || fixed[2] != fixed[3] || fixed[0] == fixed[2])
        throw Failure(Status::InconsistentBoundary,
                      "Opposite boundaries do not define an analytic rectangle");
    for (unsigned j = 0; j < 2; ++j)
    {
        if ((hi[j] - lo[j]) * metric[j] <= Roundoff(metric[j]))
            throw Failure(Status::Degenerate, "Analytic patch has zero extent");
        if (std::abs(paths[1].back()[j] - paths[3].back()[j]) * metric[j] > 8 * o.distanceTolerance)
            throw Failure(Status::InconsistentBoundary,
                          "Analytic boundary loop does not close consistently");
    }
    if ((periodU && hi[0] - lo[0] > 2 * Pi + o.angularTolerance) ||
        (periodV && hi[1] - lo[1] > 2 * Pi + o.angularTolerance))
        throw Failure(Status::Unsupported, "Multiple-turn patches are not supported");
    mapping.swappedUV = fixed[2] == 0;
    mapping.reversedU = paths[2].back()[1 - fixed[2]] < paths[2].front()[1 - fixed[2]];
    mapping.reversedV = paths[0].back()[1 - fixed[0]] < paths[0].front()[1 - fixed[0]];
    std::array<UV, 4> corners{paths[2].front(), paths[2].back(), paths[0].back(), paths[1].back()};
    for (unsigned i = 0; i < 4; ++i)
        mapping.analyticCorners[i] = sggk::Point2D(corners[i][0], corners[i][1]);
    range = sggk::UVRange(sggk::Interval(lo[0], hi[0]), sggk::Interval(lo[1], hi[1]));
}
} // namespace bspline2analytic::detail
