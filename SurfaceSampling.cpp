#include "ApproximateInternal.h"
namespace bspline2analytic::detail
{
static Sample At(const sggk::BSplineSurface &s, double u, double v, int edge,
                 sggk::DerivType du = sggk::DerivType::NONE,
                 sggk::DerivType dv = sggk::DerivType::NONE)
{
    Sample x;
    x.u = u;
    x.v = v;
    x.edge = edge;
    if (edge < 0)
    {
        x.point = Get(s.CalcPoint(u, v));
        return x;
    }
    if (du == sggk::DerivType::NONE)
        du = u == s.MaxParamU() ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
    if (dv == sggk::DerivType::NONE)
        dv = v == s.MaxParamV() ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
    sggk::Point3D p;
    auto d = s.CalcPartialD1(u, v, du, dv, p);
    x.point = Get(p);
    if (edge >= 0)
    {
        x.tangent = Unit(Get(d[edge < 2 ? 1 : 0]));
        x.normal = Unit(Cross(Get(d[0]), Get(d[1])));
        x.fraction = edge < 2 ? (v - s.MinParamV()) / (s.MaxParamV() - s.MinParamV())
                              : (u - s.MinParamU()) / (s.MaxParamU() - s.MinParamU());
    }
    return x;
}
static std::vector<std::pair<double, double>> Spans(const sggk::RealArray &k, double lo, double hi)
{
    std::vector<double> ends{lo};
    for (double t : k)
        if (t > lo && t < hi)
            ends.push_back(t);
    ends.push_back(hi);
    std::vector<std::pair<double, double>> r;
    for (size_t i = 1; i < ends.size(); ++i)
        if (ends[i] > ends[i - 1])
            r.push_back({ends[i - 1], ends[i]});
    return r;
}
SampleCache SampleSurface(const sggk::BSplineSurface &s, const SamplingOptions &o)
{
    if (!o.interiorPointsPerSpan || o.interiorPointsPerSpan > 8 || o.boundaryPointsPerSpan < 3 ||
        o.boundaryPointsPerSpan > 32 || !o.maxSamples)
        throw Failure(Status::InvalidInput, "Invalid sampling options");
    SampleCache cache;
    auto us = Spans(s.KnotsU(), s.MinParamU(), s.MaxParamU()),
         vs = Spans(s.KnotsV(), s.MinParamV(), s.MaxParamV());
    auto add = [&](Sample x) {
        if (cache.points.size() >= o.maxSamples)
            throw Failure(Status::VerificationFailed, "Sample budget exhausted");
        cache.points.push_back(x);
    };
    for (auto u : us)
        for (auto v : vs)
        {
            for (unsigned i = 0; i < o.interiorPointsPerSpan; ++i)
                for (unsigned j = 0; j < o.interiorPointsPerSpan; ++j)
                {
                    double x = o.interiorPointsPerSpan == 2
                                   ? (i ? 0.7886751345948129 : 0.2113248654051871)
                                   : (i + 0.5) / o.interiorPointsPerSpan;
                    double y = o.interiorPointsPerSpan == 2
                                   ? (j ? 0.7886751345948129 : 0.2113248654051871)
                                   : (j + 0.5) / o.interiorPointsPerSpan;
                    add(At(s, u.first + (u.second - u.first) * x,
                           v.first + (v.second - v.first) * y, -1));
                }
            add(At(s, (u.first + u.second) * 0.5, (v.first + v.second) * 0.5, -1));
        }
    for (int edge = 0; edge < 4; ++edge)
    {
        for (auto span : edge < 2 ? vs : us)
            for (unsigned j = 0; j < o.boundaryPointsPerSpan; ++j)
            {
                double t =
                    span.first + (span.second - span.first) * j / (o.boundaryPointsPerSpan - 1);
                auto side = j == 0 ? sggk::DerivType::RIGHT
                                   : (j + 1 == o.boundaryPointsPerSpan ? sggk::DerivType::LEFT
                                                                       : sggk::DerivType::NONE);
                auto fixed = edge % 2 ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
                add(At(s, edge < 2 ? (edge == 0 ? s.MinParamU() : s.MaxParamU()) : t,
                       edge < 2 ? t : (edge == 2 ? s.MinParamV() : s.MaxParamV()), edge,
                       edge < 2 ? fixed : side, edge < 2 ? side : fixed));
            }
    }
    return cache;
}
void AddWitnessSamples(const sggk::BSplineSurface &s, const ConversionValidationReport &r,
                       const SamplingOptions &o, SampleCache &cache)
{
    auto add = [&](const ConstraintReport &c, bool angle) {
        if (c.state == CheckState::Passed || c.state == CheckState::NotChecked ||
            cache.points.size() >= o.maxSamples)
            return;
        double u = (c.parameterRegion[0] + c.parameterRegion[1]) * 0.5,
               v = (c.parameterRegion[2] + c.parameterRegion[3]) * 0.5;
        int edge = angle ? static_cast<int>(c.boundary) : -1;
        if (!angle)
        {
            if (u == s.MinParamU())
                edge = 0;
            else if (u == s.MaxParamU())
                edge = 1;
            else if (v == s.MinParamV())
                edge = 2;
            else if (v == s.MaxParamV())
                edge = 3;
        }
        if (edge >= 0)
        {
            if (edge < 2)
                u = edge == 0 ? s.MinParamU() : s.MaxParamU();
            else
                v = edge == 2 ? s.MinParamV() : s.MaxParamV();
        }
        if (u < s.MinParamU() || u > s.MaxParamU() || v < s.MinParamV() || v > s.MaxParamV())
            return;
        cache.points.push_back(At(s, u, v, edge));
        // Include immediate neighbours: a witness should influence a region,
        // rather than just interpolate one troublesome point.
        double du = (c.parameterRegion[1] - c.parameterRegion[0]) * 0.25,
               dv = (c.parameterRegion[3] - c.parameterRegion[2]) * 0.25;
        if (cache.points.size() + 2 <= o.maxSamples)
        {
            cache.points.push_back(At(s, u - du, v - dv, edge));
            cache.points.push_back(At(s, u + du, v + dv, edge));
        }
    };
    add(r.geometry, false);
    add(r.boundaryTangent, true);
    add(r.boundaryNormal, true);
}
} // namespace bspline2analytic::detail
