/**
 * @file SurfaceSampleContext.cpp
 * @brief 单次调用内的惰性源面求值与采样分区。
 * 按参数及导数侧别缓存点，生成增强 Greville 与独立评分视图，并向下传播共享截止时间。
 */
#include "SurfaceSampleContext.h"
#include <set>
namespace bspline2analytic::detail
{
static std::vector<double> GrevilleParameters(const sggk::RealArray &knots,
    const sggk::UIntArray &mults, unsigned degree, size_t count, double lo, double hi)
{
    if (!degree || knots.size() != mults.size())
        throw Failure(Status::InvalidInput, "Invalid Greville knot data");
    std::vector<double> full;
    for (size_t i = 0; i < knots.size(); ++i)
    {
        if (!std::isfinite(knots[i]) || !mults[i] || mults[i] > degree + 1 ||
            (i && knots[i] <= knots[i - 1]))
            throw Failure(Status::InvalidInput, "Invalid Greville knots/multiplicities");
        full.insert(full.end(), mults[i], knots[i]);
    }
    if (full.size() != count + degree + 1)
        throw Failure(Status::InvalidInput, "Invalid full knot count");
    std::set<double> out{lo, hi};
    for (size_t i = 0; i < count; ++i)
    {
        double t = 0;
        for (unsigned k = 1; k <= degree; ++k)
            t += full[i + k] / degree;
        out.insert(std::clamp(t, lo, hi));
    }
    // 每个非零结点区间都包含内部见证点，一次输入也不例外。
    for (size_t i = 1; i < knots.size(); ++i)
    {
        double a = std::max(lo, knots[i - 1]), b = std::min(hi, knots[i]);
        if (a < b)
            out.insert(a + (b - a) * .5);
    }
    // 保证内部点足够构造四点圆环初值，并保留独立评分集。
    for (unsigned i = 1; i < 6; ++i)
        if (out.size() < 7)
            out.insert(lo + (hi - lo) * i / 6);
    return {out.begin(), out.end()};
}
// 按精确参数和单侧导数语义缓存；位置已存在但导数未求时按需补算。
// 在调用 SDK 前检查共享截止时间，预算计入首次求值和后续候选使用的全部采样工作。
const SurfaceSample &SurfaceSampleContext::At(double u, double v,
    sggk::DerivType sideU, sggk::DerivType sideV, bool needDerivatives)
{
    CheckDeadline("sampling");
    if (sideU == sggk::DerivType::NONE)
        sideU = u == source.MaxParamU() ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
    if (sideV == sggk::DerivType::NONE)
        sideV = v == source.MaxParamV() ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
    auto key = std::make_tuple(u, v, int(sideU), int(sideV));
    auto found = values.find(key);
    if (found != values.end() && (!needDerivatives || found->second.derivativesReady))
        return found->second;
    if (found == values.end() && values.size() >= budget)
        throw Failure(Status::VerificationFailed, "Shared source sampling budget exhausted");
    SurfaceSample x;
    x.u = u; x.v = v;
    if (found != values.end()) x=found->second;
    if (!needDerivatives)
    {
        x.point=Get(source.CalcPoint(u,v));
        if (!Finite(x.point)) throw Failure(Status::InvalidInput,"Nonfinite sampled surface point");
        return values.emplace(key,x).first->second;
    }
    sggk::Point3D p;
    auto d = source.CalcPartialD1(u, v, sideU, sideV, p);
    if (found==values.end()) x.point=Get(p);
    x.du = Get(d[0]); x.dv = Get(d[1]);x.derivativesReady=true;
    if (!Finite(x.point) || !Finite(x.du) || !Finite(x.dv))
        throw Failure(Status::InvalidInput, "Nonfinite sampled surface data");
    V cross = Cross(x.du, x.dv);
    x.area = Norm(cross);
    x.normalValid = x.area > 0 && x.area > 1e-12 * Norm(x.du) * Norm(x.dv);
    if (x.normalValid)
        x.normal = cross / x.area;
    if (found!=values.end()) {found->second=x;return found->second;}
    return values.emplace(key, x).first->second;
}
const SurfaceSample &SurfaceSampleContext::AtCurvePoint(
    double u, double v, const std::function<V()> &evaluate, bool &inserted)
{
    CheckDeadline("iso_sampling");
    auto su = u == source.MaxParamU() ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
    auto sv = v == source.MaxParamV() ? sggk::DerivType::LEFT : sggk::DerivType::RIGHT;
    auto key = std::make_tuple(u, v, int(su), int(sv));
    auto found = values.find(key);
    inserted = false;
    if (found != values.end()) return found->second;
    if (values.size() >= budget)
        throw Failure(Status::VerificationFailed, "Shared source sampling budget exhausted");
    SurfaceSample x;
    x.u = u; x.v = v; x.point = evaluate();
    if (!Finite(x.point)) throw Failure(Status::InvalidInput, "Nonfinite isocurve sample");
    inserted = true;
    return values.emplace(key, x).first->second;
}
std::vector<std::array<double, 2>> SurfaceSampleContext::ScoreParameters() const
{
    auto us = GrevilleParameters(source.KnotsU(), source.MultsU(), source.DegreeU(),
        source.ControlPoints().size(), source.MinParamU(), source.MaxParamU());
    auto vs = GrevilleParameters(source.KnotsV(), source.MultsV(), source.DegreeV(),
        source.ControlPoints().front().size(), source.MinParamV(), source.MaxParamV());
    if (us.size() > budget / vs.size())
        throw Failure(Status::VerificationFailed, "Greville grid exceeds sample budget");
    std::vector<std::array<double, 2>> out;
    for (size_t i = 1; i + 1 < us.size(); ++i)
        for (size_t j = 1; j + 1 < vs.size(); ++j)
            if ((i + j) % 2) out.push_back({us[i], vs[j]});
    return out;
}
const std::vector<SurfaceSample> &SurfaceSampleContext::Greville(bool needNormals)
{
    CheckDeadline("sampling");
    if (grevilleReady)
    {
        if (needNormals && !grevilleNormalsReady)
        {
            for (size_t i=0;i<greville.size();++i)
            {
                auto x=At(greville[i].u,greville[i].v);
                x.area*=grevilleWeights[i];greville[i]=x;
            }
            seeds.clear();
            for (size_t i:seedPartition) if (greville[i].normalValid) seeds.push_back(i);
            grevilleNormalsReady=true;
        }
        return greville;
    }
    auto us = GrevilleParameters(source.KnotsU(), source.MultsU(), source.DegreeU(),
        source.ControlPoints().size(), source.MinParamU(), source.MaxParamU());
    auto vs = GrevilleParameters(source.KnotsV(), source.MultsV(), source.DegreeV(),
        source.ControlPoints().front().size(), source.MinParamV(), source.MaxParamV());
    if (us.size() > budget / vs.size())
        throw Failure(Status::VerificationFailed, "Greville grid exceeds sample budget");
    greville.clear(); seeds.clear(); scores.clear();seedPartition.clear();grevilleWeights.clear();
    for (size_t i = 1; i + 1 < us.size(); ++i)
        for (size_t j = 1; j + 1 < vs.size(); ++j)
        {
            auto x = At(us[i], vs[j],sggk::DerivType::NONE,sggk::DerivType::NONE,needNormals);
            // 该参数单元所代表的近似面积权重。
            double weight=(us[i + 1] - us[i - 1]) * (vs[j + 1] - vs[j - 1]) * .25;
            grevilleWeights.push_back(weight);x.area*=weight;
            size_t index = greville.size();
            greville.push_back(x);
            if ((i + j) % 2) scores.push_back(index);
            else {seedPartition.push_back(index);if (x.normalValid) seeds.push_back(index);}
        }
    grevilleReady = true;
    grevilleNormalsReady=needNormals;
    return greville;
}
const std::vector<size_t> &SurfaceSampleContext::Seeds() { Greville(); return seeds; }
const std::vector<size_t> &SurfaceSampleContext::Scores() { Greville(false); return scores; }
const std::vector<size_t> &SurfaceSampleContext::CoarseScores(unsigned count)
{
    Scores();
    auto found = coarseScores.find(count);
    if (found != coarseScores.end()) return found->second;
    std::vector<size_t> selected;
    if (count >= scores.size()) selected = scores;
    else if (count && !scores.empty())
    {
        // 在归一化源 UV 中确定性地选取最远点以改善覆盖；只使用独立评分分区，不额外调用 SDK 求值。
        double du = source.MaxParamU() - source.MinParamU();
        double dv = source.MaxParamV() - source.MinParamV();
        std::vector<double> nearest(scores.size(), std::numeric_limits<double>::infinity());
        size_t next = 0;
        while (selected.size() < count)
        {
            size_t chosen = scores[next];
            selected.push_back(chosen);
            for (size_t j = 0; j < scores.size(); ++j)
            {
                double u = (greville[scores[j]].u - greville[chosen].u) / du;
                double v = (greville[scores[j]].v - greville[chosen].v) / dv;
                nearest[j] = std::min(nearest[j], u*u + v*v);
            }
            nearest[next] = -1;
            next = static_cast<size_t>(std::max_element(nearest.begin(), nearest.end()) - nearest.begin());
        }
    }
    return coarseScores.emplace(count, std::move(selected)).first->second;
}
const std::array<std::vector<SurfaceSample>, 4> &SurfaceSampleContext::Boundaries()
{
    if (boundariesReady) return boundaries;
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        boundaries[edge].clear();
        bool uFixed = edge < 2;
        auto ts = GrevilleParameters(uFixed ? source.KnotsV() : source.KnotsU(),
            uFixed ? source.MultsV() : source.MultsU(),
            uFixed ? source.DegreeV() : source.DegreeU(),
            uFixed ? source.ControlPoints().front().size() : source.ControlPoints().size(),
            uFixed ? source.MinParamV() : source.MinParamU(),
            uFixed ? source.MaxParamV() : source.MaxParamU());
        for (double t : ts)
            boundaries[edge].push_back(At(uFixed ? (edge == 0 ? source.MinParamU() : source.MaxParamU()) : t,
                uFixed ? t : (edge == 2 ? source.MinParamV() : source.MaxParamV()),
                sggk::DerivType::NONE,sggk::DerivType::NONE,false));
    }
    boundariesReady = true;
    return boundaries;
}
}
