#include "GeneralInternal.h"
#include "StageTimer.h"
namespace bspline2analytic::detail
{
static void CheckGeneralInput(const sggk::BSplineSurface &s, const GeneralConversionOptions &o)
{
    ValidateSurface(s);
    const auto &v = o.fitting.validation;
    const auto &r = o.fitting.refinement;
    auto nonnegative = [](double x) { return std::isfinite(x) && x >= 0; };
    if (!(v.distanceTolerance > 0 && std::isfinite(v.distanceTolerance)) ||
        v.maxSubdivisionDepth > 30 || !v.maxValidationPieces ||
        !s.DegreeU() || !s.DegreeV() ||
        !o.seedAttemptsPerType || o.seedAttemptsPerType > 4096 ||
        !o.candidatesPerType || o.candidatesPerType > 32 ||
        o.coarseScorePointCount > o.fitting.sampling.maxSamples ||
        o.boundarySamplesPerSpan<3 || o.boundarySamplesPerSpan>64 ||
        !o.fitting.sampling.maxSamples || r.maxIterations > 500 || r.maxRounds > 20 ||
        r.maxWarmupIterations > 500 || !r.stagnationWindow || r.stagnationWindow > 500 ||
        !r.maxConsecutiveFailures || r.maxConsecutiveFailures > 500 ||
        !(r.minTrustRadius > 0 && r.minTrustRadius <= .1) ||
        !nonnegative(r.worstImprovementAbsTol) || !nonnegative(r.worstImprovementRelTol) ||
        !nonnegative(r.costImprovementRelTol) || !nonnegative(r.smallModelChangeTolerance) ||
        !nonnegative(o.similarScoreAbsoluteRatio) || !nonnegative(o.similarScoreRelativeTolerance))
        throw Failure(Status::InvalidInput, "Invalid general conversion data/options");
    if (s.DegreeU() >= s.ControlPoints().size() ||
        s.DegreeV() >= s.ControlPoints().front().size())
        throw Failure(Status::InvalidInput, "Invalid surface degree/control dimensions");
    if (s.IsRational())
    {
        if (s.Weights().empty()) throw Failure(Status::InvalidInput, "Missing rational weights");
        for (const auto &row : s.Weights()) for (double w : row)
            if (!(w > 0)) throw Failure(Status::Unsupported, "General interval extraction requires positive weights");
    }
    for (size_t i = 1; i + 1 < s.MultsU().size(); ++i)
        if (s.MultsU()[i] > s.DegreeU()) throw Failure(Status::Unsupported, "Discontinuous U knot");
    for (size_t i = 1; i + 1 < s.MultsV().size(); ++i)
        if (s.MultsV()[i] > s.DegreeV()) throw Failure(Status::Unsupported, "Discontinuous V knot");
}
SampleCache GeneralSamples(SurfaceSampleContext &context)
{
    SampleCache samples;
    auto add = [&](const SurfaceSample &x) {
        Sample s; s.u = x.u; s.v = x.v; s.point = x.point;
        // No tangent row is generated for projected-image trims.
        samples.points.push_back(s);
    };
    for (const auto &x : context.Greville(false)) add(x);
    for (const auto &edge : context.Boundaries())
    {
        std::vector<V> positions;
        for (const auto &x : edge)
        {
            bool duplicate = false;
            for (V p : positions) duplicate |= Norm(x.point - p) <= Roundoff(Norm(p));
            if (!duplicate) { add(x); positions.push_back(x.point); }
        }
    }
    samples.RebuildIndices();
    return samples;
}
GeneralConversionResult ConvertGeneralWithSamples(const sggk::BSplineSurface &s,
    const GeneralConversionOptions &o, SurfaceSampleContext &context)
{
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    GeneralConversionResult result;
    try
    {
        CheckGeneralInput(s, o);
        GeneralCandidatePool pool;
        {
            StageTimer timer(result.timing.classification);
            pool = InitializeGeneralCandidates(context, o, result.classification);
        }
        result.message = "No general analytic candidate could be initialized";
        std::optional<SampleCache> base;
        unsigned index = 0;
        std::vector<GeneralSeed> seeds;
        size_t batchIndex = 0;
        for (;;)
        {
            if (batchIndex == seeds.size())
            {
                StageTimer timer(result.timing.classification);
                seeds = NextGeneralCandidateBatch(pool, context, o, result.classification);
                if (seeds.empty()) break;
                batchIndex = 0;
                for (const auto &seed : seeds) result.candidates.push_back(seed.score);
                result.statistics.candidateCount = static_cast<unsigned>(result.candidates.size());
            }
            const auto &seed = seeds[batchIndex++];
            const unsigned iterationsBefore=result.statistics.optimizationIterations;
            ++index;
            AnalyticPatch candidate;
            candidate.parameters = seed.parameters;
            // Valid scratch chart for reusing the legacy mother parameter update;
            // it does not define the projected-image patch's trims.
            candidate.domain = {0, 1, .25, 1};
            GeneralTrim trim;
            auto validate = [&] {
                StageTimer timer(result.timing.validation);
                result.validation = ValidateGeneral(s, candidate.parameters, o, context, &trim);
                result.statistics.validationPieces += result.validation.piecesVisited;
                result.chartAttempts += result.validation.chartAttempts;
                result.chartIterations += result.validation.chartIterations;
            };
            try
            {
                validate();
                if (!result.validation.Passed() && o.fitting.refinement.maxRounds &&
                    o.fitting.refinement.maxIterations)
                {
                    if (!base)
                    {
                        StageTimer timer(result.timing.optimization);
                        base = GeneralSamples(context);
                        ++result.statistics.baseSamplingBuilds;
                    }
                    auto samples = *base;
                    for (unsigned round = 0; round < o.fitting.refinement.maxRounds; ++round)
                    {
                        bool added = false;
                        if (result.validation.hasWitness)
                        {
                            auto uv = result.validation.witness;
                            SampleKey key{uv[0], uv[1], -1, int(sggk::DerivType::NONE), int(sggk::DerivType::NONE)};
                            if (!samples.keys.count(key) && samples.size() < o.fitting.sampling.maxSamples)
                            {
                                Sample witness; witness.u = uv[0]; witness.v = uv[1];
                                witness.point = context.At(uv[0], uv[1],sggk::DerivType::NONE,
                                    sggk::DerivType::NONE,false).point;
                                samples.points.push_back(witness); samples.RebuildIndices();
                                ++result.statistics.witnessSamplesAdded; added = true;
                            }
                        }
                        if (round && !added) break;
                        {
                            StageTimer timer(result.timing.optimization);
                            EvaluationWorkspace workspace;
                            auto problem = MotherRefinementProblem(samples, o.fitting);
                            auto refined=RefineCandidateLevel(candidate, samples, o.fitting, result.statistics,
                                ParameterLevel::MotherSurface, o.fitting.refinement.maxIterations,
                                true, workspace, &problem);
                            auto &diagnostic = result.statistics.refinementRuns.back();
                            diagnostic.candidateIndex = index; diagnostic.roundIndex = round + 1;
                            diagnostic.entryReason = round ? RefinementEntryReason::WitnessSamples : RefinementEntryReason::Initial;
                            ++result.statistics.refinementRounds;
                            result.statistics.sampleCount = std::max(result.statistics.sampleCount,
                                static_cast<unsigned>(samples.size()));
                            // The same mother and same validation options cannot
                            // repair a chart failure by repeating the same check.
                            if(!refined.improved) break;
                        }
                        validate();
                        if (result.validation.Passed()) break;
                    }
                }
                if (result.validation.Passed())
                {
                    result.patch.emplace(trim.parameters,
                        std::array<double,4>{s.MinParamU(),s.MaxParamU(),s.MinParamV(),s.MaxParamV()},
                        trim.corners,trim.edges);
                    result.status = ConversionStatus::Success;
                    result.message = result.validation.message;
                    result.attempts.push_back({seed.parameters.type,PCurveFailure::None,
                        result.validation.failedBoundary,true,result.message});
                    result.attempts.back().seedGroup=seed.score.seedGroup;
                    result.attempts.back().optimizationIterations=result.statistics.optimizationIterations-iterationsBefore;
                    break;
                }
                result.status = ConversionStatus::VerificationFailed;
                result.message = result.validation.message;
            }
            catch (const std::runtime_error &e)
            {
                result.status = ConversionStatus::VerificationFailed;
                result.message = e.what();
            }
            result.attempts.push_back({seed.parameters.type,result.validation.pcurveFailure,
                result.validation.failedBoundary,false,result.message});
            result.attempts.back().seedGroup=seed.score.seedGroup;
            result.attempts.back().optimizationIterations=result.statistics.optimizationIterations-iterationsBefore;
            // Finish the current ranked batch before drawing from the retained
            // pool. No new seed groups are generated during backup attempts.
        }
    }
    catch (const Failure &e)
    {
        result.status = e.status == Status::InvalidInput ? ConversionStatus::InvalidInput :
            e.status == Status::Unsupported ? ConversionStatus::Unsupported : ConversionStatus::VerificationFailed;
        result.message = e.what();
    }
    catch (const std::exception &e)
    {
        result.status = ConversionStatus::InvalidInput; result.message = e.what();
    }
    if (result.status != ConversionStatus::Success) result.patch.reset();
    result.timing.totalMs = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    result.timing.otherMs = result.timing.totalMs - result.timing.classification.elapsedMs -
        result.timing.optimization.elapsedMs - result.timing.validation.elapsedMs;
    return result;
}
}
namespace bspline2analytic
{
GeneralConversionResult ConvertGeneral(const sggk::BSplineSurface &s, const GeneralConversionOptions &o)
{
    detail::SurfaceSampleContext context(s, o.fitting.sampling.maxSamples);
    return detail::ConvertGeneralWithSamples(s, o, context);
}
}
