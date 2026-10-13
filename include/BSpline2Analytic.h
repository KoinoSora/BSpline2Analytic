#pragma once

#include <GeomBase/Dir3D.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <array>
#include <limits>
#include <optional>
#include <string>

namespace bspline2analytic
{
enum class Status
{
    Success,
    InvalidInput,
    Degenerate,
    IllConditioned,
    InconsistentBoundary,
    Ambiguous,
    VerificationFailed,
    Unsupported
};
enum class BoundaryId
{
    Umin,
    Umax,
    Vmin,
    Vmax
};
enum class BoundaryKind
{
    Unknown,
    Line,
    Circle,
    Degenerate
};
enum class ClassificationFailureReason
{
    None,
    DegenerateBoundary,
    NonLineOrCircleBoundary,
    InconsistentGeometry,
    NoCandidate
};
enum class VerificationLevel
{
    None,
    Sampled,
    BernsteinBound
};

struct RecognitionOptions
{
    // Length units must match those of the input SDK model (normally mm).
    double distanceTolerance = 1e-5;
    double rankTolerance = 1e-12;
    double angularTolerance = 1e-8;
    unsigned maxSubdivisionDepth = 12;
    unsigned maxValidationPieces = 20000;
    unsigned samplesPerSpan = 4;
    bool requireBound = true;
};

struct LineData
{
    sggk::Point3D origin;
    sggk::Dir3D direction{1, 0, 0};
    double start = 0, end = 0;
};

struct CircleData
{
    sggk::Point3D center;
    sggk::Dir3D normal{0, 0, 1};
    sggk::Dir3D xDirection{1, 0, 0};
    double radius = 0, startAngle = 0, sweepAngle = 0;
};

struct BoundaryInfo
{
    BoundaryId id = BoundaryId::Umin;
    BoundaryKind kind = BoundaryKind::Unknown;
    sggk::Curve3DPtr curve;
    LineData line;
    CircleData circle;
    double errorBound = std::numeric_limits<double>::infinity();
    std::string message;
};

struct BoundaryAnalysisResult
{
    Status status = Status::InvalidInput;
    std::array<BoundaryInfo, 4> boundaries;
    std::string message;
    ClassificationFailureReason failureReason = ClassificationFailureReason::None;
};

// A canonical geometric frame; it need not reproduce the source parameter origin.
struct AnalyticParameters
{
    sggk::SurfaceType type = sggk::SurfaceType::Unknown;
    sggk::Point3D origin; // plane origin, axis origin, cone apex, sphere/torus center
    sggk::Dir3D axis{0, 0, 1};
    sggk::Dir3D xDirection{1, 0, 0};
    double radius = 0;      // sphere/cylinder radius or torus tube radius
    double majorRadius = 0; // signed for SGK Lemon surfaces
    double halfAngle = 0;
};

struct ParameterCorrespondence
{
    bool swappedUV = false;
    bool reversedU = false;
    bool reversedV = false;
    // These are orientation flags, not an assumed affine NURBS-to-angle mapping.
    // Corner order: (umin,vmin), (umax,vmin), (umin,vmax), (umax,vmax).
    std::array<sggk::Point2D, 4> analyticCorners;
};

struct ValidationReport
{
    VerificationLevel level = VerificationLevel::None;
    double maxSampledDistance = 0;
    double maxBoundaryBound = 0;
    double surfaceDistanceBound = std::numeric_limits<double>::infinity();
    unsigned piecesVisited = 0;
    // BernsteinBound includes a floating-point guard, not directed interval arithmetic.
    std::string message;
};

struct StageTiming
{
    double elapsedMs = 0;
    unsigned calls = 0;
};
struct ConversionTiming
{
    StageTiming classification, optimization, validation;
    double otherMs = 0, totalMs = 0;
};

struct RecognitionResult
{
    Status status = Status::InvalidInput;
    AnalyticParameters parameters;
    sggk::SurfacePtr surface;
    std::array<BoundaryInfo, 4> boundaries;
    std::optional<sggk::UVRange> domain;
    ParameterCorrespondence mapping;
    ValidationReport validation;
    ConversionTiming timing;
    std::string message;
    ClassificationFailureReason classificationFailure = ClassificationFailureReason::None;
};

BoundaryAnalysisResult AnalyzeBoundaries(const sggk::BSplineSurface &surface,
                                         const RecognitionOptions &options = {});
RecognitionResult Recognize(const sggk::BSplineSurface &surface,
                            const RecognitionOptions &options = {});
const char *StatusName(Status status);
const char *TypeName(sggk::SurfaceType type);
} // namespace bspline2analytic
