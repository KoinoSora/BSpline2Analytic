// Independent test-only audits, adapted from the current specialized and rational suites.
double Map(const PairedParameterMap &m, double t)
{
    if (t == m.source.back()) return m.target.back();
    auto it = std::upper_bound(m.source.begin(), m.source.end(), t);
    Require(it != m.source.begin() && it != m.source.end(), "Audit map coverage");
    size_t i = size_t(it - m.source.begin() - 1);
    return m.target[i] + (t - m.source[i]) / (m.source[i + 1] - m.source[i]) * (m.target[i + 1] - m.target[i]);
}
// Audit in translated coordinates to avoid cancellation in SDK derivatives
// and world-point subtraction at 1e6/1e9 translations. This does not alter
// algorithm inputs or its measured candidate/certificates.
sggk::BSplineSurface CenteredAuditSurface(const sggk::BSplineSurface &s, V anchor)
{
    sggk::BSplineSurface out(s);
    const auto &cp = s.ControlPoints();
    for (unsigned i=0;i<cp.size();++i)
        for (unsigned j=0;j<cp[i].size();++j)
            out.SetControlPoint(Point(Get(cp[i][j])-anchor),i,j);
    return out;
}
void AuditPosition(const sggk::BSplineSurface &s, const AnalyticParameters &p, double bound, Row &row)
{
    V anchor = Get(s.ControlPoints().front().front());
    auto centered = CenteredAuditSurface(s,anchor);
    auto local = p; local.origin = Point(Get(p.origin)-anchor);
    SurfaceSampleContext samples(centered);
    double maximum = 0; bool valid = true;
    for (unsigned i = 0; i <= 12; ++i)
        for (unsigned j = 0; j <= 12; ++j)
        {
            double u = s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 12;
            double v = s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * j / 12;
            double d = Distance(local, Point(samples.At(u, v, sggk::DerivType::NONE, sggk::DerivType::NONE, false).point));
            valid &= std::isfinite(d) && d < 1e-5 && d <= bound + 1e-10;
            maximum = std::max(maximum, d);
        }
    Set(row, "audit_position_max", maximum);
    Require(valid,"Independent centered mother-distance audit");
}
void AuditEdges(const sggk::BSplineSurface &s, const AnalyticPatch &p,
                const MotherDistanceValidationReport &r, Row &row)
{
    auto chart = paired::BuildChart(p, false);
    V anchor = Get(s.ControlPoints().front().front());
    auto centered = CenteredAuditSurface(s,anchor);
    auto local = p.parameters; local.origin = Point(Get(local.origin)-anchor);
    auto frame = FrameOf(local);
    SurfaceSampleContext samples(centered);
    double maxPosition = 0, maxU = 0, maxV = 0, maxJoint = 0;
    bool positionOK = true, tangentOK = true;
    for (unsigned edge = 0; edge < 4; ++edge)
    {
        for (unsigned i = 0; i <= 16; ++i)
        {
            double u = edge < 2 ? (edge ? s.MaxParamU() : s.MinParamU())
                               : s.MinParamU() + (s.MaxParamU() - s.MinParamU()) * i / 16;
            double v = edge >= 2 ? (edge == 3 ? s.MaxParamV() : s.MinParamV())
                                : s.MinParamV() + (s.MaxParamV() - s.MinParamV()) * i / 16;
            double q[2] = {Map(r.uMap, u), Map(r.vMap, v)}, angle[2]{};
            if (p.mapping.reversedU) q[0] = 1 - q[0];
            if (p.mapping.reversedV) q[1] = 1 - q[1];
            if (p.mapping.swappedUV) std::swap(q[0], q[1]);
            for (unsigned d = 0; d < 2; ++d)
            {
                unsigned n = unsigned(chart.angles[d].size() - 1);
                double x = std::clamp(q[d], 0., 1.) * n;
                unsigned k = std::min(n - 1, unsigned(x));
                double t = x - k, a = chart.angles[d][k], b = chart.angles[d][k + 1];
                angle[d] = chart.angular[d] ? a + 2 * std::atan2(t * std::sin((b - a) / 2),
                    1 - t + t * std::cos((b - a) / 2)) : a + t * (b - a);
            }
            const auto &sample = samples.At(u, v);
            V target = frame.global(CanonicalPoint(p.parameters, angle[0], angle[1]));
            const double cu = std::cos(angle[0]), su = std::sin(angle[0]);
            const double cv = std::cos(angle[1]), sv = std::sin(angle[1]);
            const auto &a = p.parameters;
            V du{}, dv{};
            switch (a.type)
            {
            case sggk::SurfaceType::Plane: du = {1, 0, 0}; dv = {0, 1, 0}; break;
            case sggk::SurfaceType::Cylinder:
                du = {-a.radius * su, a.radius * cu, 0}; dv = {0, 0, 1}; break;
            case sggk::SurfaceType::Cone: {
                double k = std::tan(a.halfAngle);
                du = {-angle[1] * k * su, angle[1] * k * cu, 0}; dv = {k * cu, k * su, 1}; break;
            }
            case sggk::SurfaceType::Sphere:
                du = {-a.radius * sv * su, a.radius * sv * cu, 0};
                dv = {a.radius * cv * cu, a.radius * cv * su, a.radius * sv}; break;
            case sggk::SurfaceType::Torus:
                du = {-(a.majorRadius + a.radius * cv) * su, (a.majorRadius + a.radius * cv) * cu, 0};
                dv = {-a.radius * sv * cu, -a.radius * sv * su, a.radius * cv}; break;
            default: throw std::runtime_error("Unsupported audit mother");
            }
            if (p.mapping.swappedUV) std::swap(du, dv);
            du = frame.vector(du) * (p.mapping.reversedU ? -1. : 1.);
            dv = frame.vector(dv) * (p.mapping.reversedV ? -1. : 1.);
            double distance = Norm(sample.point - target);
            auto angleOf = [](V x, V y) {
                Require(Norm(x) > 0 && Norm(y) > 0, "Degenerate audit tangent");
                return std::atan2(Norm(Cross(x, y)), Dot(x, y));
            };
            double au = angleOf(sample.du, du), av = angleOf(sample.dv, dv);
            double joint = std::asin(std::min(1., std::sqrt(std::pow(std::sin(au), 2) +
                                                          std::pow(std::sin(av), 2)))) * 180 / Pi;
            positionOK &= std::isfinite(distance) && distance < 1e-5 && distance <= r.boundaryPosition.upperBound + 1e-10;
            tangentOK &= std::isfinite(au) && std::isfinite(av) && au < Pi / 2 && av < Pi / 2 &&
                au * 180 / Pi <= r.boundaryTangentU.upperBound + 1e-7 &&
                av * 180 / Pi <= r.boundaryTangentV.upperBound + 1e-7 &&
                joint < .1 && joint <= r.boundaryJointTangent.upperBound + 1e-7;
            maxPosition = std::max(maxPosition, distance);
            maxU = std::max(maxU, au * 180 / Pi); maxV = std::max(maxV, av * 180 / Pi);
            maxJoint = std::max(maxJoint, joint);
        }
    }
    Set(row, "audit_edge_max", maxPosition); Set(row, "audit_tangent_u_max", maxU);
    Set(row, "audit_tangent_v_max", maxV); Set(row, "audit_joint_tangent_max", maxJoint);
    Require(positionOK, "Independent centered paired edge-position audit");
    Require(tangentOK, "Independent centered joint U/V tangent audit");
}

double MotherDistance(const AnalyticParameters &a,V point)
{
    V p=FrameOf(a).local(point);double rho=std::hypot(p.x,p.y);
    using T=sggk::SurfaceType;
    switch(a.type)
    {
    case T::Plane:return std::abs(p.z);
    case T::Cylinder:return std::abs(rho-a.radius);
    case T::Cone:return std::abs(rho*std::cos(a.halfAngle)-p.z*std::sin(a.halfAngle));
    case T::Sphere:return std::abs(Norm(p)-a.radius);
    default:return std::abs(std::hypot(rho-a.majorRadius,p.z)-a.radius);
    }
}
double AuditGeneral(const sggk::BSplineSurface &s,const GeneralRationalBoundaryReport &r,double motherTolerance=1e-8)
{
    Require(r.Passed(),"Cannot audit failed output");double maximum=0;const auto &patch=*r.patch;
    SurfaceSampleContext samples(s);
    for(unsigned e=0;e<4;++e)
    {
        const auto &curve=*patch.boundaries[e];
        for(unsigned endpoint=0;endpoint<2;++endpoint)
        {
            auto p=Get(curve.ControlPoints()[endpoint?curve.ControlPoints().size()-1:0]);
            auto v=Get(patch.vertices[rb::Ends[e][endpoint]]);
            Require(p.x==v.x&&p.y==v.y&&p.z==v.z,"Output wire endpoint mismatch");
        }
        for(unsigned i=0;i<=32;++i)
        {
            double f=double(i)/32;
            Require(MotherDistance(patch.mother,Get(curve.CalcPoint(f)))<motherTolerance,"Exported boundary off mother");
            double u=e<2?(e==0?s.MinParamU():s.MaxParamU()):s.MinParamU()+(s.MaxParamU()-s.MinParamU())*f;
            double v=e<2?s.MinParamV()+(s.MaxParamV()-s.MinParamV())*f:(e==2?s.MinParamV():s.MaxParamV());
            auto point=samples.At(u,v);double parameter=f;
            auto nearest=curve.CalcNearestPoint(Point(point.point),parameter);
            double distance=Norm(point.point-Get(nearest));maximum=std::max(maximum,distance);
            Require(distance<1.001e-5,"Independent SDK boundary distance audit failed");
            bool exempt=false;
            for(const auto &piece:patch.tangentExemptions)
                if(piece.edge==e&&(e<2?v:u)>=piece.sourceFirst&&(e<2?v:u)<=piece.sourceLast)exempt=true;
            if(!patch.tangentNotApplicable[e]&&!exempt&&i>0&&i<32)
            {
                V a=e<2?point.dv:point.du,b=Get(curve.CalcDeriv1(parameter));
                Require(Dot(a,b)>0&&Norm(Cross(a,b))/Dot(a,b)<std::tan(.1*Pi/180),"Independent tangent audit failed");
            }
        }
    }
    return maximum;
}
