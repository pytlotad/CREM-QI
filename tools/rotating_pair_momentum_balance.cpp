// Is the third-order residual of audit 142 a defect, or the momentum the
// FIELD and the RADIATION carry (audit section 564)?
//
// The synthetic circular pair of audits 122-142 (1 r*, tilt 0, exact circular
// kinematics in the history, moments along the axis) is rigidly rotating:
// every field, the hidden momentum and the interaction field momentum turn
// with angular velocity omega z.  Momentum conservation for the mutual
// (cross) part of the fields then says, with the engine force on a moment
// already carrying -d p_hid/dt (audit 123),
//
//   F1 + F2 = oint T_cross . n dA - omega z x (P_V + p_hid),
//
// P_V = eps0 int_{|x|<R} (E1 x B2 + E2 x B1) dV the interaction field
// momentum inside the sphere R about the centre, T_cross the cross terms of
// the Maxwell stress on that sphere, p_hid = sum m_i x E_j(r_i)/c^2.  The
// right side is computed from the model's own field functions -- no new
// physics -- and must not depend on R.  If it equals the engine's F1 + F2,
// the residual is physical (the field's and the radiation's momentum) and
// nothing is missing; if not, the difference is a defect of the force law.
//
// Quadrature: a smooth partition of unity -- a C-infinity bump of radius rho
// about each particle (integrated on a grid centred on that particle, which
// absorbs the 1/r^2 of the point charge field), and the rest on a grid
// centred on the centre of mass, cut at R.
//
// Usage: rotating_pair_momentum_balance [beta] [history nodes] [span r/c]
//        [R in separations] [resolution multiplier] [separation in r*]
// Build from the repository root (OpenMP for the volume integral):
// g++ -std=c++20 -O2 -fopenmp -I . $(root-config --cflags)
//     tools/rotating_pair_momentum_balance.cpp -o /tmp/rotating
//     $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
struct Rule { std::vector<double> x,w; };
// Gauss-Legendre on [a,b].
Rule gauss(int n,double a,double b) {
    Rule r; r.x.resize(n); r.w.resize(n);
    for(int i=0;i<n;++i) {
        double x=std::cos(pi*(i+0.75)/(n+0.5)), d=0.0;
        for(int it=0;it<100;++it) {
            double p0=1.0,p1=0.0;
            for(int k=0;k<n;++k){ const double p2=p1; p1=p0;
                p0=((2.0*k+1.0)*x*p1-k*p2)/(k+1.0); }
            d=n*(x*p0-p1)/(x*x-1.0);
            const double s=p0/d; x-=s; if(std::abs(s)<1e-15) break;
        }
        double p0=1.0,p1=0.0;
        for(int k=0;k<n;++k){ const double p2=p1; p1=p0;
            p0=((2.0*k+1.0)*x*p1-k*p2)/(k+1.0); }
        d=n*(x*p0-p1)/(x*x-1.0);
        r.x[i]=0.5*(b-a)*x+0.5*(b+a);
        r.w[i]=0.5*(b-a)*2.0/((1.0-x*x)*d*d);
    }
    return r;
}
Rule composite(const std::vector<double>& edges,int perSegment) {
    Rule all;
    for(size_t s=0;s+1<edges.size();++s) {
        const Rule g=gauss(perSegment,edges[s],edges[s+1]);
        all.x.insert(all.x.end(),g.x.begin(),g.x.end());
        all.w.insert(all.w.end(),g.w.begin(),g.w.end());
    }
    return all;
}
// C-infinity step: 1 for u <= 1/2, 0 for u >= 1.
double bump(double u) {
    if(u<=0.5) return 1.0;
    if(u>=1.0) return 0.0;
    const double t=(u-0.5)/0.5;
    const auto f=[](double y){ return y>0.0?std::exp(-1.0/y):0.0; };
    return f(1.0-t)/(f(1.0-t)+f(t));
}
}

static Vec3 gDistributedDifference; static bool gDistributedValid=false;
int main(int argc,char** argv) {
    // PAIR (audit 571): e.g. proton,electron; arms stay equal (synthetic kinematics)
    if(const char* pair=std::getenv("PAIR")) applyPairFromOption(pair);
    const double beta=argc>1?atof(argv[1]):0.1;
    const int nodes=argc>2?atoi(argv[2]):65;
    const double span=argc>3?atof(argv[3]):16.0;
    const double outerFactor=argc>4?atof(argv[4]):5.0;
    const double resolution=argc>5?atof(argv[5]):1.0;
    // argument 6: separation in r* (audit 564, uzupelnienie 1); R and the
    // bumps scale with it, the history span is in separation/c already
    const double separation=(argc>6?atof(argv[6]):1.0)*comptonBarrierRadius;
    const Vec3 axis{0.0,0.0,1.0};
    const double armRadius=0.5*separation;
    const double angular=beta*c/armRadius;
    // MOMENT_SCALE multiplies both moments (audit 571)
    const double momentScale=std::getenv("MOMENT_SCALE")?std::atof(std::getenv("MOMENT_SCALE")):1.0;
    const auto build=[&](double time) {
        State s{};
        s.time=time;
        const double phase=angular*time;
        const Vec3 offset{armRadius*std::cos(phase),armRadius*std::sin(phase),0.0};
        const Vec3 velocity{-armRadius*angular*std::sin(phase),
                             armRadius*angular*std::cos(phase),0.0};
        const Vec3 acceleration=offset*(-angular*angular);
        s.firstPosition=offset; s.secondPosition=offset*-1.0;
        s.firstVelocity=velocity; s.secondVelocity=velocity*-1.0;
        s.firstAcceleration=acceleration; s.secondAcceleration=acceleration*-1.0;
        s.firstProperDipole=axis*(firstMagneticMoment*momentScale);
        s.secondProperDipole=axis*(secondMagneticMoment*momentScale);
        synchronizeCovariantDipoles(s);
        return s;
    };
    StateHistory history;
    const double lightTime=separation/c;
    for(int i=0;i<nodes;++i)
        history.push_back(build(-span*lightTime*(1.0-static_cast<double>(i)/nodes)));
    const State present=build(0.0);
    const MutualForces forces=retardedExternalForces(present,history);
    const Vec3 net=forces.first+forces.second;
    const Vec3 omega=axis*angular;
    // the engine's own hidden-momentum rate force, and p_hid with the field
    // the engine uses for it (matched softening, audit 126)
    const Vec3 rateForce=hiddenMomentumRateForce(present,history,true)
        +hiddenMomentumRateForce(present,history,false);
    const auto electricAt=[&](bool targetIsFirst) {
        return fieldFromOtherParticleAt(targetIsFirst?present.firstPosition
            :present.secondPosition,present.time,present,history,
            targetIsFirst,1.0e-5,matchedMomentSoftening()).electric;
    };
    const Vec3 hidden=(cross(present.firstDipole,electricAt(true))
        +cross(present.secondDipole,electricAt(false)))*(1.0/(c*c));
    // fields of each source at a point (charge LW + moment), point charge
    const auto sources=[&](const Vec3& p,ElectromagneticField& one,
                           ElectromagneticField& two) {
        one=fieldFromOtherParticleAt(p,present.time,present,history,false);
        two=fieldFromOtherParticleAt(p,present.time,present,history,true);
    };
    const auto density=[&](const Vec3& p) {
        ElectromagneticField one,two;
        sources(p,one,two);
        return (cross(one.electric,two.magnetic)
               +cross(two.electric,one.magnetic))*epsilon0;
    };
    const double outer=outerFactor*separation;
    const double rho=0.4*separation;
    // AUDIT 572: the momentum source of fields that break Faraday or div B,
    //   S = eps0 (E1 x K2 + E2 x K1) - (1/mu0)(b1 B2 + b2 B1),
    //   K_i = curl E_i + dB_i/dt, b_i = div B_i, from differences on the ONE
    // history (all times <= 0, evaluated at t = -2 dt).  densityKind selects
    // what the volume quadrature integrates (0 momentum, 1 S, 2 b-part only).
    int densityKind=0;
    const double sourceStep=1.0e-4*separation, sourceTimeStep=sourceStep/c;
    const auto sourceDensity=[&](const Vec3& p)->Vec3{
        const double te=-2.0*sourceTimeStep;
        const auto fieldOf=[&](bool fromFirst,const Vec3& x,double t){
            return fieldFromOtherParticleAt(x,t,present,history,!fromFirst); };
        Vec3 E[2],B[2],K[2],A[2]; double b[2],a[2];
        const Vec3 e[3]={{sourceStep,0,0},{0,sourceStep,0},{0,0,sourceStep}};
        for(int i=0;i<2;++i){
            const bool fromFirst=i==0;
            const ElectromagneticField f0=fieldOf(fromFirst,p,te);
            Vec3 dE[3],dB[3];
            for(int k=0;k<3;++k){ const ElectromagneticField a=fieldOf(fromFirst,p+e[k],te), m=fieldOf(fromFirst,p-e[k],te);
                dE[k]=(a.electric-m.electric)*(1.0/(2*sourceStep)); dB[k]=(a.magnetic-m.magnetic)*(1.0/(2*sourceStep)); }
            const ElectromagneticField fp=fieldOf(fromFirst,p,te+sourceTimeStep), fm=fieldOf(fromFirst,p,te-sourceTimeStep);
            const Vec3 curlE{dE[1].z-dE[2].y,dE[2].x-dE[0].z,dE[0].y-dE[1].x};
            K[i]=curlE+(fp.magnetic-fm.magnetic)*(1.0/(2*sourceTimeStep));
            b[i]=dB[0].x+dB[1].y+dB[2].z;
            const Vec3 curlB{dB[1].z-dB[2].y,dB[2].x-dB[0].z,dB[0].y-dB[1].x};
            A[i]=curlB-(fp.electric-fm.electric)*(1.0/(2*sourceTimeStep*c*c));
            a[i]=dE[0].x+dE[1].y+dE[2].z;
            E[i]=f0.electric; B[i]=f0.magnetic;
        }
        const Vec3 kPart=(cross(E[0],K[1])+cross(E[1],K[0]))*epsilon0;
        const Vec3 bPart=(B[1]*b[0]+B[0]*b[1])*(-1.0/mu0);
        // Ampere-Maxwell and div E residuals with J = 0, rho = 0 assumed: they
        // carry the smeared moment's physical sources as well (audit 572, report)
        const Vec3 aPart=(cross(A[0],B[1])+cross(A[1],B[0]))*(-1.0/mu0);
        const Vec3 ePart=(E[1]*a[0]+E[0]*a[1])*(-epsilon0);
        if(densityKind==3) return kPart+bPart+aPart+ePart;
        // audit 573: Lorentz force of the partner's field on moment k's own
        // smeared sources (rho = eps0 div E_m, J = curl B_m/mu0 - eps0 dE_m/dt)
        if(densityKind==6||densityKind==7) {
            const bool momentFirst=densityKind==6;
            const auto momentField=[&](const Vec3& x,double t){
                return retardedMagneticDipoleField(x,t,history,present,momentFirst); };
            const ElectromagneticField m0=momentField(p,te);
            Vec3 dE[3],dB[3];
            for(int k=0;k<3;++k){ const ElectromagneticField u=momentField(p+e[k],te), d=momentField(p-e[k],te);
                dE[k]=(u.electric-d.electric)*(1.0/(2*sourceStep)); dB[k]=(u.magnetic-d.magnetic)*(1.0/(2*sourceStep)); }
            const ElectromagneticField mp=momentField(p,te+sourceTimeStep), mm=momentField(p,te-sourceTimeStep);
            const double rho=epsilon0*(dE[0].x+dE[1].y+dE[2].z);
            const Vec3 curlB{dB[1].z-dB[2].y,dB[2].x-dB[0].z,dB[0].y-dB[1].x};
            const Vec3 J=curlB*(1.0/mu0)-(mp.electric-mm.electric)*(epsilon0/(2*sourceTimeStep));
            const ElectromagneticField other=fieldFromOtherParticleAt(p,te,present,history,momentFirst);
            (void)m0;
            return other.electric*rho+cross(J,other.magnetic);
        }
        if(densityKind==4) return aPart;
        if(densityKind==5) return ePart;
        return densityKind==2?bPart:kPart+bPart; };
    const auto spherical=[&](const Vec3& centre,const Rule& radial,int nPolar,
                             int nAzimuth,auto weightOf) {
        const Rule polar=gauss(nPolar,-1.0,1.0);
        std::vector<Vec3> partial(radial.x.size());
        #pragma omp parallel for schedule(dynamic)
        for(int ir=0;ir<static_cast<int>(radial.x.size());++ir) {
            const double r=radial.x[ir];
            Vec3 sum;
            for(int ip=0;ip<nPolar;++ip) {
                const double ct=polar.x[ip], st=std::sqrt(std::max(0.0,1.0-ct*ct));
                for(int ia=0;ia<nAzimuth;++ia) {
                    const double ph=2.0*pi*(ia+0.5)/nAzimuth;
                    const Vec3 p=centre+Vec3{r*st*std::cos(ph),r*st*std::sin(ph),r*ct};
                    const double w=weightOf(p);
                    if(w==0.0) continue;
                    const Vec3 d=densityKind==0?density(p):sourceDensity(p);
                    if(isFinite(d)) sum+=d*(w*polar.w[ip]);
                }
            }
            partial[ir]=sum*(radial.w[ir]*r*r*2.0*pi/nAzimuth);
        }
        Vec3 total; for(const Vec3& v:partial) total+=v;
        return total;
    };
    const auto bumpOf=[&](const Vec3& p,const Vec3& centre){
        return bump((p-centre).norm()/rho); };
    const int nr=static_cast<int>(24*resolution), np=static_cast<int>(32*resolution),
              na=static_cast<int>(64*resolution);
    const Rule nearRadial=composite({0.0,0.25*rho,0.5*rho,rho},nr/2);
    const Vec3 nearOne=spherical(present.firstPosition,nearRadial,np,na,
        [&](const Vec3& p){ return bumpOf(p,present.firstPosition); });
    const Vec3 nearTwo=spherical(present.secondPosition,nearRadial,np,na,
        [&](const Vec3& p){ return bumpOf(p,present.secondPosition); });
    std::vector<double> edges{0.0,0.1*separation,0.3*separation,0.45*separation,
        0.55*separation,0.7*separation,0.9*separation,1.2*separation,
        2.0*separation};
    for(double e=3.0*separation;e<outer;e*=1.5) edges.push_back(e);
    edges.push_back(outer);
    const Rule farRadial=composite(edges,nr/2);
    const Vec3 farPart=spherical(Vec3{},farRadial,2*np,2*na,
        [&](const Vec3& p){ return 1.0-bumpOf(p,present.firstPosition)
                                     -bumpOf(p,present.secondPosition); });
    const Vec3 fieldMomentum=nearOne+nearTwo+farPart;
    // momentum flux of the cross terms through the sphere: oint T . n dA
    Vec3 flux;
    {
        const Rule polar=gauss(2*np,-1.0,1.0);
        const int nAz=2*na;
        for(int ip=0;ip<2*np;++ip) {
            const double ct=polar.x[ip], st=std::sqrt(std::max(0.0,1.0-ct*ct));
            for(int ia=0;ia<nAz;++ia) {
                const double ph=2.0*pi*(ia+0.5)/nAz;
                const Vec3 n{st*std::cos(ph),st*std::sin(ph),ct};
                ElectromagneticField one,two;
                sources(n*outer,one,two);
                const auto stress=[&](const Vec3& a,const Vec3& b,double k){
                    return (a*dot(b,n)+b*dot(a,n)-n*dot(a,b))*k; };
                const Vec3 t=stress(one.electric,two.electric,epsilon0)
                    +stress(one.magnetic,two.magnetic,1.0/mu0);
                flux+=t*(polar.w[ip]*outer*outer*2.0*pi/nAz);
            }
        }
    }
    const Vec3 expected=flux-cross(omega,fieldMomentum+hidden);
    std::printf("beta %.4f  nodes %d  span %.0f r/c  R %.2f r*  resolution %.2f"
                "  separation %.2f r*\n",beta,nodes,span,
                outer/comptonBarrierRadius,resolution,
                separation/comptonBarrierRadius);
    std::printf("  engine F1+F2              %+.6e %+.6e %+.6e  |.| %.6e N"
                "  (|F1| %.6e N, ratio %.4e)\n",net.x,net.y,net.z,net.norm(),
                forces.first.norm(),net.norm()/forces.first.norm());
    std::printf("  engine hidden rate force  %+.6e %+.6e %+.6e  |.| %.6e N\n",
                rateForce.x,rateForce.y,rateForce.z,rateForce.norm());
    const Vec3 rotHidden=cross(omega,hidden);
    std::printf("  -omega x p_hid            %+.6e %+.6e %+.6e  (rate-force check |diff|/|.| %.3e)\n",
                -rotHidden.x,-rotHidden.y,-rotHidden.z,
                (rateForce+rotHidden).norm()/std::max(rotHidden.norm(),1e-300));
    std::printf("  p_hid                     %+.6e %+.6e %+.6e\n",hidden.x,hidden.y,hidden.z);
    std::printf("  P_V (near1 near2 far)     %+.6e %+.6e %+.6e  (|near| %.3e %.3e |far| %.3e)\n",
                fieldMomentum.x,fieldMomentum.y,fieldMomentum.z,
                nearOne.norm(),nearTwo.norm(),farPart.norm());
    std::printf("  |P_V + p_hid|/|p_hid|     %.6e\n",
                (fieldMomentum+hidden).norm()/std::max(hidden.norm(),1e-300));
    std::printf("  -omega x (P_V + p_hid)    %+.6e %+.6e %+.6e\n",
                -cross(omega,fieldMomentum+hidden).x,-cross(omega,fieldMomentum+hidden).y,
                -cross(omega,fieldMomentum+hidden).z);
    std::printf("  flux oint T.n dA          %+.6e %+.6e %+.6e\n",flux.x,flux.y,flux.z);
    std::printf("  expected                  %+.6e %+.6e %+.6e  |.| %.6e N\n",
                expected.x,expected.y,expected.z,expected.norm());
    std::printf("  |F1+F2 - expected|/|F1+F2| %.4e\n",
                (net-expected).norm()/std::max(net.norm(),1e-300));
    // AUDIT 570: the expectation with the engine's covariant spin momenta --
    // hidden m_0 x E'/c^2 and Thomas q (E' x S_0)/(m c^2), both rest-frame,
    // boosted as four-vectors.  Rigid rotation: dp/dt = omega x p, dp^0/dt = 0,
    // and the engine's force is -dp/dtau projected orthogonal to u, so
    //   F_spin = -omega x p - v gamma^2 ((omega x p).v)/c^2.
    {
        Vec3 spinForce, spinMomentum;
        for(int which=0;which<2;++which) {
            const bool first=which==0;
            const Vec3 x0=first?present.firstPosition:present.secondPosition;
            const Vec3 v=first?present.firstVelocity:present.secondVelocity;
            const Vec3 m0=first?present.firstProperDipole:present.secondProperDipole;
            const double gyro=first?firstGyromagneticRatioOf():secondGyromagneticRatioOf();
            const double q=first?firstCharge:secondCharge, mass=first?firstMass:secondMass;
            const double gu=1.0/std::sqrt(1.0-v.squaredNorm()/(c*c));
            const auto restField=[&](double floor){
                const ElectromagneticField f=fieldFromOtherParticleAt(x0,present.time,present,history,first,1.0e-5,floor);
                return (f.electric+cross(v,f.magnetic))*gu-v*(gu*gu/(gu+1.0)*dot(v,f.electric)/(c*c)); };
            const Vec3 rest=cross(m0,restField(matchedMomentSoftening()))*(1.0/(c*c))
                +cross(restField(0.0),m0*(1.0/gyro))*(q/(mass*c*c));
            const double vn=v.norm(); const Vec3 vh=vn>0?v*(1.0/vn):Vec3{};
            const Vec3 ps=rest+vh*((gu-1.0)*dot(rest,vh));
            spinMomentum+=ps;
            spinForce+=cross(omega,ps)*-1.0-v*(gu*gu*dot(cross(omega,ps),v)/(c*c));
        }
        const Vec3 expectedCov=flux-cross(omega,fieldMomentum)+spinForce;
        std::printf("  [570] spin momenta (cov)   %+.6e %+.6e %+.6e\n",spinMomentum.x,spinMomentum.y,spinMomentum.z);
        std::printf("  [570] expected (covariant) %+.6e %+.6e %+.6e  |.| %.6e N\n",expectedCov.x,expectedCov.y,expectedCov.z,expectedCov.norm());
        std::printf("  [570] |F1+F2 - expected_cov|/|F1+F2| %.4e\n",(net-expectedCov).norm()/std::max(net.norm(),1e-300));
        // AUDIT 571: split the residual into the spin-momentum part (engine
        // stencil against the rigid-rotation formula) and the field part
        const Vec3 engineSpin=covariantSpinMomentumForce(present,history,true)
            +covariantSpinMomentumForce(present,history,false);
        const Vec3 total=net-expectedCov;
        const Vec3 rSpin=engineSpin-spinForce;
        const Vec3 rField=(net-engineSpin)-(flux-cross(omega,fieldMomentum));
        const double tn=std::max(total.norm(),1e-300);
        std::printf("  [571] residual %+.6e %+.6e %+.6e |.| %.6e\n",total.x,total.y,total.z,total.norm());
        std::printf("  [571] R_spin   %+.6e %+.6e %+.6e  share %.4f\n",rSpin.x,rSpin.y,rSpin.z,dot(rSpin,total)/(tn*tn));
        std::printf("  [571] R_field  %+.6e %+.6e %+.6e  share %.4f\n",rField.x,rField.y,rField.z,dot(rField,total)/(tn*tn));
        std::printf("  [571] identity |R_spin + R_field - residual|/|residual| %.3e\n",(rSpin+rField-total).norm()/tn);
        std::printf("  [571] |engine spin| %.6e  |analytic spin| %.6e  |F1+F2| %.6e\n",engineSpin.norm(),spinForce.norm(),net.norm());
        if(std::getenv("DISTRIBUTED_FORCE")) {
            // whole sphere R with the partition of unity of P_V (uzupelnienie 1)
            const auto whole=[&](int kind){
                densityKind=kind;
                const Vec3 a=spherical(present.firstPosition,nearRadial,np,na,
                    [&](const Vec3& q){ return bumpOf(q,present.firstPosition); });
                const Vec3 b2=spherical(present.secondPosition,nearRadial,np,na,
                    [&](const Vec3& q){ return bumpOf(q,present.secondPosition); });
                const Vec3 f=spherical(Vec3{},farRadial,2*np,2*na,
                    [&](const Vec3& q){ return 1.0-bumpOf(q,present.firstPosition)-bumpOf(q,present.secondPosition); });
                densityKind=0; return a+b2+f; };
            const Vec3 d1=whole(6), d2=whole(7);
            const Vec3 g1=covariantDipoleGradientForce(present,history,true), g2=covariantDipoleGradientForce(present,history,false);
            const Vec3 diff=(g1+g2)-(d1+d2);
            const double rn=std::max(rField.norm(),1e-300);
            std::printf("  [573] F_grad sum  %+.6e %+.6e %+.6e  |.| %.6e\n",(g1+g2).x,(g1+g2).y,(g1+g2).z,(g1+g2).norm());
            std::printf("  [573] F_distr sum %+.6e %+.6e %+.6e  |.| %.6e\n",(d1+d2).x,(d1+d2).y,(d1+d2).z,(d1+d2).norm());
            std::printf("  [573] per particle |F_grad - F_distr|/|F_grad|: %.4e %.4e\n",(g1-d1).norm()/std::max(g1.norm(),1e-300),(g2-d2).norm()/std::max(g2.norm(),1e-300));
            std::printf("  [573] sum(F_grad - F_distr) %+.6e %+.6e %+.6e  |.|/|sum F_grad| %.4e  projection on R_field %.4f  |R - D|/|R| %.4f\n",
                diff.x,diff.y,diff.z,diff.norm()/std::max((g1+g2).norm(),1e-300),dot(diff,rField)/(rn*rn),(rField-diff).norm()/rn);
            gDistributedDifference=diff; gDistributedValid=true;
        }
        if(std::getenv("FARADAY_SOURCE")) {
            const auto integrate=[&](int kind){
                densityKind=kind;
                const Vec3 a=spherical(present.firstPosition,nearRadial,np,na,
                    [&](const Vec3& q){ return bumpOf(q,present.firstPosition); });
                const Vec3 b2=spherical(present.secondPosition,nearRadial,np,na,
                    [&](const Vec3& q){ return bumpOf(q,present.secondPosition); });
                const Vec3 f=spherical(Vec3{},farRadial,2*np,2*na,
                    [&](const Vec3& q){ return 1.0-bumpOf(q,present.firstPosition)-bumpOf(q,present.secondPosition); });
                densityKind=0; return a+b2+f; };
            const Vec3 S=integrate(1);
            const double rn=std::max(rField.norm(),1e-300);
            std::printf("  [572] int S      %+.6e %+.6e %+.6e  projection on R_field %.4f  |R_field - S|/|R_field| %.4f\n",
                S.x,S.y,S.z,dot(S,rField)/(rn*rn),(rField-S).norm()/rn);
            if(gDistributedValid) {
                const Vec3 P=gDistributedDifference+S;
                std::printf("  [573] prediction D + int S %+.6e %+.6e %+.6e  projection on R_field %.4f  |R_field - P|/|R_field| %.4f\n",
                    P.x,P.y,P.z,dot(P,rField)/(rn*rn),(rField-P).norm()/rn);
            }
            if(std::getenv("FULL_SOURCE")) {
                const Vec3 Sa=integrate(4), Se=integrate(5), St=S+Sa+Se;
                std::printf("  [572] int S_A    %+.6e %+.6e %+.6e  (Ampere-Maxwell part, J = 0)\n",Sa.x,Sa.y,Sa.z);
                std::printf("  [572] int S_divE %+.6e %+.6e %+.6e  (div E part, rho = 0)\n",Se.x,Se.y,Se.z);
                std::printf("  [572] int S_all  %+.6e %+.6e %+.6e  projection on R_field %.4f  |R_field - S_all|/|R_field| %.4f\n",
                    St.x,St.y,St.z,dot(St,rField)/(rn*rn),(rField-St).norm()/rn);
            }
            if(std::getenv("FARADAY_SOURCE_SPLIT")) {
                const Vec3 Sb=integrate(2);
                std::printf("  [572] int S_b    %+.6e %+.6e %+.6e  (div B part)\n",Sb.x,Sb.y,Sb.z);
            }
        }
        // candidates for R_field
        Vec3 c1,c3,c4;
        const double hs=1.0e-4*separation;
        for(int which=0;which<2;++which){
            const bool first=which==0;
            const Vec3 x0=first?present.firstPosition:present.secondPosition;
            const Vec3 v=first?present.firstVelocity:present.secondVelocity;
            const Vec3 p=first?present.firstElectricDipole:present.secondElectricDipole;
            const Vec3 m=first?present.firstDipole:present.secondDipole;
            const auto F=[&](const Vec3& x){ return fieldFromOtherParticleAt(x,present.time,present,history,first,1.0e-5,matchedMomentSoftening()); };
            Vec3 dE[3],dB[3]; const Vec3 e[3]={{hs,0,0},{0,hs,0},{0,0,hs}};
            for(int k=0;k<3;++k){ const ElectromagneticField a=F(x0+e[k]), b=F(x0-e[k]);
                dE[k]=(a.electric-b.electric)*(1.0/(2*hs)); dB[k]=(a.magnetic-b.magnetic)*(1.0/(2*hs)); }
            const Vec3 curlE{dE[1].z-dE[2].y,dE[2].x-dE[0].z,dE[0].y-dE[1].x};
            const Vec3 curlB{dB[1].z-dB[2].y,dB[2].x-dB[0].z,dB[0].y-dB[1].x};
            const Vec3 pGradE=dE[0]*p.x+dE[1]*p.y+dE[2]*p.z, pGradB=dB[0]*p.x+dB[1]*p.y+dB[2]*p.z;
            Vec3 gradPE{dot(p,dE[0]),dot(p,dE[1]),dot(p,dE[2])};
            const ElectromagneticField here=F(x0);
            c1+=pGradE+cross(cross(omega,p),here.magnetic)+cross(v,pGradB)-gradPE;
            c3+=cross(m,curlB)*-1.0; c4+=cross(p,curlE)*-1.0;
        }
        const auto share=[&](const char* n,const Vec3& k){
            std::printf("  [571] cand %-6s %+.6e %+.6e %+.6e  share %.4f  |R-K|/|R| %.4f\n",n,k.x,k.y,k.z,
                dot(k,rField)/std::max(rField.squaredNorm(),1e-300),(rField-k).norm()/std::max(rField.norm(),1e-300)); };
        std::printf("  [571] R_field abs %.6e N  (|F1| %.6e)\n",rField.norm(),forces.first.norm());
        share("K1",c1); share("K3",c3); share("K4",c4); share("K3+K4",c3+c4); share("K1+K3",c1+c3);
    }
    // CANDIDATES for the missing third-order term (audit 566).  Delta is
    // what the balance asks the forces to add.
    const Vec3 delta=expected-net;
    const double h=1.0e-4*separation;
    const auto partner=[&](bool targetIsFirst,const Vec3& x){
        return fieldFromOtherParticleAt(x,present.time,present,history,
            targetIsFirst,1.0e-5,matchedMomentSoftening()); };
    Vec3 k1, k2, k2a, k2b;
    for(int which=0;which<2;++which) {
        const bool first=which==0;
        const Vec3 x0=first?present.firstPosition:present.secondPosition;
        const Vec3 v=first?present.firstVelocity:present.secondVelocity;
        const Vec3 p=first?present.firstElectricDipole:present.secondElectricDipole;
        const Vec3 mLab=first?present.firstDipole:present.secondDipole;
        const Vec3 m0=first?present.firstProperDipole:present.secondProperDipole;
        Vec3 dE[3], dB[3];   // d/dx_k of E and B
        const Vec3 e[3]={{h,0,0},{0,h,0},{0,0,h}};
        for(int k=0;k<3;++k){
            const ElectromagneticField plus=partner(first,x0+e[k]), minus=partner(first,x0-e[k]);
            dE[k]=(plus.electric-minus.electric)*(1.0/(2*h));
            dB[k]=(plus.magnetic-minus.magnetic)*(1.0/(2*h)); }
        const Vec3 curlE{dE[1].z-dE[2].y,dE[2].x-dE[0].z,dE[0].y-dE[1].x};
        const Vec3 pGradB=dB[0]*p.x+dB[1]*p.y+dB[2]*p.z;
        const ElectromagneticField here=partner(first,x0);
        const Vec3 pRate=cross(omega,p);
        k1+=cross(p,curlE)*-1.0+cross(pRate,here.magnetic)+cross(v,pGradB);
        const Vec3 hidLab=cross(mLab,here.electric)*(1.0/(c*c));
        const Vec3 hidCov=cross(m0,here.electric+cross(v,here.magnetic))*(1.0/(c*c));
        k2+=cross(omega,hidCov-hidLab)*-1.0;
        k2a+=cross(omega,cross(m0,cross(v,here.magnetic)))*(-1.0/(c*c));
        k2b+=cross(omega,cross(m0-mLab,here.electric))*(-1.0/(c*c));
    }
    const auto report=[&](const char* name,const Vec3& k){
        std::printf("  candidate %-6s %+.6e %+.6e %+.6e  |Delta-K|/|Delta| %.4f  K.x/Delta.x %.4f\n",
            name,k.x,k.y,k.z,(delta-k).norm()/std::max(delta.norm(),1e-300),
            delta.x!=0.0?k.x/delta.x:0.0); };
    std::printf("  Delta = expected - F1+F2   %+.6e %+.6e %+.6e\n",delta.x,delta.y,delta.z);
    report("K1",k1); report("K2",k2); report("K1+K2",k1+k2);
    report("K2a",k2a); report("K2b",k2b);
    // the Thomas back-reaction's own net force (audit 567): L_T depends on
    // r_p - r_o, so F_p + F_o = -dA_T/dt, not zero
    {
        const MutualForces th=thomasBackReactionForces(present,
            thomasBmtDipoleDerivatives(present,history));
        const Vec3 tn=th.first+th.second;
        std::printf("  Thomas back-reaction net  %+.6e %+.6e %+.6e  |.|/|F1+F2| %.4f\n",
            tn.x,tn.y,tn.z,tn.norm()/std::max(net.norm(),1e-300));
    }
    std::printf("  |m_lab|/|m_0| - 1  %.6e   v.B check |v||B_partner|/(|E|) %.6e\n",
        present.firstDipole.norm()/present.firstProperDipole.norm()-1.0,
        present.firstVelocity.norm()*partner(true,present.firstPosition).magnetic.norm()
            /partner(true,present.firstPosition).electric.norm());
}
