// Audit 242: F^mu u_mu = 0 on trajectory states, and the quantity
// that test cannot see.
//
// The objection: the four-force from a scalar coupling
// U0 = -mu0.B0 has a component along u^mu, so the rest mass is not
// constant, m(tau) = m0 + U0/c^2, while the code holds m fixed.
//
// The proposed probe passes trivially.  fourForce() builds the time
// component FROM the space component as gamma*(F.v)/c, so
// F^mu u_mu = gamma^2 (F.v - F.v) = 0 identically.  Orthogonality is
// imposed, not measured, and the longitudinal piece is discarded
// before anything can test it.  Both are reported: the trivial test
// to show it is trivial, and the discarded piece, which is the real
// quantity.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    std::printf("%9s %14s %14s %14s %14s\n",
                "r/a_pair","F^mu u_mu","/(|F||u|)","U0/(m c^2)",
                "dU0/U0 orbita");
    for(double f:{0.25,0.1,0.05,0.02,0.01}){
        const double a=f*A;
        const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}*firstMagneticMoment;
        const State start=osculatingPeriapsisState(
            el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start; bool ok=true;
        const int n=400;
        double u0lo=1e300,u0hi=-1e300,worst=0.0,worstRel=0.0,u0mid=0.0;
        for(int i=0;i<n&&ok;++i){
            ok=engine.advance(s,period/n);
            const MutualForces F=retardedExternalForces(s,engine.history());
            // The code's own four-force, and its Minkowski product
            // with the four-velocity.
            const FourVector K=fourForce(s.firstVelocity,F.first);
            const FourVector U=fourVelocity(s.firstVelocity);
            const double prod=minkowskiDot(K,U);
            const double scale=std::sqrt(K.time*K.time
                +K.space.squaredNorm())
                *std::sqrt(U.time*U.time+U.space.squaredNorm());
            worst=std::max(worst,std::abs(prod));
            if(scale>0.0) worstRel=std::max(worstRel,std::abs(prod)/scale);
            // The rest-frame coupling U0 = -mu0.B0, the quantity whose
            // variation the constant rest mass discards.
            const ElectromagneticField lab=fieldFromOtherParticleAt(
                s.firstPosition,s.time,s,engine.history(),true);
            const Vec3 v=s.firstVelocity;
            const double g=gamma(v);
            const double v2=v.squaredNorm();
            const Vec3 B0=v2>0.0
                ?(lab.magnetic-cross(v,lab.electric)*(1.0/(c*c)))*g
                  -v*((g-1.0)*dot(lab.magnetic,v)/v2)
                :lab.magnetic;
            const double U0=-dot(s.firstProperDipole,B0);
            u0lo=std::min(u0lo,U0); u0hi=std::max(u0hi,U0);
            if(i==n/2) u0mid=U0;
        }
        if(!ok){ std::printf("%9.3f failed\n",f); continue; }
        std::printf("%9.3f %14.6e %14.6e %14.6e %14.6e\n",
            f,worst,worstRel,u0mid/(electron.mass*c*c),
            (u0hi-u0lo)/(electron.mass*c*c));
    }
    std::printf("# kolumna 2-3: test zaproponowany, zero z konstrukcji\n");
    std::printf("# kolumna 4-5: U0/(m c^2), czyli dm/m ktore kod pomija\n");
}
