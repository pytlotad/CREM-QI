// Audit 249: how does 248's irreversibility floor accumulate?
//
// 248e left the decisive number unmeasured: 1.0e-08 per quarter
// orbit is 1.6e-05 over a full cascade if it random-walks and
// 2.5e-02 if it accumulates linearly.  Measured here by running the
// same forward/reverse test over 1, 2, 4, 8, 16 and 32 quarter
// orbits at a fixed tolerance, and fitting the exponent.
//   slope 0.5 -> random walk, harmless
//   slope 1.0 -> linear, percent-level over a cascade
//   slope 1.5 -> worse than either
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    std::printf("%10s %14s %12s\n","cwiartki","|dr|/r","nachylenie");
    double prevErr=0.0; double prevQ=0.0;
    for(int q:{1,2,4,8,16,32}){
        const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}*firstMagneticMoment;
        State start;
        const double w1=secondMass/(firstMass+secondMass);
        const double w2=-firstMass/(firstMass+secondMass);
        const double v=std::sqrt(k/a);
        start.firstPosition=Vec3{a,0,0}*w1;
        start.secondPosition=Vec3{a,0,0}*w2;
        start.firstVelocity=Vec3{0,v,0}*w1;
        start.secondVelocity=Vec3{0,v,0}*w2;
        start.firstDipole=m1;
        start.secondDipole=m1*(secondMagneticMoment/firstMagneticMoment);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-12; acc.maximumDepth=30;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        acc.useRetardedExternalForces=false;
        const int n=200*q;
        const double dt=0.25*q*period/n;
        ClassicalTrajectoryEngine fwd(start,acc);
        State s=start; bool ok=true;
        for(int i=0;i<n&&ok;++i) ok=fwd.advance(s,dt);
        if(!ok){ std::printf("%10d  failed\n",q); continue; }
        State back=s;
        back.firstVelocity=s.firstVelocity*-1.0;
        back.secondVelocity=s.secondVelocity*-1.0;
        back.firstDipole=s.firstDipole*-1.0;
        back.secondDipole=s.secondDipole*-1.0;
        back.firstProperDipole=Vec3{};
        back.secondProperDipole=Vec3{};
        ClassicalTrajectoryEngine rev(back,acc);
        State r=back;
        for(int i=0;i<n&&ok;++i) ok=rev.advance(r,dt);
        if(!ok){ std::printf("%10d  failed wstecz\n",q); continue; }
        const Vec3 dr=(r.firstPosition-r.secondPosition)
                     -(start.firstPosition-start.secondPosition);
        const double rel=dr.norm()
            /(start.firstPosition-start.secondPosition).norm();
        const double slope=prevErr>0.0
            ?std::log(rel/prevErr)/std::log(q/prevQ):0.0;
        std::printf("%10d %14.6e %12.3f\n",q,rel,slope);
        prevErr=rel; prevQ=q;
    }
    std::printf("# 0.5 blazenie losowe, 1.0 liniowo, 1.5 gorzej\n");
}
