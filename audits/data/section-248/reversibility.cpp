// Audit 248: is the PRODUCTION integrator time-reversible?
//
// The objection: the step is explicit predictor-corrector, the
// magnetic and Darwin forces depend on v, and an explicit kick with
// velocity-dependent forces is neither symplectic nor reversible --
// it wants a Boris-type implicit treatment.  Reversibility is
// checked in validation with relativisticBorisPush on a synthetic
// uniform field (maxwell_validation.hpp:1220-1225), not on the
// production path.
//
// The test: integrate N steps forward, apply time reversal
// (v -> -v, mu -> -mu, both T-odd), integrate N steps back, and
// compare with the start.  Retardation is switched OFF and the
// reaction disabled, because a retarded system is not reversible for
// physical reasons and would confound the numerical question.  Under
// those settings the dynamics IS reversible, so whatever is left is
// the scheme.
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
    std::printf("%10s %10s %14s %14s %10s\n",
                "momenty","tolerancja","|dr|/r","|dv|/v","spadek");
    for(int moments=1;moments>=0;--moments){
        const double scale=moments?1.0:1.0e-6;
        double previous=0.0;
        for(double tol:{1e-9,1e-10,1e-11,1e-12,1e-13,1e-14}){
            const int n=200;
            const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}
                *(firstMagneticMoment*scale);
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
            acc.relativeTolerance=tol; acc.maximumDepth=30;
            acc.reactionModel=ChargeRadiationReactionModel::disabled;
            acc.computeOutwardFlux=false;
            acc.useRetardedExternalForces=false;
            const double dt=0.25*period/n;
            ClassicalTrajectoryEngine forward(start,acc);
            State s=start; bool ok=true;
            for(int i=0;i<n&&ok;++i) ok=forward.advance(s,dt);
            if(!ok){ std::printf("%10s %10.0e  failed\n",
                     moments?"wlaczone":"wylaczone",tol); continue; }
            // Time reversal: v and mu are both T-odd.
            State back=s;
            back.firstVelocity=s.firstVelocity*-1.0;
            back.secondVelocity=s.secondVelocity*-1.0;
            back.firstDipole=s.firstDipole*-1.0;
            back.secondDipole=s.secondDipole*-1.0;
            back.firstProperDipole=Vec3{};
            back.secondProperDipole=Vec3{};
            ClassicalTrajectoryEngine reverse(back,acc);
            State r=back;
            for(int i=0;i<n&&ok;++i) ok=reverse.advance(r,dt);
            if(!ok){ std::printf("%10s %10.0e  failed wstecz\n",
                     moments?"wlaczone":"wylaczone",tol); continue; }
            const Vec3 dr=(r.firstPosition-r.secondPosition)
                         -(start.firstPosition-start.secondPosition);
            const Vec3 dv=(r.firstVelocity*-1.0-start.firstVelocity);
            const double rel=dr.norm()
                /(start.firstPosition-start.secondPosition).norm();
            const double relv=dv.norm()/start.firstVelocity.norm();
            const double order=previous>0.0
                ?std::log(previous/rel)/std::log(2.0):0.0;
            std::printf("%10s %10.0e %14.6e %14.6e %10.3f\n",
                        moments?"wlaczone":"wylaczone",tol,rel,relv,
                        previous>0.0?previous/rel:0.0);
            previous=rel;
        }
    }
    std::printf("# schemat odwracalny dalby blad na poziomie zaokraglenia\n");
}
