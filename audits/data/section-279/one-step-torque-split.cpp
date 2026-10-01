// Audit 279: separate the mutual and self torques WITHIN ONE STEP.
//
// 278e's 53/50 split put a number measured by SUBTRACTING TWO RUNS
// (reaction on minus reaction off) next to a number read directly out
// of one run.  That is not a decomposition: switching the reaction
// off does not remove one term from a fixed trajectory, it changes
// the trajectory, and with it the bound-field content.  278g recorded
// that as a limit; it is really an objection to the section's main
// claim.
//
// Here all three torques come out of ONE state and ONE call:
//   retardedExternalForces  -- the production mutual force assembly,
//   particleMultipoleRadiation(...).chargeReaction -- the self
//     reaction force, computed in the same call from the same state,
//   ...outwardFlux.angularMomentum -- the flux, likewise.
// Nothing is subtracted between runs.  Averaged over a whole orbit,
// because the instantaneous torques oscillate and 278c showed half an
// orbit is transient-dominated.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const OsculatingElements el{-kk/(2.0*A),std::sqrt(kk*A)};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
    const State start=osculatingPeriapsisState(el,kk,m1,
        m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    FarFieldSampling sampling;
    sampling.directionCount=194;
    sampling.controlRadius=1.0e6*bohrRadius;
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
    acc.reactionModel=
        ChargeRadiationReactionModel::individualLandauLifshitz;
    acc.computeOutwardFlux=true;
    acc.farFieldSampling=sampling;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start; bool ok=true;
    // Let the startup transient pass: one orbit in before sampling.
    for(int i=0;i<200&&ok;++i) ok=engine.advance(st,period/200);
    if(!ok){ std::printf("nieudane\n"); return 1; }
    Vec3 mutualSum{},selfSum{},dipoleSum{},fluxSum{};
    int samples=0;
    for(int i=0;i<200&&ok;++i){
        ok=engine.advance(st,period/200);
        if(!ok) break;
        const StateHistory& h=engine.history();
        const MutualForces mutual=retardedExternalForces(st,h);
        const ParticleMultipoleRadiation radiation=
            particleMultipoleRadiation(st,mutual,h,true,
                ChargeRadiationReactionModel::individualLandauLifshitz,
                true,sampling);
        mutualSum=mutualSum
            +cross(st.firstPosition,mutual.first)
            +cross(st.secondPosition,mutual.second);
        selfSum=selfSum
            +cross(st.firstPosition,radiation.chargeReaction.first)
            +cross(st.secondPosition,radiation.chargeReaction.second);
        dipoleSum=dipoleSum
            +radiation.firstDipoleTorque+radiation.secondDipoleTorque;
        fluxSum=fluxSum+radiation.outwardFlux.angularMomentum;
        ++samples;
    }
    if(samples<1){ std::printf("brak probek\n"); return 1; }
    const double n=samples;
    const Vec3 mutualTorque=mutualSum*(1.0/n);
    const Vec3 selfTorque=selfSum*(1.0/n);
    const Vec3 dipoleTorque=dipoleSum*(1.0/n);
    const Vec3 flux=fluxSum*(1.0/n);
    std::printf("# jeden obieg, %d probek, wszystko z jednego wywolania\n",
        samples);
    std::printf("%26s %14s %10s\n","czlon","z-skladowa","/|strumien|");
    const double scale=std::abs(flux.z)>0.0?std::abs(flux.z):1.0;
    std::printf("%26s %14.6e %10.4f\n","moment sily WZAJEMNEJ",
        mutualTorque.z,mutualTorque.z/scale);
    std::printf("%26s %14.6e %10.4f\n","moment sily WLASNEJ",
        selfTorque.z,selfTorque.z/scale);
    std::printf("%26s %14.6e %10.4f\n","moment sily DIPOLOWEJ",
        dipoleTorque.z,dipoleTorque.z/scale);
    std::printf("%26s %14.6e %10.4f\n","strumien pola",flux.z,
        flux.z/scale);
    const double residual=mutualTorque.z+selfTorque.z+flux.z;
    std::printf("%26s %14.6e %10.4f\n","RESIDUUM (wzaj+wl+str)",
        residual,residual/scale);
    std::printf("%26s %14.6e %10.4f\n","z dipolowym",
        residual+dipoleTorque.z,(residual+dipoleTorque.z)/scale);
}
