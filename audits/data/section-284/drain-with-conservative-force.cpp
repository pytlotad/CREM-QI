// Audit 283: does the stochastic emission mode double-count the
// radiated energy?
//
// crem code at electrodynamics.hpp:2032 states that
// stochasticElectricDipole carries ZERO continuous reaction force,
// "same as disabled: its whole point is that nothing drags the orbit
// between photons".  But retardedExternalForces is unchanged, and
// 281b measured what it drains: the acceleration term -0.795 of the
// radiated energy plus the velocity term's non-conservative part
// +0.263, i.e. about -0.53 continuously.  Photons are then drawn
// against quantizedPower, the full analytic sum.
//
// If both happen, the emission budget counts about half the radiated
// energy twice.  Measured here as the ratio of the mechanical energy
// the orbit actually loses to the energy the field actually carries,
// with the continuous reaction off -- which is the same force
// configuration the stochastic mode runs in between photons.
//
// Prediction, written before the run: about 0.53 in energy, matching
// the 0.5360 that 277d/280h measured in angular momentum.
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
    std::printf("%26s %14s %14s %10s\n","konfiguracja","dE_mech",
        "E_wypromien.","iloraz");
    struct Case{const char* name;ChargeRadiationReactionModel model;
                bool retarded;};
    for(const Case& cs:{
            Case{"stochastic, retardowana",
                 ChargeRadiationReactionModel::stochasticElectricDipole,true},
            Case{"stochastic, ZACHOWAWCZA",
                 ChargeRadiationReactionModel::stochasticElectricDipole,false},
            Case{"LL, retardowana",
                 ChargeRadiationReactionModel::individualLandauLifshitz,true},
            Case{"LL, ZACHOWAWCZA",
                 ChargeRadiationReactionModel::individualLandauLifshitz,
                 false}}){
        const State start=osculatingPeriapsisState(el,kk,m1,
            m1*(secondMagneticMoment/firstMagneticMoment),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
        acc.reactionModel=cs.model;
        acc.useRetardedExternalForces=cs.retarded;
        acc.computeOutwardFlux=true;
        acc.farFieldSampling.directionCount=194;
        acc.farFieldSampling.controlRadius=1.0e6*bohrRadius;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<200&&ok;++i) ok=engine.advance(st,period/400);
        if(!ok){ std::printf("%26s przerwane\n",cs.name); continue; }
        const double weight=retardedDipoleSectorWeight(st);
        const double mechanicalFirst=conservativeParticleEnergy(st,weight);
        const double radiatedFirst=st.radiatedEnergy;
        for(int i=0;i<400&&ok;++i) ok=engine.advance(st,period/400);
        if(!ok){ std::printf("%26s przerwane\n",cs.name); continue; }
        const double weightLast=retardedDipoleSectorWeight(st);
        const double deltaMechanical=
            conservativeParticleEnergy(st,weightLast)-mechanicalFirst;
        const double radiated=st.radiatedEnergy-radiatedFirst;
        std::printf("%26s %14.6e %14.6e %10.4f\n",cs.name,
            deltaMechanical,radiated,
            radiated!=0.0?-deltaMechanical/radiated:0.0);
    }
    std::printf("# iloraz 1 = orbita traci dokladnie tyle, ile pole unosi\n");
    std::printf("# iloraz 0.53 przy disabled/stochastic = reszta musi "
        "przyjsc od fotonu,\n#   ale foton jest losowany na pelna moc "
        "-> podwojne liczenie\n");
}
