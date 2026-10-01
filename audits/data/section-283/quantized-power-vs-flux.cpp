// Audit 283, brakujace ogniwo: ile wynosi quantizedPower wzgledem
// ZMIERZONEGO strumienia?
//
// Bez tego nie wolno mowic o podwojnym liczeniu.  Orbita traci 53.05
// procenta strumienia ciagle (zmierzone), a foton jest losowany
// przeciw quantizedPower = leadingElectricDipolePower +
// magneticDipoleFlux.energy + electricQuadrupolePower, czyli sumie
// ANALITYCZNEJ.  Jesli ta suma rowna sie strumieniowi, budzet wynosi
// 153 procent i podwojne liczenie jest realne.  Jesli wynosi 0.47
// strumienia, budzet domyka sie dokladnie i problemu nie ma.
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
        ChargeRadiationReactionModel::stochasticElectricDipole;
    acc.computeOutwardFlux=true;
    acc.farFieldSampling=sampling;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start; bool ok=true;
    for(int i=0;i<200&&ok;++i) ok=engine.advance(st,period/400);
    if(!ok){ std::printf("nieudane\n"); return 1; }
    double quantized=0.0,leading=0.0,dipole=0.0,quadrupole=0.0,flux=0.0;
    int n=0;
    for(int i=0;i<200&&ok;++i){
        ok=engine.advance(st,period/400);
        if(!ok) break;
        const StateHistory& h=engine.history();
        const MutualForces mutual=retardedExternalForces(st,h);
        const ParticleMultipoleRadiation r=
            particleMultipoleRadiation(st,mutual,h,true,
                ChargeRadiationReactionModel::stochasticElectricDipole,
                true,sampling);
        leading+=r.leadingElectricDipolePower;
        dipole+=r.magneticDipoleFlux.energy;
        quadrupole+=r.electricQuadrupolePower;
        quantized+=r.leadingElectricDipolePower
            +r.magneticDipoleFlux.energy+r.electricQuadrupolePower;
        flux+=r.outwardFlux.energy;
        ++n;
    }
    if(n<1){ std::printf("brak probek\n"); return 1; }
    const double inv=1.0/n;
    std::printf("# %d probek, moce usrednione [W]\n",n);
    std::printf("%34s %14.6e %9.4f\n","leadingElectricDipolePower",
        leading*inv,leading/flux);
    std::printf("%34s %14.6e %9.4f\n","magneticDipoleFlux.energy",
        dipole*inv,dipole/flux);
    std::printf("%34s %14.6e %9.4f\n","electricQuadrupolePower",
        quadrupole*inv,quadrupole/flux);
    std::printf("%34s %14.6e %9.4f\n","quantizedPower (suma)",
        quantized*inv,quantized/flux);
    std::printf("%34s %14.6e %9.4f\n","outwardFlux.energy (zmierzony)",
        flux*inv,1.0);
    std::printf("\n# budzet trybu stochastycznego:\n");
    std::printf("%34s %9.4f\n","ciagly dren (zmierzony)",0.5305);
    std::printf("%34s %9.4f\n","fotony (quantizedPower)",quantized/flux);
    std::printf("%34s %9.4f\n","RAZEM",0.5305+quantized/flux);
}
