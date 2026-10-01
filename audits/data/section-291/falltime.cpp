// Audit 291c: czemu zmierzony czas spadku przekracza analityczny o 3.5e-05?
// PREDYKCJA ZAREJESTROWANA PRZED URUCHOMIENIEM: to zmiekczanie Plummera.
// Opoznienie zbiera sie wewnatrz r <~ r_podlogi, gdzie sila jest slabsza,
// wiec skaluje sie jak r_podlogi^1.5.  Obnizenie podlogi 4x ma zmniejszyc
// nadwyzke 8x: 3.5e-05 -> 4.4e-06.  Jesli nadwyzka nie drgnie, to nie
// zmiekczanie, a wtedy kandydatem jest blad calkowania.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double rStar=comptonBarrierRadius;
    const double a=0.5*A;
    const double exact=pi*std::sqrt(a*a*a/kk);
    // reszta drogi pod odcieciem, odjeta analitycznie
    const double rc=nuclearCutoff;
    const double E=std::acos(1.0-rc/a);
    const double tail=std::sqrt(a*a*a/kk)*(E-std::sin(E));
    std::printf("# podloga=%.5f r*, analityczny polokres=%.6e s,"
        " ogon pod odcieciem=%.3e s (%.2e polokresu)\n",
        separationFloor()/rStar,exact,tail,tail/exact);
    const double w1=secondMass/(firstMass+secondMass);
    const double w2=-firstMass/(firstMass+secondMass);
    State start;
    start.firstPosition=Vec3{A,0,0}*w1;
    start.secondPosition=Vec3{A,0,0}*w2;
    SimulationOptions options;
    options.collectFrames=false;
    options.radiatedEnergyBookkeeping=false;
    const MechanicalTrajectoryResult run=runMechanicalTrajectory(
        start,1.2*exact,nuclearCutoff,options,
        ChargeRadiationReactionModel::disabled);
    const double measured=run.elapsedTime+tail;  // ekstrapolacja do zera
    std::printf("czas=%.9e s  +ogon=%.9e s  nadwyzka=%+.4e\n",
        run.elapsedTime,measured,measured/exact-1.0);
}
