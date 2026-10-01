// Audit 291, czesc druga: pytamy SILNIK o jego minimumSeparation.
//
// Moja sonda wolala ClassicalTrajectoryEngine wprost i probkowala na
// STALEJ siatce, wiec dno przeskakiwala -- przy 0.084 r* predkosc to
// 0.83 c, a jeden krok pokrywal dwa razy wiecej drogi niz cale dno.
// Gorzej: omijala TERMINACJE.  runMechanicalTrajectory konczy przebieg
// na trajectoryCutoff (crem_trajectory.hpp:695) z wynikiem
// ReachedCutoff, a produkcja przekazuje tam nuclearCutoff, czyli
// 10 fm = 0.0517 r*.  To jest ponizej moich probkowanych 0.0838 r*,
// wiec mierzylem przebieg, ktorego produkcja nigdy nie wykonuje.
//
// Tu minimum czyta sie z MechanicalTrajectoryResult, gdzie jest
// zapisywane na KAZDYM zatwierdzonym kroku adaptacyjnym, a adaptacja
// zageszcza krok przy przelocie.  Skan odciecia sprawdza, czy minimum
// je sledzi.
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
    const double halfPeriod=pi*std::sqrt(std::pow(0.5*A,3)/kk);
    std::printf("# nuclearCutoff=%.4e m = %.5f r*,"
        " separationFloor=%.4e m = %.5f r*\n",
        nuclearCutoff,nuclearCutoff/rStar,
        separationFloor(),separationFloor()/rStar);
    std::printf("%14s %12s %13s %9s %13s\n","odciecie/r*","min/r*",
        "min/odciecie","wynik","czas [s]");
    for(double factor:{1.0,0.3,0.1,0.03}){
        const double cutoff=nuclearCutoff*factor;
        const double w1=secondMass/(firstMass+secondMass);
        const double w2=-firstMass/(firstMass+secondMass);
        State start;
        start.firstPosition=Vec3{A,0,0}*w1;
        start.secondPosition=Vec3{A,0,0}*w2;
        SimulationOptions options;
        options.collectFrames=false;
        options.radiatedEnergyBookkeeping=false;
        const MechanicalTrajectoryResult run=runMechanicalTrajectory(
            start,1.2*halfPeriod,cutoff,options,
            ChargeRadiationReactionModel::disabled);
        std::printf("%14.5f %12.5f %13.5f %9d %13.6e\n",
            cutoff/rStar,run.minimumSeparation/rStar,
            run.minimumSeparation/cutoff,
            static_cast<int>(run.outcome),run.elapsedTime);
    }
    std::printf("# min/odciecie = 1 oznacza, ze minimum wyznacza"
        " POWIERZCHNIA ODCIECIA\n");
}
