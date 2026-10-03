// Audit 316: IMPORT jawnie oznaczony -- "anihilacja przy pierwszym wejsciu w
// r <= r*" (terminalSeparation = r*), dynamika zachowawcza (reakcja
// wylaczona), start n = 1, L = J0 hbar (CREM_INITIAL_ANGULAR_MOMENTUM).
// Czas zycia = czas pierwszego wejscia; brak wejscia w oknie = cenzura.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
    const std::uint64_t seed=std::strtoull(argv[1],nullptr,10);
    const int ph=std::atoi(argv[2]);
    const double orbits=argc>3?std::atof(argv[3]):20.0;
    gRadiationReactionModel=ChargeRadiationReactionModel::disabled;
    if(std::getenv("CREM_PROBE_SPIN_QUANT")) gSpinQuantization=true;
    SimulationOptions o; o.collectFrames=true; o.frameCount=2;
    o.radiatedEnergyBookkeeping=false;
    const double period=3.0397e-16;
    o.observationTime=(orbits+0.05)*period;
    o.terminalSeparation=comptonBarrierRadius;      // regula anihilacji
    const SimulationResult r=simulate(seed,ph,o);
    const auto& f=r.frames.front();
    const double c12=dot(f.firstDipole,f.secondDipole)/(f.firstDipole.norm()*f.secondDipole.norm());
    // 317: waga 2 gamma w chwili zdarzenia (ostatnia klatka = klatka
    // przeciecia r* przy outcome 0; przy cenzurze = koniec okna).
    const auto& g=r.frames.back();
    const double mu=0.5*(g.firstDipole.norm()+g.secondDipole.norm());
    const double w=std::pow((g.firstDipole+g.secondDipole).norm()/(2.0*mu),2);
    const double c12e=dot(g.firstDipole,g.secondDipole)/(g.firstDipole.norm()*g.secondDipole.norm());
    const double w0=std::pow((f.firstDipole+f.secondDipole).norm()/(f.firstDipole.norm()+f.secondDipole.norm()),2);
    std::printf("%llu %d %d %.6e %.4f %.6f %+.5f %.6f %.6f %+.5f\n",(unsigned long long)seed,ph,
        static_cast<int>(r.outcome),r.elapsedTime,r.elapsedTime/period,
        r.minimumSeparation/comptonBarrierRadius,c12,w0,w,c12e);
}
