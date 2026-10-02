// Audit 311: odtworzenie "--diagnose" z README dla obu kanalow, z/bez
// kwantyzacji spinu i z/bez sily dipolowej (CREM_NO_DIPOLE_FORCE czytane
// przy starcie procesu, wiec kazda kombinacja to osobny proces).
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
    const int level=argc>1?std::atoi(argv[1]):2;
    gInitialPrincipalLevel=level;
    if(std::getenv("CREM_PROBE_SPIN_QUANT")) gSpinQuantization=true;
    const double eV=1.602176634e-19;
    for(int ph=1;ph<=2;++ph){
        SimulationOptions o;                  // jak --diagnose: okno zjawiska
        const SimulationResult r=simulate(42,ph,o);
        if(r.frames.empty()){ std::printf("%s brak klatek\n",ph==1?"para":"orto"); continue; }
        const auto& f0=r.frames.front(); const auto& f1=r.frames.back();
        const Vec3 m=f0.firstDipole+f0.secondDipole;
        const double mu=f0.firstDipole.norm();
        const double cosm=dot(f0.firstDipole,f0.secondDipole)/(f0.firstDipole.norm()*f0.secondDipole.norm());
        std::printf("%s  E_rad=%.7e eV  |P_rad|=%.3e kg m/s  |mu1+mu2|/mu=%.6f  cos=%+.6f  t=%.4e s  wynik=%d\n",
            ph==1?"para":"orto",(f1.radiatedEnergy-f0.radiatedEnergy)/eV,
            (f1.radiatedMomentum-f0.radiatedMomentum).norm(),m.norm()/mu,cosm,
            r.elapsedTime,static_cast<int>(r.outcome));
    }
}
