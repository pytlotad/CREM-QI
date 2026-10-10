#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 396: hydrogen (PAIR) from a Langer start LEVEL/J0, photons by type (audit 397: PHEN 1 p-Ps, 2 o-Ps; BUDGET s); argv = seed indices.  Env as for h367 (section 372).
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  gInitialPrincipalLevel=std::atoi(std::getenv("LEVEL")); gInitialAngularMomentumFraction=std::atof(std::getenv("J0"));
  gMicrocanonicalStart=false; gSpinQuantization=true;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+i),std::atoi(std::getenv("PHEN")),std::atof(std::getenv("BUDGET")));
    std::printf("END idx %d stop=%d photons=%lld E2=%lld M1=%lld refusedM1=%lld flipE=%.6e lifetime=%.6e L_end/hbar=%.6f\n",i,(int)e.stopCause,
      (long long)e.emittedPhotonCount,(long long)e.quadrupolePhotonCount,(long long)e.magneticPhotonCount,(long long)e.refusedMagneticFlips,e.magneticFlipEnergyJoules,e.lifetimeSeconds,e.terminalAngularMomentum);
    std::fflush(stdout); }
}
