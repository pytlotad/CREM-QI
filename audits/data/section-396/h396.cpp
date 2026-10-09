#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 396: hydrogen (PAIR) from a Langer start LEVEL/J0, photons by type; argv = seed indices.  Env as for h367 (section 372).
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  gInitialPrincipalLevel=std::atoi(std::getenv("LEVEL")); gInitialAngularMomentumFraction=std::atof(std::getenv("J0"));
  gMicrocanonicalStart=false; gSpinQuantization=true;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+i),1,3000.0);
    std::printf("END idx %d stop=%d photons=%lld E2=%lld M1=%lld lifetime=%.6e L_end/hbar=%.6f\n",i,(int)e.stopCause,
      (long long)e.emittedPhotonCount,(long long)e.quadrupolePhotonCount,(long long)e.magneticPhotonCount,e.lifetimeSeconds,e.terminalAngularMomentum);
    std::fflush(stdout); }
}
