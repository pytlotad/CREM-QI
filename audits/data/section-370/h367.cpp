#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 367: first-photon time from a Langer start, for a chosen pair, level and J0 (env PAIR, LEVEL, J0); prints CREM_REACH.
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  gInitialPrincipalLevel=std::atoi(std::getenv("LEVEL")); gInitialAngularMomentumFraction=std::atof(std::getenv("J0"));
  gMicrocanonicalStart=false; gSpinQuantization=true;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    std::printf("BEGIN %d\n",i); std::fflush(stdout);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+i),1,3000.0);
    std::printf("END idx %d stop=%d photons=%lld lifetime=%.6e\n",i,(int)e.stopCause,(long long)e.emittedPhotonCount,e.lifetimeSeconds);
    std::fflush(stdout); }
}
