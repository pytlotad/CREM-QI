#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 500: Langer 2p start (n = 2, L = 1.5 hbar) -- time to the first photon (2 -> 1) from the emission ledger; env PAIR; argv: run indices.
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  gSpinQuantization=false; gInitialAngularMomentumFraction=0.75;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    gInitialSeparationOverride=4.0*pairBohrRadius(activePair);
    std::fprintf(stderr,"START idx %d\n",i); std::fflush(stderr);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(500ULL+i),1,1200.0);
    std::fprintf(stderr,"END idx %d stop %d photons %lld\n",i,(int)e.stopCause,(long long)e.emittedPhotonCount); std::fflush(stderr);
  }
}
