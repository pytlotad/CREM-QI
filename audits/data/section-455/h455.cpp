#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 424: path-memory wave lock (CREM_WAVE_LOCK) from a non-integer start.
// env PAIR, START_N (real n; a = n^2 a_pair through gInitialSeparationOverride),
// L_HBAR (orbital L in hbar; J0 = L/n).  Prints CREM_REACH per photon and a
// summary per run.
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  const double n0=std::atof(std::getenv("START_N")), l=std::atof(std::getenv("L_HBAR"));
  gInitialSeparationOverride=n0*n0*pairBohrRadius(activePair);
  gInitialAngularMomentumFraction=l/n0;
  gMicrocanonicalStart=false; gSpinQuantization=true;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    std::printf("BEGIN %d\n",i); std::fflush(stdout);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+i),std::getenv("PHEN")?std::atoi(std::getenv("PHEN")):1,3000.0);
    std::printf("END idx %d stop=%d photons=%lld lifetime=%.6e lockWork_eV=%.9e lockShift=%.9e lockN=%llu Lwave=%.6f Lspins=%.6f fallbacks=%llu\n",i,(int)e.stopCause,
      (long long)e.emittedPhotonCount,e.lifetimeSeconds,e.waveLockWorkJoules/1.602176634e-19,e.waveLockLargestShift,e.waveLockApplications,e.waveLockAngularHbar,e.waveLockAngularToSpinsHbar,e.waveLockSpinFallbacks);
    std::fflush(stdout); }
}
