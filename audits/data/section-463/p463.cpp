#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 463: frozen wave lock (CREM_WAVE_LOCK defaults of 424/453) from phase-space-uniform starts.
// env PAIR, NC (rung n), argv: run indices.  n0 drawn with density ~ n0^2 on [NC - 1/2, NC + 1/2] (phase-space measure of
// the energy shell), L from the engine's microcanonical start (L^2 uniform on [0, (n0 hbar)^2]).  The lock's trace
// (CREM_DEBUG_WAVE_LOCK) goes to stderr; START lines mark each run there.
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  const double nc=std::atof(std::getenv("NC"));
  gMicrocanonicalStart=true; gSpinQuantization=true;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    std::mt19937_64 r(463000ULL+i); std::uniform_real_distribution<double> u(0,1);
    const double lo=nc-0.5, hi=nc+0.5, n0=std::cbrt(lo*lo*lo+u(r)*(hi*hi*hi-lo*lo*lo));
    gInitialSeparationOverride=n0*n0*pairBohrRadius(activePair);
    std::fprintf(stderr,"START idx %d n0 %.9f\n",i,n0); std::fflush(stderr);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+i),1,600.0);
    std::fprintf(stderr,"END idx %d stop %d photons %lld\n",i,(int)e.stopCause,(long long)e.emittedPhotonCount); std::fflush(stderr);
  }
}
