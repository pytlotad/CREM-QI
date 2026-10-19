#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 512 (from 490): photon as a fast event in the current configuration (imports 20, Larmor, split transport); env BEXT [T], random direction per run.
// Audit 490: photon emission as a fast event -- cascade from a non-quantized start, ledger before/after each photon.
// env PAIR, NC (rung n), argv: run indices.  n0 drawn with density ~ n0^2 on [NC - 1/2, NC + 1/2] (phase-space measure of
// the energy shell), L from the engine's microcanonical start (L^2 uniform on [0, (n0 hbar)^2]).  The lock's trace
// (CREM_DEBUG_WAVE_LOCK) goes to stderr; START lines mark each run there.
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  const double nc=std::atof(std::getenv("NC"));
  gMicrocanonicalStart=true; gSpinQuantization=false;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    std::mt19937_64 r(512000ULL+i); std::uniform_real_distribution<double> u(0,1);
    const double lo=nc-0.5, hi=nc+0.5, n0=std::cbrt(lo*lo*lo+u(r)*(hi*hi*hi-lo*lo*lo));
    gInitialSeparationOverride=n0*n0*pairBohrRadius(activePair);
    { std::uniform_real_distribution<double> ub(0,1); const double cz=2*ub(r)-1, ph=2*pi*ub(r), sz=std::sqrt(1-cz*cz); const double B=std::getenv("BEXT")?std::atof(std::getenv("BEXT")):0.0;
      gExternalMagneticField=Vec3{sz*std::cos(ph),sz*std::sin(ph),cz}*B; std::fprintf(stderr,"BFIELD idx %d %.12e %.12e %.12e\n",i,gExternalMagneticField.x,gExternalMagneticField.y,gExternalMagneticField.z); }
    std::fprintf(stderr,"START idx %d n0 %.9f\n",i,n0); std::fflush(stderr);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+i),1,600.0);
    std::fprintf(stderr,"END idx %d stop %d photons %lld\n",i,(int)e.stopCause,(long long)e.emittedPhotonCount); std::fflush(stderr);
  }
}
