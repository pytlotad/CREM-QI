#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 503 (from 500): Langer 2p start in an external field BEXT [T] along a random direction per run (seeded) -- photon energy,
// direction, helicity and angular momentum relative to B from the emission ledger.  env PAIR, BEXT; argv: run indices.
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  gSpinQuantization=false; gInitialAngularMomentumFraction=0.75;
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    gInitialSeparationOverride=4.0*pairBohrRadius(activePair);
    { std::mt19937_64 r(503000ULL+i); std::uniform_real_distribution<double> u(0,1); const double cz=2*u(r)-1, ph=2*pi*u(r), sz=std::sqrt(1-cz*cz);
      const double B=std::getenv("BEXT")?std::atof(std::getenv("BEXT")):0.0; gExternalMagneticField=Vec3{sz*std::cos(ph),sz*std::sin(ph),cz}*B;
      std::fprintf(stderr,"BFIELD idx %d %.12e %.12e %.12e\n",i,gExternalMagneticField.x,gExternalMagneticField.y,gExternalMagneticField.z); }
    std::fprintf(stderr,"START idx %d\n",i); std::fflush(stderr);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(500ULL+i),1,1200.0);
    std::fprintf(stderr,"END idx %d stop %d photons %lld\n",i,(int)e.stopCause,(long long)e.emittedPhotonCount); std::fflush(stderr);
  }
}
