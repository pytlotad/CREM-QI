// Audit 641: estimator bookkeeping fields of one cascade (seed splitMix64(623000 + i), as probe 623): energy credited to the orbit from the spin coupling,
// pending, clamps, wave lock, hidden-J bookings, and the orbit energy against the Bohr ladder at the start.  env PAIR, PHEN, NLEV, LFRAC.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
static double envd(const char* k,double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gSpinQuantization=true; gDeterministicEmission=true; gMeasureCollapseTransit=false;
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
  gInitialAngularMomentumFraction=envd("LFRAC",0.75);
  const double nl=envd("NLEV",2.0); gInitialSeparationOverride=nl*nl*pairBohrRadius(activePair);
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(623000ULL+i),(int)envd("PHEN",1),3000.0);
    std::printf("RUN %d photons %lld spinCredited %.6e J (%.4f GHz) pending %.6e J clamps %llu waveLockWork %.3e hiddenBookings %llu\n",i,
      (long long)e.emittedPhotonCount,e.spinEnergyCreditedJoules,e.spinEnergyCreditedJoules/(2*pi*hbar)/1e9,e.spinEnergyPendingJoules,
      (unsigned long long)e.spinEnergyClamps,e.waveLockWorkJoules,(unsigned long long)e.hiddenAngularBookings); std::fflush(stdout); }
}
