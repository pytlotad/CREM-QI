// Audit 636: one o-Ps estimator cascade per argv index (seed splitMix64(623000 + i), as probe 623), printing what the validation check
// excited-start-consistency reads: stop cause, photon count, tau = 1/Gamma(stop) and the source-frame energy of the last photon.
// env NLEV (start level), LFRAC (L / n hbar).
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
static double envd(const char* k,double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
int main(int argc,char**argv){
  gSpinQuantization=true; gDeterministicEmission=true; gMeasureCollapseTransit=false;
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
  gInitialAngularMomentumFraction=envd("LFRAC",0.5);
  const double nl=envd("NLEV",1.0); gInitialSeparationOverride=nl*nl*pairBohrRadius(activePair);
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    const CremCollapseEstimate e=estimateCremCollapse(splitMix64(623000ULL+i),2,3000.0);
    double last=0.0;
    if(!e.labFramePhotons.empty()){ const LabFramePhoton& p=e.labFramePhotons.back(); last=p.energyJoules;
      if(p.sourceBeta>1.0e-9&&std::isfinite(p.angleFromRecoilAxisRadians))
        last*=(1.0-p.sourceBeta*std::cos(p.angleFromRecoilAxisRadians))/std::sqrt(1.0-p.sourceBeta*p.sourceBeta); }
    std::printf("RUN %d nlev %g stop %d photons %lld tau %.12e lastEV %.12f w %.6e n %.12e\n",i,nl,(int)e.stopCause,(long long)e.emittedPhotonCount,
                1.0/e.annihilationRateAtStop,last/1.602176634e-19,e.annihilationTwoPhotonWeightAtStop,e.annihilationContactDensityAtStop);
    std::fflush(stdout); }
}
