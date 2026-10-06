#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
static const char* name(CollapseStopCause c){ switch(c){case CollapseStopCause::ComptonBarrier:return "barrier";
 case CollapseStopCause::RetardationLimit:return "retard"; case CollapseStopCause::GroundStateFloor:return "floor";
 case CollapseStopCause::EmissionChannelClosed:return "closed"; default:return "none";} }
int main(int argc,char**argv){
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
  gInitialPrincipalLevel=std::atoi(std::getenv("LEVEL")); gMicrocanonicalStart=false; gSpinQuantization=true;
  const double mu=firstMass*secondMass/(firstMass+secondMass), K=pairCoulombStrength/mu, a1=pairBohrRadius(activePair);
  const double psi0=1.0/(pi*a1*a1*a1), eps=magneticDipoleRadius();
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    for(int ph=1;ph<=1;++ph){ std::printf("BEGIN %d %d\n",i,ph); std::fflush(stdout);
      const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+i),ph,1200.0);
      const double A=e.terminalSemiMajorAxis, L=e.terminalAngularMomentum*hbar/mu;   // specific
      const double e2=std::max(0.0,1.0-L*L/(K*A)), ecc=std::sqrt(e2), rp=A*(1-ecc);
      const double n=orbitAveragedContactDensity(A,std::min(ecc,0.999999));
      std::printf("END idx %2d %s stop=%s photons=%lld nE=%.4f L=%.4f e=%.5f rp/eps=%.4e n/psi0=%.4e ET=%.6e S=%.6f Gf=%.6e check=%.6e\n",
        i,ph==1?"para ":"ortho",name(e.stopCause),(long long)e.emittedPhotonCount,std::sqrt(a1/A),e.terminalAngularMomentum,ecc,rp/eps,n/psi0,
        e.annihilationMeanLifetimeSeconds,e.annihilationSurvivalAtStop,e.annihilationRateAtStop,
        e.annihilationRateAtStop>0?e.annihilationSurvivalAtStop/e.annihilationRateAtStop:0.0);
      std::printf("CAP events=%lld largest=%.3e hbar\n",(long long)e.orbitalCapEvents,e.orbitalCapLargestExcessHbar);
      std::fflush(stdout);
    } }
}
