// Audit 624: mechanical photons on the production rules.  Trajectory mode (simulate) from a Langer state (NLEV, LHBAR = l)
// in the photon mode; the first photons are read from the engine's CREM_DEBUG_LADDER / RULES lines on stderr.  At the start the
// probe prints the rule shares of that (n, l) from the same functions (EXPECT lines), for the comparison in check624.py.
// env PAIR, NLEV, LHBAR, OBS [orbits], CREM_MECHANICAL_HAZARD_SCALE; argv: run indices.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
static double envd(const char* k,double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static void expect(const char* tag,const std::vector<double>& w,int n,int firstFinal,int step){
  double t=0.0; for(double x: w) t+=x;
  std::printf("EXPECT %s total %.12e",tag,t);
  for(size_t i=0;i<w.size();++i) std::printf(" n'=%d:%.9e",firstFinal+step*static_cast<int>(i),t>0.0?w[i]/t:0.0);
  std::printf("\n"); (void)n; }
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
  gSpinQuantization=true;
  const int n=(int)envd("NLEV",4), l=(int)envd("LHBAR",1);
  gInitialSeparationOverride=n*n*pairBohrRadius(activePair);
  gInitialAngularMomentumFraction=(l+0.5)/n;
  const double level=n, L=l+0.5;
  // rule shares of the start state (the same functions the engine calls)
  if(l==0){ expect("sE1",sStateChannelRates(level),n,n-1,-1);
            const auto q=sStateQuadrupoleChannelRates(level); expect("sE2",q,n,n-1,-1); }
  else { const auto d=downChannelShares(n,l); std::vector<double> rev(d.rbegin(),d.rend());
         expect("down",rev,n,n-1,-1);   // index 0 -> n' = n - 1
         expect("up",upChannelRates(level,L),n,n-1,-1); }
  const double period=2.0*pi*std::sqrt(std::pow(n*n*pairBohrRadius(activePair),3)*pairReducedMass/pairCoulombStrength);
  std::printf("PERIOD %.9e\n",period); std::fflush(stdout);
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    SimulationOptions o; o.frameCount=2; o.observationTime=envd("OBS",0.3)*period;
    std::fprintf(stderr,"START idx %d\n",i); std::fflush(stderr);
    const SimulationResult r=simulate(splitMix64(624000ULL+i),1,o);
    std::fprintf(stderr,"END idx %d outcome %d elapsed %.6e\n",i,(int)r.outcome,r.elapsedTime); std::fflush(stderr);
  }
}
