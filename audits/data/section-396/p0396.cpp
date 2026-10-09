#include <cmath>
#include <cstdio>
#include <random>
#include "modules/multipole_photons.hpp"
// Audit 396: checks of the multipole module -- E2 sum rule against Peters' f(e), circle, angular sampler, Delta l table.
int main(){
  double worst=0;
  for(double e: {0.0,0.1,0.3,0.5,0.6,0.7,0.8,0.9}){ const auto& s=keplerQuadrupoleSpectrum(e);
    const double ee=std::sqrt(std::round(e*e*4096)/4096.0); const double ratio=s.totalPower/(288.0*petersEccentricityFactor(ee));
    worst=std::max(worst,std::abs(ratio-1));
    int m2=0,m0=0,mm2=0; double p2=0,p0=0,pm=0; for(const auto& h: s.entries){ if(h.m==2){++m2;p2+=h.power;} else if(h.m==0){++m0;p0+=h.power;} else {++mm2;pm+=h.power;} }
    std::printf("e=%.3f sum/288f = %.10f  countFactor %.5f  shares m=+2 %.4f m=0 %.4f m=-2 %.4f  entries %zu\n",ee,ratio,s.countFactor,p2/s.totalPower,p0/s.totalPower,pm/s.totalPower,s.entries.size()); }
  std::printf("H396a worst |sum/288f - 1| (e<=0.9) = %.3e\n",worst);
  const auto& c0=keplerQuadrupoleSpectrum(0.0); std::printf("circle entries:"); for(const auto& h: c0.entries) std::printf(" (k=%d,m=%d,P=%.9f)",h.k,h.m,h.power); std::printf(" [%zu entries]",c0.entries.size()); std::printf("\n");
  std::mt19937_64 r(396); std::uniform_real_distribution<double> U(0,1);
  for(int m: {2,0}){ const int N=1000000; double s=0,s2=0; for(int i=0;i<N;++i){ const double c=sampleQuadrupoleCosTheta(m,U(r)); s+=c*c; s2+=c*c*c*c; }
    const double mean=s/N, sd=std::sqrt((s2/N-mean*mean)/N); const double exact=m==2?5.0/21.0:3.0/7.0;
    std::printf("H396b |m|=%d <cos^2> = %.6f +/- %.6f (exact %.6f, %.2f sd)\n",m,mean,sd,exact,(mean-exact)/sd); }
  std::printf("Delta l table (E2): l, m -> transfer hbar\n");
  for(int l=0;l<=3;++l){ std::printf("  l=%d:",l); for(int m: {2,0,-2}) std::printf(" m=%+d -> %+.0f",m,photonOrbitalTransferHbar(PhotonMultipole::E2,m,l+0.5,1.0)); std::printf("\n"); }
  std::printf("M1 transfer at l=1: %+.0f; E1 passes its sign: %+.0f\n",photonOrbitalTransferHbar(PhotonMultipole::M1,0,1.5,1.0),photonOrbitalTransferHbar(PhotonMultipole::E1,0,1.5,-1.0));
}
