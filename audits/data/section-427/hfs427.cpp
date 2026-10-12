#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 427: HFS of 1s, 2s (L = hbar/2) from spinCouplingEnergy with the virtual-annihilation term; env PAIR.
int main(){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gSingletChannelEligible=true;
  const double a1=pairBohrRadius(activePair), mu1=firstMagneticMoment, mu2=secondMagneticMoment, h=2*pi*hbar, L=0.5*hbar;
  const Vec3 ax{0,0,1}, P{1,0,0}, z{0,0,1};
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1;
  std::printf("pair: m1=%.6e m2=%.6e antiparticle=%d\n",firstMass,secondMass,(int)particleAntiparticlePair());
  for(int n: {1,2}){ const double a=a1*n*n, e=std::sqrt(1.0-0.25/(n*n));
    auto U=[&](const Vec3& m1,const Vec3& m2){ return spinCouplingEnergy(a,e,ax,P,m1,m2,L); };
    // spins along +z / -z; moment = sign(charge) mu along spin
    const double singlet=U(z*(s1*mu1),z*(-s2*mu2)), triplet=U(z*(s1*mu1),z*(s2*mu2));
    std::printf("n=%d: U_singlet=%.6f GHz U_triplet=%.6f GHz HFS=%.6f GHz\n",n,singlet/h*1e-9,triplet/h*1e-9,(triplet-singlet)/h*1e-9); }
  const double nqr=quiggaRosnerContactDensity(a1,L), Gp=4*pi*classicalElectronRadius*classicalElectronRadius*c*nqr;
  std::printf("n=1: Gamma_p=%.6e /s (tau %.4f ps); hbar Gamma_p/(2 alpha) = %.6f GHz; alpha^4 m c^2/(4h) = %.6f GHz\n",
    Gp,1e12/Gp,hbar*Gp/(2*fineStructureConstant)/h*1e-9,std::pow(fineStructureConstant,4)*electronMass*c*c/4/h*1e-9);
}
