#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 426: HFS of 1s and 2s (L = hbar/2) from spinCouplingEnergy, prepared-spin pair (gSingletChannelEligible).
int main(){
  gSingletChannelEligible=true;
  const double a1=pairBohrRadius(activePair), mu=firstMagneticMoment, h=2*pi*hbar, L=0.5*hbar;
  const Vec3 ax{0,0,1}, P{1,0,0}, z{0,0,1};
  for(int n: {1,2}){ const double a=a1*n*n, e=std::sqrt(1.0-0.25/(n*n));
    auto U=[&](const Vec3& m1,const Vec3& m2){ return spinCouplingEnergy(a,e,ax,P,m1,m2,L); };
    // e-: moment opposite to spin; singlet = antiparallel spins = parallel moments.
    const double singlet=U(z*mu,z*mu), triplet=U(z*mu,z*(-mu));
    std::printf("n=%d: U_singlet=%.6f GHz U_triplet=%.6f GHz HFS=%.6f GHz\n",n,singlet/h*1e-9,triplet/h*1e-9,(triplet-singlet)/h*1e-9); }
}
