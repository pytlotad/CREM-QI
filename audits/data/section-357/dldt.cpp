#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
// Audit 357: d|L|/dt = -Lhat . dS/dt from the secular rates (L = J - S), o-Ps and p-Ps at n = 1, L = hbar/2,
// raw vs s-state-isotropic rates; decomposition: orbital part only, partner part.
int main(){
  const double a1=pairBohrRadius(activePair), mu=firstMagneticMoment, mred=firstMass*secondMass/(firstMass+secondMass);
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  std::mt19937_64 rng(357); std::normal_distribution<double> N(0,1);
  const Vec3 Lh{0,0,1}, P{1,0,0}; const Vec3 Lv=Lh*(0.5*hbar);
  for(int ch=1;ch<=2;++ch){
    double mxRaw=0,mxIso=0,mxOrb=0,mxW=0;
    for(int k=0;k<50;++k){
      Vec3 s{N(rng),N(rng),N(rng)}; s=s/s.norm();
      const Vec3 m1=s*mu, m2=s*(ch==1?mu:-mu);
      auto rate=[&](const OrbitAveragedBmtAngularVelocities& r){ const Vec3 dS=cross(r.first,m1)/g1+cross(r.second,m2)/g2; return -dot(Lh,dS)/hbar*1e-9; }; // hbar per ns
      auto raw=orbitAveragedBmtAngularVelocities(a1,Lv,m1,m2,mred,0.0,P);
      auto iso=raw; applySStateIsotropy(iso,a1,Lv,m1,m2,mred,0.0,P);
      auto orb=orbitAveragedBmtAngularVelocities(a1,Lv,m1*1e-12,m2*1e-12,mred,0.0,P);
      mxRaw=std::max(mxRaw,std::abs(rate(raw))); mxIso=std::max(mxIso,std::abs(rate(iso))); mxOrb=std::max(mxOrb,std::abs(rate(orb)));
      // does omega depend on the OWN moment?  omega_1 with partner tiny, own moment full vs own tiny
      auto own=orbitAveragedBmtAngularVelocities(a1,Lv,m1,m2*1e-12,mred,0.0,P);
      mxW=std::max(mxW,(own.first-orb.first).norm());
    }
    std::printf("%s: max |d|L|/dt| [hbar/ns]: raw %.3e  isotropic %.3e  orbital-only %.3e ; own-moment dependence of omega1 %.3e rad/s\n",
      ch==1?"p-Ps":"o-Ps",mxRaw,mxIso,mxOrb,mxW);
  }
}
