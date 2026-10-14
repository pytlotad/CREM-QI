// Audit 448: spin (moment) precession rates vs orbital frequency and vs the de Broglie clock rate J omega.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  for(const char* pair: {"electron,positron","proton,electron"}){
    applyPairFromOption(pair);
    const double mu=firstMass*secondMass/(firstMass+secondMass), aB=pairBohrRadius(activePair), k=pairCoulombStrength;
    const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1;
    std::printf("pair %s\n",pair);
    for(auto [n,l]: {std::pair{1,0},{2,0},{2,1},{3,2}}){
      const double a=aB*n*n, L=(l+0.5)*hbar, T=2*pi*std::sqrt(mu*a*a*a/k), w=2*pi/T;
      const Vec3 z{0,0,1}, x{1,0,0};
      // spins along x (in plane) -- both moments; triplet-like parallel spins
      const auto r=orbitAveragedBmtAngularVelocities(a,z*L,x*(s1*firstMagneticMoment),x*(s2*secondMagneticMoment),mu,0.0,x);
      const double w1=r.first.norm(), w2=r.second.norm(), hfs=2*pi*1e9*0;   // HFS separately below
      std::printf("  n=%d l=%d: omega_orb %.4e rad/s; |Omega_1| %.4e (x %.3e), |Omega_2| %.4e (x %.3e); clock J omega = %.4e (spin/clock %.3e)\n",
        n,l,w,w1,w1/w,w2,w2/w,n*w,std::max(w1,w2)/(n*w));
    }
  }
  // mutual (HFS) precession of Ps 1s from the engine level splitting
  applyPairFromOption("electron,positron");
  const double aB=pairBohrRadius(activePair), hP=2*pi*hbar; gSingletChannelEligible=true;
  const double e1=std::sqrt(1-0.25), mu=firstMagneticMoment; const Vec3 z{0,0,1}, x{1,0,0};
  const double hfs=spinCouplingEnergy(aB,e1,z,x,z*(-mu),z*mu,0.5*hbar)-spinCouplingEnergy(aB,e1,z,x,z*(-mu),z*(-mu),0.5*hbar);
  const double mred=firstMass/2, T=2*pi*std::sqrt(mred*aB*aB*aB/pairCoulombStrength);
  std::printf("Ps 1s HFS %.4f GHz -> 2 pi HFS/h = %.4e rad/s, omega_orb = %.4e (ratio %.3e)\n",hfs/hP*1e-9,2*pi*hfs/hP,2*pi/T,(2*pi*hfs/hP)/(2*pi/T));
}
