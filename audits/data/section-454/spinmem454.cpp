// Audit 454: can spin motion be the memory / reservoir of the Langer attractor?  Spin precession per orbit vs apsidal precession per
// orbit (pi alpha^2 / L^2 for the Darwin Hamiltonian, audit 449) on Langer orbits of Ps (engine BMT rates).
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("electron,positron");
  const double mu=firstMass*secondMass/(firstMass+secondMass), aB=pairBohrRadius(activePair), k=pairCoulombStrength, alpha=k/(hbar*c);
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1; const Vec3 z{0,0,1}, x{1,0,0};
  std::printf("# n l  L/hbar  spin-prec/orbit [rad] (in-plane spins)  apsidal/orbit pi a^2/L^2  ratio\n");
  for(int n=1;n<=5;++n) for(int l=0;l<n;++l){
    const double a=aB*n*n, L=(l+0.5)*hbar, T=2*pi*std::sqrt(mu*a*a*a/k);
    const auto r=orbitAveragedBmtAngularVelocities(a,z*L,x*(s1*firstMagneticMoment),x*(s2*secondMagneticMoment),mu,0.0,x);
    const double sp=std::max(r.first.norm(),r.second.norm())*T, ap=pi*alpha*alpha/((l+0.5)*(l+0.5));
    std::printf("%d %d %.1f  %.6e  %.6e  %.5f\n",n,l,l+0.5,sp,ap,sp/ap); }
}
