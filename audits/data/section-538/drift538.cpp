#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 538: invariants of the two-spin Ps 1S system in a field along the engine's transport: U/K = n1.n2 + (Z/K)(n1z - n2z) and n1z + n2z, every 0.1 ns.
int main(int argc,char**argv){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false;
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), h=2*pi*hbar;
  const double B=std::atof(argv[1]), muB=9.2740100783e-24, g=std::abs(firstGFactor), Z=(g/2)*muB*B, K=204.6577e9*h/1.748264/4*sStateContactTorqueFactor();
  const Vec3 z{0,0,1}, x{1,0,0}; gExternalMagneticField=z*B;
  SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x; s.firstDipole=x*(-firstMagneticMoment); s.secondDipole=x*secondMagneticMoment;
  for(int i=0;i<=4000;++i){
    const Vec3 n1=s.firstDipole*(-1.0/s.firstDipole.norm()), n2=s.secondDipole*(1.0/s.secondDipole.norm());
    if(i%400==0) std::printf("t %.2f ns: U/K %+.6f  n1z+n2z %+.6f  w %.4e\n",i*2.5e-4,dot(n1,n2)+(Z/K)*(n1.z-n2.z),n1.z+n2.z,contactTwoPhotonWeight(s.firstDipole,s.secondDipole));
    s=advanceCoupledSecularSpinOrbit(s,a,mu,2.5e-13).state; }
}
