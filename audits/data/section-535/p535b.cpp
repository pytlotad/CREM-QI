#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 535b: o-Ps 1S (triplet, parallel spins at angle theta to B), B = 0.1 T: <w>(theta) over 1 ns; argv[1] = theta in degrees.
int main(int argc,char**argv){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false;
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  const Vec3 z{0,0,1}, x{1,0,0}; gExternalMagneticField=z*0.1;
  const double th=std::atof(argv[1])*pi/180; const Vec3 sh{std::sin(th),0,std::cos(th)};
  SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x; s.firstDipole=sh*(-firstMagneticMoment); s.secondDipole=sh*secondMagneticMoment;
  const double dt=2.5e-13, T=1e-9; double wsum=0; long n=0;
  for(double t=0;t<T-1e-18;t+=dt){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state; wsum+=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); ++n; }
  std::printf("theta %5.1f deg: <w> %.6e\n",std::atof(argv[1]),wsum/n);
}
