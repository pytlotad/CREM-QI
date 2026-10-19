#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 537 (side question): o-Ps 1S "m = 0" in a field -- does the orbit absorb the change of S (L = J - S) and leave |L| = hbar/2?
int main(int argc,char**argv){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false;
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  const Vec3 z{0,0,1}, x{1,0,0}; const double B=std::atof(argv[1]); gExternalMagneticField=z*B;
  SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x; s.firstDipole=x*(-firstMagneticMoment); s.secondDipole=x*secondMagneticMoment;
  const double dt=2.5e-13, T=1e-9; double wsum=0, Lmin=9, Lmax=0, Smin=9, Jmin=9, Jmax=0; long n=0;
  for(double t=0;t<T-1e-18;t+=dt){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state; wsum+=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); ++n;
    const double l=s.orbitalAngularMomentum.norm()/hbar; Lmin=std::min(Lmin,l); Lmax=std::max(Lmax,l);
    const Vec3 S=(s.firstDipole*(-1.0/s.firstDipole.norm())+s.secondDipole*(1.0/s.secondDipole.norm()))*0.5; Smin=std::min(Smin,S.norm());
    const double J=(s.orbitalAngularMomentum+S*hbar).norm()/hbar; Jmin=std::min(Jmin,J); Jmax=std::max(Jmax,J); }
  std::printf("B %.1f T (CREM_NO_SPIN_CONTACT_LEVEL %s): <w> %.4e; |L|/hbar in [%.4f, %.4f]; min |S|/hbar %.4f; |J|/hbar in [%.4f, %.4f]\n",B,std::getenv("CREM_NO_SPIN_CONTACT_LEVEL")?"set":"-",wsum/n,Lmin,Lmax,Smin,Jmin,Jmax);
}
