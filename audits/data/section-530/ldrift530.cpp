#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true;
  const double aB=pairBohrRadius(activePair), a=4*aB, L=1.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  for(double c:{-0.34,0.86}){ const double sn=std::sqrt(1-c*c); const Vec3 Sh{sn,0,c};
    SecularSpinOrbitState s; s.orbitalAngularMomentum=Vec3{0,0,1}*L; s.firstDipole=Sh*(-firstMagneticMoment); s.secondDipole=Sh*secondMagneticMoment; s.periapsisDirection=Vec3{1,0,0};
    double t=0; const double dt=1e-10; double Lmin=9,Lmax=0;
    for(;t<3.35e-9;t+=dt){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state; const double x=s.orbitalAngularMomentum.norm()/hbar; Lmin=std::min(Lmin,x); Lmax=std::max(Lmax,x); }
    std::printf("cos %+.2f: |L|/hbar after 3.35 ns %.6f (min %.6f max %.6f)\n",c,s.orbitalAngularMomentum.norm()/hbar,Lmin,Lmax); }
}
