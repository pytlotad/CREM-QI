#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 534: p-Ps 2p singlet in the secular transport, 2 ns, B = 0: w = (|mu1 + mu2|/2 mu)^2, s1.s2, rotation of the pair axis.
int main(){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true;
  gSingletSpinTransport=std::getenv("CREM_NO_SINGLET_TRANSPORT")==nullptr;
  const double aB=pairBohrRadius(activePair), a=4*aB, L=1.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  std::printf("singlet transport %d, CREM_NO_PS_SINGLET_RULE %s\n",(int)gSingletSpinTransport,std::getenv("CREM_NO_PS_SINGLET_RULE")?"set":"-");
  for(double c:{0.3,-0.7}){
    const Vec3 sh{std::sqrt(1-c*c),0,c};
    SecularSpinOrbitState s; s.orbitalAngularMomentum=Vec3{0,0,1}*L; s.firstDipole=sh*(-firstMagneticMoment); s.secondDipole=(sh*(-1.0))*secondMagneticMoment; s.periapsisDirection=Vec3{1,0,0};
    double wmin=2, ssmax=-2, Lmin=9, Lmax=0; const double dt=3e-12, T=2e-9;
    for(double t=0;t<T-1e-15;t+=dt){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state;
      const double w=std::pow((s.firstDipole+s.secondDipole).norm()/(2*firstMagneticMoment),2); wmin=std::min(wmin,w);
      const Vec3 s1=s.firstDipole*(-1.0/s.firstDipole.norm()), s2=s.secondDipole*(1.0/s.secondDipole.norm()); ssmax=std::max(ssmax,dot(s1,s2));
      Lmin=std::min(Lmin,s.orbitalAngularMomentum.norm()/hbar); Lmax=std::max(Lmax,s.orbitalAngularMomentum.norm()/hbar); }
    const Vec3 s1=s.firstDipole*(-1.0/s.firstDipole.norm());
    std::printf("cos(L,s1) %+.1f: min w %.12f, max s1.s2 %+.12f, |L| in [%.6f, %.6f], pair axis turned %.3f deg in 2 ns\n",c,wmin,ssmax,Lmin,Lmax,std::acos(std::clamp(dot(s1,sh),-1.0,1.0))*180/pi);
  }
}
