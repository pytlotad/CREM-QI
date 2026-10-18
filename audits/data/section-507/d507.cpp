#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 507: radiation driven by the dipole force alone -- hypothetical lifetime if it triggered the photon.  argv: pair.
int main(int argc,char**argv){
  applyPairFromOption(argv[1]);
  const double aB=pairBohrRadius(activePair), M=firstMass+secondMass, mu=firstMass*secondMass/M, a=4*aB, Lm=1.5*hbar, k=pairCoulombStrength;
  const double nmot=std::sqrt(k/(mu*a*a*a)), ecc=std::sqrt(1-std::pow(Lm/(2*hbar),2)), ec=std::sqrt(1-ecc*ecc);
  const double fe=std::pow(1-ecc*ecc,1.5)/(1+0.5*ecc*ecc)*coRotatingOverNetCount(ecc), Eg=k/(2*aB)*0.75;
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  const Vec3 P{1,0,0}, Q{0,1,0}, Z{0,0,1};
  struct Cfg{const char* n; Vec3 s1,s2;};
  const Cfg cf[]={{"both along L",Z,Z},{"both against L",Z*-1,Z*-1},{"both perpendicular",P,P},{"antiparallel, 1 along L",Z,Z*-1},{"antiparallel, perpendicular",P,P*-1}};
  const double kL=1.0/(6*pi*8.8541878128e-12*c*c*c);
  std::printf("%s: n 2, L 1.5 hbar, e %.4f, f(e) %.4f, E_gamma %.4e J, hbar omega %.4e J, larmorOrbitAveragedPower %.5e W\n",argv[1],ecc,fe,Eg,hbar*nmot,larmorOrbitAveragedPower(a,ecc));
  for(const Cfg& q: cf){
    double sC=0,sD=0,sX=0,sw=0; const int N=2048;
    for(int j=0;j<N;++j){ const double E=2*pi*(j+0.5)/N, w=1-ecc*std::cos(E);
      const Vec3 rr=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E)), vv=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w);
      State s{}; s.firstPosition=rr*(secondMass/M); s.secondPosition=rr*(-firstMass/M); s.firstVelocity=vv*(secondMass/M); s.secondVelocity=vv*(-firstMass/M);
      s.firstProperDipole=q.s1*(g1*0.5*hbar); s.secondProperDipole=q.s2*(g2*0.5*hbar); synchronizeCovariantDipoles(s);
      const StateHistory h{State{s}};
      gDipoleForceEnabled=true; const MutualForces fOn=retardedExternalForces(s,h);
      gDipoleForceEnabled=false; const MutualForces fOff=retardedExternalForces(s,h); gDipoleForceEnabled=true;
      const Vec3 dC=fOff.first*(firstCharge/firstMass)+fOff.second*(secondCharge/secondMass);
      const Vec3 dD=(fOn.first-fOff.first)*(firstCharge/firstMass)+(fOn.second-fOff.second)*(secondCharge/secondMass);
      sC+=dC.squaredNorm()*w; sD+=dD.squaredNorm()*w; sX+=2*dot(dC,dD)*w; sw+=w; }
    const double PC=kL*sC/sw, PD=kL*sD/sw, PX=kL*sX/sw;
    std::printf("  [%s] P_E1 %.5e W, P_dip %.4e W (ratio %.3e), P_cross/P_E1 %+.3e; hypothetical tau_E = E_gamma/P_dip %.3e s, tau_L = hbar omega/(f P_dip) %.3e s\n",
      q.n,PC,PD,PD/PC,PX/PC,Eg/PD,hbar*nmot/(fe*PD));
  }
}
