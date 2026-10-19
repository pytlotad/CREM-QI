#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 529: Ps 2^3P beat in the secular transport -- rotation of S^ (parallel triplet spins) about the conserved J_e = L + S.
int main(){
  applyPairFromOption("electron,positron");
  gSingletChannelEligible=std::getenv("PS_CHANNEL")!=nullptr;   // audit 530: the triplet channel of the level energy
  const double aB=pairBohrRadius(activePair), a=4*aB, L=1.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  std::printf("Ps 2p, channel %d, CREM_RULE_DIPOLE_TENSOR %s, CREM_NO_PS_SPIN_RULE %s\n",(int)gSingletChannelEligible,std::getenv("CREM_RULE_DIPOLE_TENSOR")?"on":"off",std::getenv("CREM_NO_PS_SPIN_RULE")?"set":"-");
  for(double Jt:{0.6,1.0,1.5,2.0,2.4}){
    const double c=(Jt*Jt-3.25)/3.0, sn=std::sqrt(1-c*c); const Vec3 Sh{sn,0,c};
    SecularSpinOrbitState s; s.orbitalAngularMomentum=Vec3{0,0,1}*L;
    s.firstDipole=Sh*(-firstMagneticMoment); s.secondDipole=Sh*(secondMagneticMoment);   // s = sign(q) m^: electron -, positron +
    s.periapsisDirection=Vec3{1,0,0};
    const double dt=3e-12, T=2e-9; double ang=0, pmin=2, Jmin=1e9, Jmax=0; Vec3 prev; bool have=false;
    for(double t=0;t<T-1e-15;t+=dt){
      s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state;
      const Vec3 s1=s.firstDipole*(-1.0/s.firstDipole.norm()), s2=s.secondDipole*(1.0/s.secondDipole.norm());
      pmin=std::min(pmin,dot(s1,s2)); const Vec3 S=(s1+s2)*(0.5*hbar), Je=s.orbitalAngularMomentum+S, Jh=Je*(1.0/Je.norm());
      Jmin=std::min(Jmin,Je.norm()/hbar); Jmax=std::max(Jmax,Je.norm()/hbar);
      Vec3 p=S-Jh*dot(S,Jh); if(p.norm()>1e-12*hbar){ p=p*(1.0/p.norm()); if(have) ang+=std::atan2(dot(cross(prev,p),Jh),dot(prev,p)); prev=p; have=true; } }
    const double f=std::abs(ang)/(2*pi*T)*1e-9;
    std::printf("|J_e| %.2f (cos(L,S) %+.4f): beat %.4f GHz; |J_e| in [%.4f, %.4f]; min s1.s2 %.6f%s\n",Jt,c,f,Jmin,Jmax,pmin,
      Jt==1.0?"   <- midpoint of P0/P1: measured P1-P0 5.487 GHz":(Jt==2.0?"   <- midpoint of P1/P2: measured P2-P1 4.388 GHz":""));
  }
}
