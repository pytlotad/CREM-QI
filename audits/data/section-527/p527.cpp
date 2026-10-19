#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 527: fine-structure beat in the secular transport -- hydrogen 2p, proton moment x 1e-12, rotation of s^ about the conserved J_e = L + s.
int main(){
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), a=4*aB, L=1.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  for(double c:{-0.9,-0.75,-0.5,0.0,0.5,0.9}){
    const double sn=std::sqrt(1-c*c); const Vec3 Lh{0,0,1}, sh{sn,0,c};
    SecularSpinOrbitState s; s.orbitalAngularMomentum=Lh*L; s.secondDipole=sh*(-secondMagneticMoment); s.firstDipole=Vec3{0.3,0.5,0.81}*(firstMagneticMoment*1e-12);
    s.periapsisDirection=Vec3{1,0,0};
    const double dt=3e-12, T=2e-9; double ang=0; Vec3 prev; bool have=false; double Jmin=1e9,Jmax=0,cmn=2,cmx=-2;
    for(double t=0;t<T-1e-15;t+=dt){
      s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state;
      const Vec3 se=s.secondDipole*(-1.0/s.secondDipole.norm()), Je=s.orbitalAngularMomentum+se*(0.5*hbar), Jh=Je*(1.0/Je.norm());
      Jmin=std::min(Jmin,Je.norm()/hbar); Jmax=std::max(Jmax,Je.norm()/hbar); const double cc=dot(se,s.orbitalAngularMomentum)/s.orbitalAngularMomentum.norm(); cmn=std::min(cmn,cc); cmx=std::max(cmx,cc);
      Vec3 p=se-Jh*dot(se,Jh); if(p.norm()>1e-9){ p=p*(1.0/p.norm()); if(have) ang+=std::atan2(dot(cross(prev,p),Jh),dot(prev,p)); prev=p; have=true; } }
    const double Jn=std::sqrt(2.25+0.25+1.5*c), f=std::abs(ang)/(2*pi*T)*1e-6, pred=10968.566*Jn/1.5;
    std::printf("cos(L,s) %+.2f (basin j = %s): |J_e| %.4f [%.4f, %.4f] hbar, cos in [%+.5f, %+.5f]; beat %.1f MHz; xi|J_e| prediction %.1f; /interval 10968.70: %.4f; measured 10969.045\n",
      c,c<-1.0/3?"1/2":"3/2",Jn,Jmin,Jmax,cmn,cmx,f,pred,f/10968.70);
  }
}
