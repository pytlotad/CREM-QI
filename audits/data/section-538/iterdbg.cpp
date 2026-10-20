#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false; setenv("CREM_NO_SPIN_CONTACT_LEVEL","1",1);
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  const Vec3 z{0,0,1}, x{1,0,0}; gExternalMagneticField=z*1.0;
  const Vec3 n1{1,0,0}, n2{1,0,0};
  SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x; s.firstDipole=(n1*(1.0/n1.norm()))*(-firstMagneticMoment); s.secondDipole=(n2*(1.0/n2.norm()))*secondMagneticMoment;
  auto R=[&](const Vec3& a1,const Vec3& a2,const Vec3& LL){ return secularRatesWithSStateRule(a,LL,a1,a2,mu,0.0,x); };
  auto step=[&](Vec3& m1,Vec3& m2,const Vec3& LL,double dt){
    auto r0=R(m1,m2,LL); m1=rotateDipoleByAngularVelocity(m1,r0.firstExternal,dt/2); m2=rotateDipoleByAngularVelocity(m2,r0.secondExternal,dt/2);
    auto r1=R(m1,m2,LL); Vec3 h1=rotateDipoleByAngularVelocity(m1,r1.first-r1.firstExternal,dt/2), h2=rotateDipoleByAngularVelocity(m2,r1.second-r1.secondExternal,dt/2);
    auto r2=R(h1,h2,LL); m1=rotateDipoleByAngularVelocity(m1,r2.first-r2.firstExternal,dt); m2=rotateDipoleByAngularVelocity(m2,r2.second-r2.secondExternal,dt);
    m1=rotateDipoleByAngularVelocity(m1,r2.firstExternal,dt/2); m2=rotateDipoleByAngularVelocity(m2,r2.secondExternal,dt/2); };
  const double h=2*pi*hbar, Z=(std::abs(firstGFactor)/2)*9.2740100783e-24, K=204.6577e9*h/1.748264/4;
  auto U=[&](const Vec3& a1,const Vec3& a2){ const Vec3 q1=a1*(-1.0/a1.norm()), q2=a2*(1.0/a2.norm()); return dot(q1,q2)+(Z/K)*(q1.z-q2.z); };
  Vec3 m1=s.firstDipole, m2=s.secondDipole;
  for(int k=0;k<=300;++k){
    if(k%20==0) std::printf("call %4d: U engine %+.8f manual %+.8f diff %.3e Lhat %+.2e %+.2e %+.9f\n",k,U(s.firstDipole,s.secondDipole),U(m1,m2),(s.firstDipole-m1).norm()/m1.norm(),s.orbitalAngularMomentum.x/s.orbitalAngularMomentum.norm(),s.orbitalAngularMomentum.y/s.orbitalAngularMomentum.norm(),s.orbitalAngularMomentum.z/s.orbitalAngularMomentum.norm());
    auto r0=R(s.firstDipole,s.secondDipole,s.orbitalAngularMomentum); const double sp=std::max(r0.first.norm(),r0.second.norm());
    const double total=0.06/sp, dt1=0.05/sp;
    s=advanceCoupledSecularSpinOrbit(s,a,mu,total).state;
    auto q0=R(m1,m2,z*L); const double spm=std::max(q0.first.norm(),q0.second.norm());
    step(m1,m2,z*L,0.05/spm); step(m1,m2,z*L,0.06/spm-0.05/spm); }
}
