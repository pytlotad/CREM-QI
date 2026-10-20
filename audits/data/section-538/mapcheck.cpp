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
  double worst=0; int worstk=-1;
  for(int k=0;k<300;++k){
    auto r0=R(s.firstDipole,s.secondDipole,s.orbitalAngularMomentum); const double sp=std::max(r0.first.norm(),r0.second.norm());
    const double total=0.06/sp, dt1=0.05/sp;
    Vec3 m1=s.firstDipole, m2=s.secondDipole; const Vec3 LL=s.orbitalAngularMomentum;
    step(m1,m2,LL,dt1); step(m1,m2,LL,total-dt1);
    const auto ad=advanceCoupledSecularSpinOrbit(s,a,mu,total);
    const double d=(ad.state.firstDipole-m1).norm()/m1.norm();
    if(d>worst){ worst=d; worstk=k; }
    if(d>1e-8) std::printf("call %d: map discrepancy %.3e (substeps %d), cos12 %.6f\n",k,d,ad.substeps,dot(s.firstDipole,s.secondDipole)/(s.firstDipole.norm()*s.secondDipole.norm()));
    s=ad.state; }
  std::printf("max per-call discrepancy %.3e at call %d\n",worst,worstk);
}
