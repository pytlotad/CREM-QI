#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false; setenv("CREM_NO_SPIN_CONTACT_LEVEL","1",1);
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), h=2*pi*hbar;
  const Vec3 z{0,0,1}, x{1,0,0}; gExternalMagneticField=z*1.0;
  const double Z=(std::abs(firstGFactor)/2)*9.2740100783e-24, K=204.6577e9*h/1.748264/4;
  Vec3 m1=x*(-firstMagneticMoment), m2=x*secondMagneticMoment; const Vec3 Lv=z*L;
  auto R=[&](const Vec3& a1,const Vec3& a2){ return secularRatesWithSStateRule(a,Lv,a1,a2,mu,0.0,x); };
  auto U=[&](){ const Vec3 n1=m1*(-1.0/m1.norm()), n2=m2*(1.0/m2.norm()); return dot(n1,n2)+(Z/K)*(n1.z-n2.z); };
  const double big=std::atof(argv[1]), small=std::atof(argv[2]);
  for(int k=0;k<=40000;++k){
    if(k%8000==0) std::printf("pattern %.3f/%.3f step %5d: U %+.6f\n",big,small,k,U());
    for(double ang:{big,small}){ if(ang<=0) continue;
      auto r0=R(m1,m2); const double sp=std::max(r0.first.norm(),r0.second.norm()); const double dt=ang/sp;
      m1=rotateDipoleByAngularVelocity(m1,r0.firstExternal,dt/2); m2=rotateDipoleByAngularVelocity(m2,r0.secondExternal,dt/2);
      auto r1=R(m1,m2); Vec3 h1=rotateDipoleByAngularVelocity(m1,r1.first-r1.firstExternal,dt/2), h2=rotateDipoleByAngularVelocity(m2,r1.second-r1.secondExternal,dt/2);
      auto r2=R(h1,h2); m1=rotateDipoleByAngularVelocity(m1,r2.first-r2.firstExternal,dt); m2=rotateDipoleByAngularVelocity(m2,r2.second-r2.secondExternal,dt);
      m1=rotateDipoleByAngularVelocity(m1,r2.firstExternal,dt/2); m2=rotateDipoleByAngularVelocity(m2,r2.secondExternal,dt/2); } }
}
