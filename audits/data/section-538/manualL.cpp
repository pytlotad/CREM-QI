#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false; setenv("CREM_NO_SPIN_CONTACT_LEVEL","1",1);
  const double aB=pairBohrRadius(activePair), a=aB, mu=firstMass*secondMass/(firstMass+secondMass), h=2*pi*hbar;
  const Vec3 z{0,0,1}, x{1,0,0}; gExternalMagneticField=z*1.0;
  const double Z=(std::abs(firstGFactor)/2)*9.2740100783e-24, K=204.6577e9*h/1.748264/4;
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  auto Sv=[&](const Vec3& a1,const Vec3& a2){ return a1/g1+a2/g2; };
  Vec3 m1=x*(-firstMagneticMoment), m2=x*secondMagneticMoment; Vec3 Lv=z*(0.5*hbar);
  const bool updateL=std::atoi(argv[1])!=0;
  auto R=[&](const Vec3& a1,const Vec3& a2,const Vec3& LL){ return secularRatesWithSStateRule(a,LL,a1,a2,mu,0.0,x); };
  auto U=[&](){ const Vec3 n1=m1*(-1.0/m1.norm()), n2=m2*(1.0/m2.norm()); return dot(n1,n2)+(Z/K)*(n1.z-n2.z); };
  for(int k=0;k<=6000;++k){
    if(k%1000==0) std::printf("updateL %d step %4d: U %+.6f Lhat.x %+.3e\n",(int)updateL,k,U(),Lv.x/Lv.norm());
    auto r0=R(m1,m2,Lv); const double sp=std::max(r0.first.norm(),r0.second.norm()); const double dt=0.05/sp;
    m1=rotateDipoleByAngularVelocity(m1,r0.firstExternal,dt/2); m2=rotateDipoleByAngularVelocity(m2,r0.secondExternal,dt/2);
    Vec3 J=Lv+Sv(m1,m2);
    auto r1=R(m1,m2,Lv); Vec3 h1=rotateDipoleByAngularVelocity(m1,r1.first-r1.firstExternal,dt/2), h2=rotateDipoleByAngularVelocity(m2,r1.second-r1.secondExternal,dt/2);
    const Vec3 Lmid=updateL?J-Sv(h1,h2):Lv;
    auto r2=R(h1,h2,Lmid); m1=rotateDipoleByAngularVelocity(m1,r2.first-r2.firstExternal,dt); m2=rotateDipoleByAngularVelocity(m2,r2.second-r2.secondExternal,dt);
    if(updateL) Lv=J-Sv(m1,m2);
    m1=rotateDipoleByAngularVelocity(m1,r2.firstExternal,dt/2); m2=rotateDipoleByAngularVelocity(m2,r2.secondExternal,dt/2); }
}
