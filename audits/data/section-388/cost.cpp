#include "modules/crem_collapse.hpp"
#include <chrono>
#include <cstdio>
// Audit 388: wall time of the per-substep pieces at Ps 3s (a = 9 a_pair, L = hbar/2).
template<class F> static double timeIt(F f,int n){ const auto t0=std::chrono::steady_clock::now(); for(int i=0;i<n;++i) f(); return std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count()/n; }
int main(){
  const double mred=firstMass*secondMass/(firstMass+secondMass), a=9.0*pairBohrRadius(activePair);
  const Vec3 L{0,0,0.5*hbar}, m1=Vec3{0.3,0.4,0.866}*firstMagneticMoment, m2=Vec3{-0.5,0.1,0.86}*secondMagneticMoment, p{1,0,0};
  volatile double sink=0;
  const double tB=timeIt([&]{ sink+=orbitAveragedBmtAngularVelocities(a,L,m1,m2,mred,0.0,p).first.x; },5);
  const double tT=timeIt([&]{ sink+=orbitAveragedBmtAngularVelocities(a,L,m1*1e-12,m2*1e-12,mred,0.0,p).first.x; },5);
  const double tA=timeIt([&]{ sink+=orbitAveragedApsidalRate(a,L,m1,m2,mred,p).magnetic; },5);
  const auto r=orbitAveragedBmtAngularVelocities(a,L,m1,m2,mred,0.0,p); const auto q=orbitAveragedApsidalRate(a,L,m1,m2,mred,p);
  std::printf("BMT average (%d nodes) %.4e s; same with tiny dipoles %.4e s; apsidal rate (%d nodes) %.4e s\n",r.phaseNodes,tB,tT,q.phaseNodes,tA);
}
