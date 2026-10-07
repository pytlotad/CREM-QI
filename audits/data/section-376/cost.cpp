#include "modules/crem_collapse.hpp"
#include <chrono>
#include <cstdio>
int main(){
  const double mred=firstMass*secondMass/(firstMass+secondMass), a=pairBohrRadius(activePair);
  const Vec3 L{0,0,0.5*hbar};
  const Vec3 d1=Vec3{0.6,0,0.8}*firstMagneticMoment, d2=Vec3{0,0.6,0.8}*secondMagneticMoment;
  SecularSpinOrbitState st{L,d1,d2,0.0,Vec3{1,0,0}};
  const auto r=orbitAveragedBmtAngularVelocities(a,L,d1,d2,mred,0.0,Vec3{1,0,0});
  std::printf("rates |w1| %.4e |w2| %.4e rad/s\n",r.first.norm(),r.second.norm());
  for(double dt: {1e-12,1e-11}){
    auto t0=std::chrono::steady_clock::now();
    const auto x=advanceCoupledSecularSpinOrbit(st,a,mred,dt);
    double s=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::printf("dt %.1e: completed %d, %.3f s, w %.6f\n",dt,x.completed,s,contactTwoPhotonWeight(x.state.firstDipole,x.state.secondDipole));
  }
}
