#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 358: the validation's para libration (radius 0.2 a, L = sqrt(0.2) hbar, moments tilted 0.6,0,0.8), raw transport
// (s-state isotropy off), cos range over 48 samples to 24/speed, for several substeps; run under both schemes via env.
int main(){
  gSStateIsotropyEnabled=false;
  const double mred=firstMass*secondMass/(firstMass+secondMass), radius=0.2*pairBohrRadius(activePair);
  const Vec3 orbital{0,0,std::sqrt(mred*pairCoulombStrength*radius)}, tilt{0.6,0,0.8};
  const SecularSpinOrbitState st{orbital,tilt*(firstMagneticMoment/tilt.norm()),tilt*(secondMagneticMoment/tilt.norm()),0.0,Vec3{1,0,0}};
  const auto r=orbitAveragedBmtAngularVelocities(radius,orbital,st.firstDipole,st.secondDipole,mred,0.0,Vec3{1,0,0});
  const double speed=std::max(r.first.norm(),r.second.norm());
  for(double sub: {0.05,0.0125,0.003125}){
    double cmin=2,cmax=-2,Lmin=1e9,Lmax=-1;
    for(int step=1;step<=48;++step){
      const auto a=advanceCoupledSecularSpinOrbit(st,radius,mred,0.5*step/speed,sub);
      if(!a.completed){ std::printf("incomplete\n"); continue; }
      const double c=dot(a.state.firstDipole,a.state.secondDipole)/(a.state.firstDipole.norm()*a.state.secondDipole.norm());
      cmin=std::min(cmin,c); cmax=std::max(cmax,c);
      const double L=a.state.orbitalAngularMomentum.norm()/hbar; Lmin=std::min(Lmin,L); Lmax=std::max(Lmax,L);
    }
    std::printf("substep %.6f: cos %.6f .. %.6f   |L| %.6f .. %.6f hbar\n",sub,cmin,cmax,Lmin,Lmax);
  }
}
