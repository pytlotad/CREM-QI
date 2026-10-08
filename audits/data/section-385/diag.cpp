#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 385 diagnostic: exact s-state flow (step fraction f of pi) vs old substepped path, validation state; argv[1]=elapsed multiplier.
static double diff(const SecularSpinOrbitState& a,const SecularSpinOrbitState& b,const SecularSpinOrbitState& s0){
  const double g1=std::abs(firstGyromagneticRatioOf()),g2=std::abs(secondGyromagneticRatioOf());
  const double sc=s0.orbitalAngularMomentum.norm()+s0.firstDipole.norm()/g1+s0.secondDipole.norm()/g2;
  return ((a.orbitalAngularMomentum-b.orbitalAngularMomentum).norm()+(a.firstDipole-b.firstDipole).norm()/g1+(a.secondDipole-b.secondDipole).norm()/g2)/sc; }
int main(int argc,char**argv){
  const double mult=argc>1?std::atof(argv[1]):1.0;
  const double a=pairBohrRadius(activePair);
  const Vec3 d1{0.31,-0.47,0.826619}, d2{-0.62,0.19,0.761249}, L{0.08*hbar,-0.05*hbar,0.78*hbar};
  const SecularSpinOrbitState s0{L,d1*(firstMagneticMoment/d1.norm()),d2*(secondMagneticMoment/d2.norm()),0.0,orbitPlaneDirection(L,Vec3{0.71,0.43,-0.19})};
  const auto r=orbitAveragedBmtAngularVelocities(a,s0.orbitalAngularMomentum,s0.firstDipole,s0.secondDipole,pairReducedMass,0.0,s0.periapsisDirection);
  const double T=mult*0.2/std::max(r.first.norm(),r.second.norm());
  gExactSStateSpinFlow=false;
  std::vector<SecularSpinOrbitAdvance> old;
  for(int k=0;k<5;++k){ const double b=0.05/std::pow(4.0,k); old.push_back(advanceCoupledSecularSpinOrbit(s0,a,pairReducedMass,T,b,1<<26));
    std::printf("old bound %.3e substeps %d cos %.15f\n",b,old.back().substeps,dot(old.back().state.firstDipole,old.back().state.secondDipole)/(firstMagneticMoment*secondMagneticMoment)); std::fflush(stdout); }
  for(int k=0;k<4;++k) std::printf("old k=%d vs k=%d: %.3e\n",k,k+1,diff(old[k].state,old[k+1].state,s0));
  gExactSStateSpinFlow=true;
  for(int k=0;k<6;++k){ gExactPiFraction=std::pow(0.25,k); const auto x=advanceCoupledSecularSpinOrbit(s0,a,pairReducedMass,T,0.05,1<<26);
    std::printf("exact frac %.3e substeps %d vs oldfinest %.3e cos %.15f\n",gExactPiFraction,x.substeps,diff(x.state,old.back().state,s0),dot(x.state.firstDipole,x.state.secondDipole)/(firstMagneticMoment*secondMagneticMoment)); std::fflush(stdout);}
}
