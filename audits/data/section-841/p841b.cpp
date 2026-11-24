// Audit 841 (static probe P; two-body, from p839): Delta = Omega_ZD - Omega_Z = -<LRL rotation under a^Z(X) - a^ZD(X)> on the Kepler
// orbit of an 838 node (same frame / spins as p839).  ZD = CREM_SPIN_ORBIT_DIRAC (gSpinOrbitDirac toggled); parts: B = proper dipole
// only (switch off, dipoles replaced), F = factor only (switch off + spinOrbitDiracCorrection).  No trajectory integration.
// argv: pair rp e m1x m1y m1z m2x m2y m2z N
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  if(argc<11){ std::fprintf(stderr,"usage\n"); return 1; }
  applyPairFromOption(argv[1]); gMeasurementConservativeForce=true;
  const double rp=std::atof(argv[2]), ecc=std::atof(argv[3]);
  const Vec3 m1{std::atof(argv[4]),std::atof(argv[5]),std::atof(argv[6])}, m2{std::atof(argv[7]),std::atof(argv[8]),std::atof(argv[9])};
  const int N=std::atoi(argv[10]);
  const double M=firstMass+secondMass, mu=firstMass*secondMass/M, k=pairCoulombStrength, a=rp/(1-ecc);
  const double nmot=std::sqrt(k/(mu*a*a*a)), ec=std::sqrt(1-ecc*ecc);
  const Vec3 P{1,0,0}, Q{0,1,0};
  const auto acc=[&](const State& s,int form){   // 0 Z, 1 ZD, 2 B, 3 F
    State t=s; if(form==2){ t.firstDipole=s.firstProperDipole; t.secondDipole=s.secondProperDipole; }
    gSpinOrbitDirac=(form==1); MutualForces F=allExternalForces(t); gSpinOrbitDirac=false;
    if(form==3){ const MutualForces c3=spinOrbitDiracCorrection(s); F.first=F.first+c3.first; F.second=F.second+c3.second; }
    return relativisticAcceleration(s.firstVelocity,F.first,firstMass)-relativisticAcceleration(s.secondVelocity,F.second,secondMass); };
  double om[4]={0,0,0,0}, sw=0;
  for(int j=0;j<N;++j){
    const double E=2*pi*(j+0.5)/N, w=1-ecc*std::cos(E);
    const Vec3 r=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E)), v=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w);
    State s{}; s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M); s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
    s.firstProperDipole=m1; s.secondProperDipole=m2; s.firstDipole=m1; s.secondDipole=m2; synchronizeCovariantDipoles(s);
    const Vec3 a0=acc(s,0), hv=cross(r,v);
    for(int f=1;f<4;++f){ const Vec3 R=a0-acc(s,f), dA=cross(R,hv)+cross(v,cross(r,R)); om[f]+=-dot(Q,dA)/(k*ecc/mu)*w; }
    sw+=w;
  }
  std::printf("DELTA pair %s e %.12f ZD %.12e B %.12e F %.12e\n",argv[1],ecc,om[1]/sw,om[2]/sw,om[3]/sw);
}
