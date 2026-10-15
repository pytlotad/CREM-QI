// Audit 460: spin dephasing by spin-orbit (common precession) and by the partner's field (differential) vs state lifetime, Ps.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("electron,positron");
  const double mu=firstMass*secondMass/(firstMass+secondMass), aB=pairBohrRadius(activePair), k=pairCoulombStrength;
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, M1=firstMagneticMoment, M2=secondMagneticMoment;
  // hydrogen QM lifetimes from audit 432 (infinite mass), Ps x 2
  struct St{int n,l; double tauH;}; St st[]={{2,1,1.595e-9},{3,1,5.269e-9},{3,2,15.459e-9},{4,1,12.298e-9},{4,2,36.130e-9},{4,3,72.488e-9},{1,0,0},{2,0,0}};
  std::printf("# n l  Omega_c[rad/s]  Omega_d(triplet)  Omega_d(singlet)  tau_Ps[ns]  N_c  N_d(trip)  N_d(sing)\n");
  for(auto s: st){
    const double a=aB*s.n*s.n, L=(s.l+0.5)*hbar; const Vec3 z{0,0,1}, x{1,0,0};
    const double th=M_PI/3; const Vec3 sh{std::sin(th),0,std::cos(th)};               // spin 60 deg from L
    // triplet: both spins along sh; moments m_i = sign(q_i) mu_i sh
    const auto rt=orbitAveragedBmtAngularVelocities(a,z*L,sh*(s1*M1),sh*(s2*M2),mu,0.0,x);
    // singlet: spins antiparallel (+sh, -sh)
    const auto rs=orbitAveragedBmtAngularVelocities(a,z*L,sh*(s1*M1),sh*(-s2*M2),mu,0.0,x);
    const double Oc=((rt.first+rt.second)*0.5).norm(), Odt=(rt.first-rt.second).norm(), Ods=(rs.first-rs.second).norm();
    const double tau=2*s.tauH;
    if(s.l==0) std::printf("%d %d  %.4e  %.4e  %.4e  (s-state: isotropy rule zeroes the orbital part in the transport)\n",s.n,s.l,Oc,Odt,Ods);
    else std::printf("%d %d  %.4e  %.4e  %.4e  %.3f  %.2f  %.2f  %.2f\n",s.n,s.l,Oc,Odt,Ods,tau*1e9,Oc*tau/(2*pi),Odt*tau/(2*pi),Ods*tau/(2*pi));
  }
}
