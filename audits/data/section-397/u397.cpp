#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 397: M1 flip energy in the 2s state straight from spinCouplingEnergy, against (4/3) mu0 m^2 n_QR(2).
int main(){
  const double mu=pairReducedMass, a=4.0*pairBohrRadius(activePair), L=0.5*hbar;
  const double e=std::sqrt(1.0-L*L/(mu*pairCoulombStrength*a));
  const Vec3 axis{0,0,1}, peri{1,0,0}, u=Vec3{0.3,0.4,std::sqrt(0.75)};
  const Vec3 m1=u*firstMagneticMoment;
  const double Uo=spinCouplingEnergy(a,e,axis,peri,m1,-u*secondMagneticMoment,L);  // antiparallel moments (o-Ps)
  const double Up=spinCouplingEnergy(a,e,axis,peri,m1,u*secondMagneticMoment,L);   // parallel (p-Ps)
  const double expect=(4.0/3.0)*mu0*firstMagneticMoment*secondMagneticMoment*quiggaRosnerContactDensity(a,L);
  const double h=2*pi*hbar;
  std::printf("2s: U_o - U_p = %.6f GHz, (4/3) mu0 m^2 n_QR(2) = %.6f GHz, ratio %.9f; p-Ps flip (U_p - U_o) = %.6f GHz\n",
    (Uo-Up)/h/1e9,expect/h/1e9,(Uo-Up)/expect,(Up-Uo)/h/1e9);
}
