#include "modules/crem_collapse.hpp"
#include "modules/qed_reference.hpp"
#include <cstdio>
// Audit 392: Bhabha (tree QED) / Rutherford at K_CM = 20 eV for e+e-, and the lab energy.
int main(){
  const double K=20.0*eCharge, m=electronMass, k=pairCoulombStrength, lC=k/(2.0*K);
  const double sqrtS=2*m*c*c+K, Tlab=(sqrtS*sqrtS)/(2*m*c*c)-2*m*c*c;
  std::printf("T_lab = %.7f eV (2K + K^2/(2mc^2) = %.7f eV)\n",Tlab/eCharge,(2*K+K*K/(2*m*c*c))/eCharge);
  double worst=0;
  for(double deg: {5.0,10.0,20.0,45.0,90.0,135.0,170.0}){ const double t=deg*pi/180;
    const double ruth=lC*lC/std::pow(std::sin(0.5*t),4)/4.0*4.0/4.0; // l_C^2/(4 sin^4) form below
    const double r=lC*lC/(4.0*std::pow(std::sin(0.5*t),4));
    const double q=qedElasticDifferentialCrossSection(K,t,m,m);
    worst=std::max(worst,std::abs(q/r-1)); (void)ruth;
    std::printf("theta_CM %6.1f deg  Bhabha/Rutherford - 1 = %+.3e\n",deg,q/r-1); }
  std::printf("worst %.3e\n",worst);
}
