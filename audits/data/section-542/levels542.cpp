#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 542: the model's 2^3P_J level energies (import 18 at the Born point, as in applyPositroniumSpinRule) relative to 1P_1 (= 0, rule B), GHz.
int main(){
  applyPairFromOption("electron,positron");
  const double a=4*pairBohrRadius(activePair), h=2*pi*hbar; const int l=1;
  const double r3=bornInverseCubeRadius(a,(l+0.5)*hbar), g=std::abs(firstGFactor);
  const double aSO=(2*g-1)*pairCoulombStrength*hbar*hbar*r3/(2*firstMass*firstMass*c*c), bT=1e-7*firstMagneticMoment*secondMagneticMoment*r3;
  for(int J=0;J<=2;++J){ const double LS=0.5*(J*(J+1.0)-l*(l+1.0)-2.0), Ll=l+0.5, L2=l*(l+1.0);
    const double P=2.0/3.0-(LS*(LS+0.5)-2.0*L2/3.0)/(2.0*(Ll-1.0)*(Ll+1.0));
    std::printf("J %d E %.6f GHz\n",J,(aSO*LS+bT*(6*P-4))/h*1e-9); }
  std::printf("a %.6e m, g %.10f, muB-ish |m1| %.6e\n",a,g,firstMagneticMoment);
}
