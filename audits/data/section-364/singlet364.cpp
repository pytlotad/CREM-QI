#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 364: 2^1P1 in the classical secular Hamiltonian -- E_rel(2P) + orientation-averaged tensor of the classical singlet.
int main(){
  const double hP=2*pi*hbar, m=firstMass, k=pairCoulombStrength, mu=firstMagneticMoment, GHz=1e-9/hP;
  const double alpha=k/(hbar*c), a4=std::pow(alpha,4)*m*c*c, aPs=pairBohrRadius(activePair), psi0=1.0/(pi*aPs*aPs*aPs);
  const int n=2; const double a=n*n*aPs;
  auto Erel=[&](double Lh){ return a4*(11.0/(64.0*n*n*n*n)-1.0/(4.0*n*n*n*Lh)); };
  const double eP=std::sqrt(1-std::pow(1.5/n,2)); const Vec3 Lh{0,0,1}, P{1,0,0};
  auto U=[&](const Vec3& s){ return orbitAveragedDipoleEnergy(a,eP,Lh,P,s*mu,s*mu); };   // singlet: moments parallel
  const double Uiso=(U({1,0,0})+U({0,1,0})+U({0,0,1}))/3.0, Ut=0.0;  // tensor trace-free: orientation average = 0
  const double E1P1=Erel(1.5)+Ut;                 // + no spin-orbit (S = 0), no contact (l = 1)
  const double E3S1=Erel(0.5)+(2*mu0/3)*mu*mu*psi0/(n*n*n);
  std::printf("E_rel(2P) = %.4f GHz; singlet isotropic part (contact, Plummer orbit) = %.2e GHz (dropped: QR contact is 0 for l = 1)\n",Erel(1.5)*GHz,Uiso*GHz);
  std::printf("E(2^3S1) = %.4f GHz, E(2^1P1) = %.4f GHz\n",E3S1*GHz,E1P1*GHz);
  const double nu=(E3S1-E1P1)*GHz, meas=11.180;
  std::printf("2^3S1 - 2^1P1: model %+.4f GHz, measured %+.3f GHz, deficit of 2^3S1 = %.4f GHz\n",nu,meas,meas-nu);
  const double model363[3]={-6.879,-9.858,-14.190}, meas363[3]={18.49965,13.01242,8.62438};
  for(int J=0;J<=2;++J) std::printf("nu_%d residual after removing the S deficit: %+.3f GHz (part attributable to 2^3P_%d spin terms)\n",J,model363[J]+(meas-nu)-meas363[J],J);
}
