#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 431: 2^3P_J energies from the engine with import 18; n = 2 intervals; env PAIR.
int main(){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gSingletChannelEligible=true;
  const double hP=2*pi*hbar, GHz=1e-9/hP, aB=pairBohrRadius(activePair), mu1=firstMagneticMoment, mu2=secondMagneticMoment;
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, Lm=1.5, Sm=1.0, a=4*aB, eP=std::sqrt(1-std::pow(Lm/2,2));
  const Vec3 z{0,0,1}, x{1,0,0};
  const double FP=levelFineStructureEnergy(a,Lm*hbar);
  double E[3];
  for(int J=0;J<3;++J){ const double Jm=J+0.5, ct=(Jm*Jm-Lm*Lm-Sm*Sm)/(2*Lm*Sm), st=std::sqrt(std::max(0.0,1-ct*ct));
    const Vec3 sh=x*st+z*ct;
    const double U=spinCouplingEnergy(a,eP,z,x,sh*(s1*mu1),sh*(s2*mu2),Lm*hbar);
    E[J]=(FP+U)*GHz; std::printf("2^3P_%d: cos %+.4f  U %+.4f GHz  E %+.4f GHz\n",J,ct,U*GHz,E[J]); }
  const double eS=std::sqrt(1-1.0/16), ES1=(levelFineStructureEnergy(a,0.5*hbar)+spinCouplingEnergy(a,eS,z,x,z*(s1*mu1),z*(s2*mu2),0.5*hbar))*GHz;
  const double meas[3]={18.49965,13.01242,8.62438};
  std::printf("E(2^3S1) %.4f GHz; P1-P0 %.4f (meas 5.4872) P2-P1 %.4f (meas 4.3880)\n",ES1,E[1]-E[0],E[2]-E[1]);
  for(int J=0;J<3;++J) std::printf("nu_%d %+.4f measured %.5f residual %+.4f GHz\n",J,ES1-E[J],meas[J],ES1-E[J]-meas[J]);
  std::printf("2^3S1-2^1P1 %+.4f measured 11.180 residual %+.4f GHz\n",ES1-FP*GHz,ES1-FP*GHz-11.180);
  const double r3=1.0/(aB*aB*aB*8*1*1.5*2);
  std::printf("a = %.5f GHz, b = %.5f GHz, g = %.8f\n",(2*std::abs(firstGFactor)-1)*pairCoulombStrength*hbar*hbar*r3/(2*firstMass*firstMass*c*c)*GHz,1e-7*mu1*mu2*r3*GHz,firstGFactor);
}
