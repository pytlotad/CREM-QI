#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 493: 2^3P_J energies from the engine with import 18 and the latched spin cone (import 19, CREM_SPIN_CONE_TENSOR);
// env PAIR and the import-19 switches.  Same states as fs431; first a code test of latchedSpinConeTensor against s12_zatrzask.txt.
int main(){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gSingletChannelEligible=true;
  const double r3q=std::sqrt(3.0), r2=std::sqrt(2.0);
  std::printf("code test (alpha 0 / 90 deg): qm J0 %+.3f %+.3f  J1 %+.3f %+.3f  langer J0 %+.3f %+.3f  J1 %+.3f %+.3f  J2 %+.3f %+.3f\n",
    latchedSpinConeTensor(r3q,r2,-1.0,0),latchedSpinConeTensor(r3q,r2,-1.0,pi/2),
    latchedSpinConeTensor(r3q,r2,-0.5,0),latchedSpinConeTensor(r3q,r2,-0.5,pi/2),
    latchedSpinConeTensor(2,1.5,(0.25-4.5)/4.5,0),latchedSpinConeTensor(2,1.5,(0.25-4.5)/4.5,pi/2),
    latchedSpinConeTensor(2,1.5,(2.25-4.5)/4.5,0),latchedSpinConeTensor(2,1.5,(2.25-4.5)/4.5,pi/2),
    latchedSpinConeTensor(2,1.5,(6.25-4.5)/4.5,0),latchedSpinConeTensor(2,1.5,(6.25-4.5)/4.5,pi/2));
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
}
