// Audit 640: Ps n = 2 levels from the engine's level functions (as ../section-431/fs431.cpp), absolute energies beyond Bohr [MHz]:
// 2^3S1, 2^3P_J (J = 0, 1, 2), 2^1P1 (F only).  env G2 = 1: g = 2 for both particles (moments recomputed); CREM_LO_QED = 1 (engine) drops one loop.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  if(std::getenv("G2")){ firstGFactor=2.0; secondGFactor=2.0;
    firstMagneticMoment=0.5*2.0*magneton(firstSpecies); secondMagneticMoment=0.5*2.0*magneton(secondSpecies); }
  gSingletChannelEligible=true;
  const double MHz=1e-6/(2*pi*hbar), aB=pairBohrRadius(activePair), mu1=firstMagneticMoment, mu2=secondMagneticMoment;
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, a=4*aB;
  const Vec3 z{0,0,1}, x{1,0,0};
  const double eS=std::sqrt(1-1.0/16), eP=std::sqrt(1-std::pow(1.5/2,2));
  const double FS=levelFineStructureEnergy(a,0.5*hbar), US=spinCouplingEnergy(a,eS,z,x,z*(s1*mu1),z*(s2*mu2),0.5*hbar);
  const double FP=levelFineStructureEnergy(a,1.5*hbar);
  std::printf("2^3S1 %.4f (F %.4f U %.4f)\n2^1P1 %.4f\n",(FS+US)*MHz,FS*MHz,US*MHz,FP*MHz);
  for(int J=0;J<3;++J){ const double Jm=J+0.5, ct=(Jm*Jm-2.25-1)/(2*1.5), st=std::sqrt(std::max(0.0,1-ct*ct)); const Vec3 sh=x*st+z*ct;
    const double U=spinCouplingEnergy(a,eP,z,x,sh*(s1*mu1),sh*(s2*mu2),1.5*hbar);
    std::printf("2^3P%d %.4f (U %.4f)\n",J,(FP+U)*MHz,U*MHz); }
}
