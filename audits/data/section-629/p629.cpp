// Audit 629: Ps 1^3S1 - 2^3S1 from the engine's level functions (as fs429.cpp), split into the leading term and the part beyond Bohr.
// env G2 = 1 sets g = 2 for both particles (moments recomputed); CREM_LO_QED = 1 (read by the engine) drops the one-loop HFS.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  if(std::getenv("G2")){ firstGFactor=2.0; secondGFactor=2.0;
    firstMagneticMoment=0.5*2.0*magneton(firstSpecies); secondMagneticMoment=0.5*2.0*magneton(secondSpecies); }
  gSingletChannelEligible=true;
  const double hP=2*pi*hbar, MHz=1e-6/hP, a1=pairBohrRadius(activePair);
  const double mu1=firstMagneticMoment, mu2=secondMagneticMoment, s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1;
  const Vec3 ax{0,0,1}, P{1,0,0}, z{0,0,1};
  double F[3], U[3];
  for(int n: {1,2}){ const double a=a1*n*n, e=std::sqrt(1.0-0.25/(n*n));
    F[n]=levelFineStructureEnergy(a,0.5*hbar);
    U[n]=spinCouplingEnergy(a,e,ax,P,z*(s1*mu1),z*(s2*mu2),0.5*hbar);   // spins both +z: triplet
    std::printf("n=%d 3S1: F %.4f MHz, U_triplet %.4f MHz\n",n,F[n]*MHz,U[n]*MHz); }
  const double lead=0.75*pairBindingEnergy(activePair)*MHz;
  const double beyond=(F[2]+U[2]-F[1]-U[1])*MHz;
  std::printf("RESULT g %.12f %.12f lead %.4f beyond %.4f total %.4f MHz (alpha %.13e, m c^2/h %.10e Hz)\n",firstGFactor,secondGFactor,lead,beyond,lead+beyond,
              fineStructureConstant,firstMass*c*c/hP);
}
