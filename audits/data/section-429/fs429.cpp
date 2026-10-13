#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 429: level fine-structure energy F (E_rel + import 17 contact) from the engine; Dirac check; 1S-2S. env PAIR.
int main(){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  gSingletChannelEligible=true;
  const double hP=2*pi*hbar, GHz=1e-9/hP, m1=firstMass, m2=secondMass, mu=m1*m2/(m1+m2);
  const double alpha=pairCoulombStrength/(hbar*c), a1=pairBohrRadius(activePair);
  const double mu1=firstMagneticMoment, mu2=secondMagneticMoment, s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1;
  const Vec3 ax{0,0,1}, P{1,0,0}, z{0,0,1};
  double FS[3], UT[3];
  for(int n: {1,2}){ const double a=a1*n*n, e=std::sqrt(1.0-0.25/(n*n));
    FS[n]=levelFineStructureEnergy(a,0.5*hbar);
    UT[n]=spinCouplingEnergy(a,e,ax,P,z*(s1*mu1),z*(s2*mu2),0.5*hbar);   // spins both +z: triplet
    const double dirac=-(mu*c*c*std::pow(alpha,4)/(2.0*std::pow(n,4)))*(n-0.75);
    std::printf("n=%d S: F %.6f GHz (Dirac reduced-mass %.6f, ratio %.6f)  U_triplet %.6f GHz\n",n,FS[n]*GHz,dirac*GHz,FS[n]/dirac,UT[n]*GHz); }
  std::printf("2P (L = 1.5 hbar): F %.6f GHz\n",levelFineStructureEnergy(4*a1,1.5*hbar)*GHz);
  const double lead=0.75*pairBindingEnergy(activePair)/hP;   // Hz
  const double alpha4=(FS[2]+UT[2]-FS[1]-UT[1])/hP;
  const bool ps=std::abs(m1-m2)<1e-9*m1;
  const double measured=ps?1233607216.4e6:2466061413187.035e3;
  std::printf("1S-2S (triplet / F=1): leading %.6f GHz, beyond Bohr %+.6f GHz, model %.6f GHz, measured %.6f GHz, residual %+.4f GHz\n",
    lead*1e-9,alpha4*1e-9,(lead+alpha4)*1e-9,measured*1e-9,(lead+alpha4-measured)*1e-9);
}
