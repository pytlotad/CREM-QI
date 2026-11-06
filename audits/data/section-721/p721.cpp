#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 721: hydrogen level energies of the engine (beyond Bohr): S band F(a, L < hbar), P levels through levelEnergyBeyondBohr
// (import 21 with lambda = |g| - 1, electron spin parallel / antiparallel to L, minus the spin-coupling energy), triplet U of nS,
// the Bohr binding, constants.  Same reading as 429 (S, U) and 516/517 (P).  Output: key value (MHz unless stated).
int main(){
  applyPairFromOption("proton,electron");
  const double hP=2*pi*hbar, MHz=1e-6/hP, aB=pairBohrRadius(activePair);
  const bool firstLight=firstMass<secondMass;
  std::printf("alpha_inv %.9f\n",1.0/fineStructureConstant);
  std::printf("me_c2_MHz %.6f\n",electronMass*c*c*MHz);
  std::printf("mp_over_me %.9f\n",(firstLight?secondMass:firstMass)/electronMass);
  std::printf("binding_MHz %.6f\n",pairBindingEnergy(activePair)*MHz);
  std::printf("lightG %.14f\n",std::abs(firstLight?firstGFactor:secondGFactor));
  const double ml=firstLight?firstMagneticMoment:secondMagneticMoment, mh=firstLight?secondMagneticMoment:firstMagneticMoment;
  const double ql=firstLight?firstCharge:secondCharge, qh=firstLight?secondCharge:firstCharge;
  const Vec3 z{0,0,1}, x{1,0,0};
  for(int n=1;n<=4;++n){ const double a=n*n*aB, e=std::sqrt(1.0-0.25/(n*n));
    std::printf("F_%dS %.6f\n",n,levelFineStructureEnergy(a,0.5*hbar)*MHz);
    // triplet (F = 1, both spins +z): s = sign(q) m^ -> m = sign(q) |mu| s
    const Vec3 m1=z*((firstCharge<0?-1.0:1.0)*firstMagneticMoment), m2=z*((secondCharge<0?-1.0:1.0)*secondMagneticMoment);
    std::printf("U_%dS_F1 %.6f\n",n,spinCouplingEnergy(a,e,z,x,m1,m2,0.5*hbar)*MHz);
    for(int l=1;l<n;++l){ const double L=(l+0.5)*hbar, el=std::sqrt(1.0-std::pow((l+0.5)/n,2));
      for(int par:{-1,1}){ // light spin antiparallel (j = l - 1/2) / parallel (j = l + 1/2) to L
        const Vec3 mL=z*(par*(ql<0?-1.0:1.0)*ml), mH=x*((qh<0?-1.0:1.0)*mh);
        const Vec3 M1=firstLight?mL:mH, M2=firstLight?mH:mL;
        std::printf("F_%d%c%d/2 %.6f\n",n,"SPDFG"[l],2*l+par,(levelEnergyBeyondBohr(a,el,z,x,M1,M2,L)-spinCouplingEnergy(a,el,z,x,M1,M2,L))*MHz); } } }
}
