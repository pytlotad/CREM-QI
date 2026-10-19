#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 517: levelFineStructureEnergy over the S band and across L = hbar for H, Mu, Ps, n = 1, 2, 3, 5.  Run three times: default
// (k = 1 in the S band), CREM_NO_SOMMERFELD_S=1 (import 17), CREM_NO_S_CONTACT_ENERGY=1 (bare E_rel(L)); ana517.py compares.
// Also the 2P1/2 level through levelEnergyBeyondBohr (import 21, electron spin antiparallel to L) for H and Mu.
int main(){
  const double Ls[]={0.05,0.2,0.35,0.5,0.65,0.75,0.8,0.95,0.999999,1.0,1.000001,1.5,2.0};
  for(const char* p:{"proton,electron","antimuon,electron","electron,positron"}){
    applyPairFromOption(p); const double aB=pairBohrRadius(activePair), MHz=1e-6/(2*pi*hbar);
    for(int n:{1,2,3,5}){ const double a=n*n*aB;
      for(double L:Ls){ if(L>n-0.5+1e-9&&L>1.000001+1e-9) continue; /* n = 1: L up to hbar+ only for the edge (formal) */ std::printf("F %s %d %.6f %.12f\n",p,n,L,levelFineStructureEnergy(a,L*hbar)*MHz); }
      if(n>=2&&!particleAntiparticlePair()){
        const Vec3 z{0,0,1},x{1,0,0}; const Vec3 mn=z*(firstMass<secondMass?-firstMagneticMoment:secondMagneticMoment)*0.0;
        const bool firstLight=firstMass<secondMass; const double ql=firstLight?firstCharge:secondCharge, mu=firstLight?firstMagneticMoment:secondMagneticMoment;
        const Vec3 ml=z*(-(ql<0?-1.0:1.0)*mu), mh=x*(firstLight?secondMagneticMoment:firstMagneticMoment);   // light spin antiparallel to L
        const Vec3 m1=firstLight?ml:mh, m2=firstLight?mh:ml; const double e=std::sqrt(1-std::pow(1.5/n,2));
        std::printf("P12 %s %d %.12f\n",p,n,(levelEnergyBeyondBohr(a,e,z,x,m1,m2,1.5*hbar)-spinCouplingEnergy(a,e,z,x,m1,m2,1.5*hbar))*MHz); (void)mn; }
    }
  }
}
