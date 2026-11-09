// Audit 763: energy probe -- j-splitting of the hydrogen level energy (sommerfeldLevelFineStructureEnergy) at k = l + 1 and k = l,
// n = 2..4, l = 1..n-1; spin parallel / antiparallel to L.  Run with and without CREM_LS_FROM_PRECESSION=1.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("proton,electron");
  const double h=2*pi*hbar, aP=pairBohrRadius(activePair); const Vec3 z{0,0,1};
  for(int n=2;n<=4;++n) for(int l=1;l<n;++l){
    const double a=n*n*aP, L=(l+0.5)*hbar; double E[2]; int K[2];
    for(int s=0;s<2;++s){ const Vec3 m=z*(s==0?1.0:-1.0)*9.2740100783e-24;
      K[s]=(int)std::lround(sommerfeldLevelAngularMomentum(z,m,m,L)/hbar); E[s]=sommerfeldLevelFineStructureEnergy(a,z,m,m,L); }
    const int hi=K[0]>K[1]?0:1;
    std::printf("n %d l %d k %d %d  E(k=l+1) - E(k=l) = %.6f MHz\n",n,l,K[hi],K[1-hi],(E[hi]-E[1-hi])/h/1e6);
  }
}
