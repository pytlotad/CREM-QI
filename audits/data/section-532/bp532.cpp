#include "modules/crem_collapse.hpp"
#include <cstdio>
// import 21 spin-orbit lambda [F(k) - F(L_l)] vs Breit-Pauli xi <L.S> with Born <r^-3>, hydrogen and muonium
int main(){
  for(const char* p:{"proton,electron","antimuon,electron"}){
    applyPairFromOption(p); const double aB=pairBohrRadius(activePair), MHz=1e-6/(2*pi*hbar);
    const double me=secondMass, mh=firstMass, g=std::abs(secondGFactor), k=pairCoulombStrength;
    double w1=0,w2=0,w3=0;
    for(int n=2;n<=8;++n) for(int l=1;l<n;++l) for(int s=-1;s<=1;s+=2){
      const double a=n*n*aB, Ll=(l+0.5)*hbar, kk=(s>0?l+1:l)*hbar, j=l+0.5*s;
      const double F21=(g-1)*(levelFineStructureEnergy(a,kk)-levelFineStructureEnergy(a,Ll));
      const double LS=0.5*(j*(j+1)-l*(l+1.0)-0.75)*hbar*hbar, r3=bornInverseCubeRadius(a,Ll);
      const double bp0=(k/(c*c))*((g-1)/(2*me*me))*r3*LS;                 // own spin-orbit with Thomas only
      const double bp1=(k/(c*c))*((g-1)/(2*me*me)+g/(2*me*mh))*r3*LS;     // + recoil (other-orbit) term, as in the transport (467)
      const double bp2=(k/(c*c))*(1.0/(2*me*me))*r3*LS;                   // g = 2 limit of own term (check of F(k)-F(L_l) itself)
      const double F21g2=(levelFineStructureEnergy(a,kk)-levelFineStructureEnergy(a,Ll));
      w1=std::max(w1,std::abs(F21/bp0-1)); w2=std::max(w2,std::abs(F21/bp1-1)); w3=std::max(w3,std::abs(F21g2/bp2-1));
      if(n==2&&p[0]=='p') std::printf("  %s n 2 l %d j %.1f: import 21 %+.4f MHz, BP own %+.4f, BP own+recoil %+.4f\n",p,l,j,F21*MHz,bp0*MHz,bp1*MHz); }
    std::printf("%s n = 2..8: max |import21/BP_own - 1| %.3e; |import21/BP_own+recoil - 1| %.3e; g = 2: |[F(k)-F(L_l)]/BP - 1| %.3e  (m_l/m_h %.3e)\n",p,w1,w2,w3,me/mh);
  }
}
