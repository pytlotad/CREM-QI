// Audit 681 smoke: pair radii and stop radii per selectable pair.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  const char* pairs[]={"electron,positron","proton,electron","muon,antimuon","proton,antiproton","electron,antimuon"};
  for(const char* p:pairs){ applyPairFromOption(p);
    std::printf("%-20s rAnn %.17e (== r_e bitwise %d) Rstop %.6e m  a_pair %.6e m  Rstop==r* %d\n",p,pairAnnihilationRadius(),
      pairAnnihilationRadius()==classicalElectronRadius,collapseStopRadius(),pairBohrRadius(activePair),collapseStopRadius()==comptonBarrierRadius);}
}
