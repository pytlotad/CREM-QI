// Audit 734: C++ EBK radial elements vs ana733 (H734a).
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  const int T[7][4]={{2,1,1,0},{3,1,1,0},{3,1,2,0},{3,0,2,1},{3,2,2,1},{4,1,1,0},{4,3,3,2}};
  const double py[7]={1.285414,0.514370,3.057997,0.943317,4.742840,0.303134,10.225021};
  for(int i=0;i<7;++i){ const double R=ebkRadialElement(T[i][0],T[i][1],T[i][2],T[i][3]);
    std::printf("%d%d -> %d%d: C++ %.6f  ana733 %.6f  rel %+.2e\n",T[i][0],T[i][1],T[i][2],T[i][3],R,py[i],R/py[i]-1); }
  applyPairFromOption("proton,electron");
  std::printf("hydrogen A(2p->1s) = %.6e 1/s -> tau %.5f ns\n",ebkE1Rate(2,1,1,0),1e9/ebkE1Rate(2,1,1,0));
}
