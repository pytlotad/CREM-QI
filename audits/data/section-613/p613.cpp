// Audit 613 T4: engine downChannelShares (normalized) for hydrogen p, d, f states.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  const int st[8][2]={{3,1},{4,1},{5,1},{6,1},{4,2},{5,2},{6,2},{6,3}};
  for(auto& x: st){ const int n=x[0], l=x[1]; auto w=downChannelShares(n,l); double t=0; for(double v: w) t+=v;
    std::printf("%d%c:",n,"spdf"[l]); for(size_t i=0;i<w.size();++i) std::printf(" n'=%zu %.8f",l+i,w[i]/t); std::printf("\n"); }
}
