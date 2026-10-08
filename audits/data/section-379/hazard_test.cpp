#include "modules/contact_annihilation.hpp"
#include <cstdio>
#include <random>
// Audit 379: drawDecayFromHazard against an exact survival.  Hazard G1 on [0, t1), 0 on [t1, t2) (a state without
// contact), G2 on [t2, ts), then Gf past the stop ts.  Exact S(t) piecewise exponential; compare P(decay in cascade),
// mean, and KS distance of 1e6 draws.
int main(){
  const double G1=1/3.36e-9, t1=1.26e-9, t2=4.0e-9, G2=1/1e-9, ts=4.62e-9, Gf=1/124.494e-12;
  const std::vector<std::array<double,3>> seg={{0,t1,G1},{t1,t2,0.0},{t2,ts,G2}};
  const auto H=[&](double t){ if(t<t1) return G1*t; if(t<t2) return G1*t1; if(t<ts) return G1*t1+G2*(t-t2); return G1*t1+G2*(ts-t2)+Gf*(t-ts); };
  const double Hs=H(ts), pCascade=1-std::exp(-Hs);
  // exact mean = int S dt
  const double mean=(1-std::exp(-G1*t1))/G1 + std::exp(-G1*t1)*(t2-t1) + std::exp(-G1*t1)*(1-std::exp(-G2*(ts-t2)))/G2 + std::exp(-Hs)/Gf;
  std::mt19937_64 r(42); std::uniform_real_distribution<double> U(0,1);
  const int N=1000000; std::vector<double> x(N); int inC=0; double sum=0;
  for(int i=0;i<N;++i){ double u=U(r); while(!(u>0)) u=U(r); bool c=false; x[i]=drawDecayFromHazard(seg,ts,Gf,u,c); inC+=c; sum+=x[i]; }
  std::sort(x.begin(),x.end()); double D=0; int inGap=0;
  for(int i=0;i<N;++i){ const double F=1-std::exp(-H(x[i])); D=std::max(D,std::max(F-(double)i/N,(i+1.0)/N-F)); if(x[i]>t1&&x[i]<t2) ++inGap; }
  std::printf("P(decay in cascade): drawn %.5f exact %.5f (sd %.5f)\n",inC/(double)N,pCascade,std::sqrt(pCascade*(1-pCascade)/N));
  std::printf("mean: drawn %.6e exact %.6e (rel %.2e)\n",sum/N,mean,sum/N/mean-1);
  std::printf("KS distance D = %.2e (1e6 draws; 1.36/sqrt(N) = %.2e at 5%%)\n",D,1.36/std::sqrt((double)N));
  std::printf("decays inside the zero-hazard gap [t1,t2): %d (exact 0)\n",inGap);
}
