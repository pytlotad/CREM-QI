// Audit 624, H624b: the rule functions at given (level, L/hbar) -- reads "n L" pairs on stdin, prints the s-state E1 total, the
// s-state E2 total and the up-channel total (pair units of the chosen PAIR), as the mechanical path computes them.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  double n, L;
  while(std::scanf("%lf %lf",&n,&L)==2){
    double sE1=0.0, sE2=0.0, up=0.0;
    if(L<1.0-1.0e-9){ for(double r: sStateChannelRates(n)) sE1+=r; for(double r: sStateQuadrupoleChannelRates(n)) sE2+=r; }
    else for(double r: upChannelRates(n,L)) up+=r;
    std::printf("%.9f %.9f %.12e %.12e %.12e\n",n,L,sE1,sE2,up);
  }
}
