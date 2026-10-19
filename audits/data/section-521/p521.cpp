#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 521: spinRadialProjectionSquared (engine) vs QM for spin 1 (<S12>_J) and spin 1/2 ((S.r^)^2 = 1/4), and the classical in-plane average.
int main(){
  double w1=0,w2=0;
  std::printf("spin 1: <S12>_J = 6 P - 4\n");
  for(int l=1;l<=6;++l) for(int J=l+1;J>=l-1;--J){
    const double SL=0.5*(J*(J+1.0)-l*(l+1.0)-2.0), qm=J==l+1?-2.0*l/(2*l+3.0):(J==l?2.0:-2.0*(l+1.0)/(2*l-1.0));
    const double rule=6*spinRadialProjectionSquared(l,2.0,SL)-4, cl=6*(2.0/3-(SL*SL-2.0*l*(l+1)/3)/(2.0*l*(l+1)))-4;
    w1=std::max(w1,std::abs(rule-qm));
    if(l<=3) std::printf("  l %d J %d: QM %+.6f rule %+.6f classical (Born squares) %+.6f\n",l,J,qm,rule,cl); }
  for(int l=1;l<=10;++l) for(double j:{l-0.5,l+0.5}){ const double SL=0.5*(j*(j+1)-l*(l+1.0)-0.75); w2=std::max(w2,std::abs(spinRadialProjectionSquared(l,0.75,SL)-0.25)); }
  std::printf("max |rule - QM| spin 1, l = 1..6: %.2e; spin 1/2 max |(S.r)^2 - 1/4|, l = 1..10: %.2e\n",w1,w2);
}
