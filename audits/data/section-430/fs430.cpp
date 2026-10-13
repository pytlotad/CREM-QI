#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 430: 2^3P_J spin part -- model (363) vs QM <r^-3> and QM tensor spectrum; fit to measured intervals.
int main(){
  const double hP=2*pi*hbar, GHz=1e-9/hP, mu=firstMagneticMoment, aPs=pairBohrRadius(activePair);
  const int n=2, l=1;
  const double r3QM=1.0/(n*n*n*aPs*aPs*aPs*l*(l+0.5)*(l+1)), r3Langer=1.0/(n*n*n*aPs*aPs*aPs*std::pow(l+0.5,3));
  // model 363: E_SO(J=2) = +3.6554 GHz at <L.S>_cl = +1.5  ->  a_mod
  const double aMod=3.6554/1.5, tMod[3]={-0.4065,0.1355,-0.4065};
  const double LSqm[3]={-2,-1,1}, LScl[3]={-1.5,-0.5,1.5}, S12[3]={-4,2,-0.4};
  const double FP=levelFineStructureEnergy(n*n*aPs,1.5*hbar)*GHz, ES1=7.4220, E1P1=FP;
  const double bQM=1e-7*mu*mu*r3QM*GHz;
  std::printf("<r^-3>: QM %.6e, Langer %.6e, ratio %.6f;  a_mod %.4f GHz, 9/8 a_mod %.4f;  b_QM = (mu0/4pi) mu^2 <r^-3>_QM = %.4f GHz\n",r3QM,r3Langer,r3Langer/r3QM,aMod,aMod*9/8,bQM);
  const double meas[3]={18.49965,13.01242,8.62438}, d10=meas[0]-meas[1], d21=meas[1]-meas[2];
  // fit: P1-P0 = a + 6 b ; P2-P1 = 2 a - 2.4 b
  const double bFit=(2*d10-d21)/14.4, aFit=d10-6*bFit;
  std::printf("fit to measured (P1-P0 %.4f, P2-P1 %.4f): a = %.4f GHz (9/8 a_mod x %.5f), b = %.4f GHz (b_QM x %.5f)\n",d10,d21,aFit,aFit/(aMod*9/8),bFit,bFit/bQM);
  for(int v=0;v<3;++v){
    double E[3];
    for(int J=0;J<3;++J) E[J]=FP+(v==0?aMod*LScl[J]+tMod[J]:(v==1?aMod*9/8*LScl[J]+tMod[J]:aMod*9/8*LSqm[J]+bQM*S12[J]));
    std::printf("V%d: P1-P0 %.3f P2-P1 %.3f GHz;",v,E[1]-E[0],E[2]-E[1]);
    for(int J=0;J<3;++J) std::printf("  nu_%d %+.3f (res %+.3f)",J,ES1-E[J],ES1-E[J]-meas[J]);
    std::printf("\n"); }
  std::printf("2^3S1-2^1P1: %+.3f measured 11.180 residual %+.3f GHz\n",ES1-E1P1,ES1-E1P1-11.180);
}
