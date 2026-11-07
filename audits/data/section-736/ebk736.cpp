// Audit 736: C++ E2 (r^2 EBK element, angular sum) and M1 (spin wave) rates vs ana735 (H736a).
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("proton,electron");
  const int T[6][4]={{3,2,1,0},{4,2,1,0},{4,2,2,0},{4,3,2,1},{5,4,3,2},{4,2,3,2}};
  const double r2[6]={1.72349,1.11905,7.87815,24.08240,111.90639,92.76444};
  const double ang[6]={1.20000,1.20000,1.20000,1.54286,1.71429,1.71429};
  const double A[6]={5.8597e+02,3.2239e+02,5.1130e+00,6.1428e+01,1.1564e+01,1.1862e+00};
  bool okR=true, okG=true, okA=true;
  for(int i=0;i<6;++i){ const double R=ebkRadialElement(T[i][0],T[i][1],T[i][2],T[i][3],2);
    const double G=quadrupoleAngularSum(T[i][1],T[i][3]); const double a=ebkE2Rate(T[i][0],T[i][1],T[i][2],T[i][3]);
    okR&=std::abs(R/r2[i]-1)<1e-5; okG&=std::abs(G-ang[i])<1e-5; okA&=std::abs(a/A[i]-1)<1e-3;
    std::printf("%d%d -> %d%d: R2 C++ %.6f ana735 %.5f rel %+.2e | G %.6f / %.5f | A_E2 C++ %.5e ana735 %.4e rel %+.2e\n",
      T[i][0],T[i][1],T[i][2],T[i][3],R,r2[i],R/r2[i]-1,G,ang[i],a,A[i],a/A[i]-1); }
  // exact G values (6 (2l'+1) (l 2 l';000)^2) for the record
  std::printf("G exact: 2->0 %.12f (6/5), 3->1 %.12f (54/35), 4->2 %.12f (12/7), 2->2 %.12f (12/7), 0->2 %.12f (6)\n",
    quadrupoleAngularSum(2,0),quadrupoleAngularSum(3,1),quadrupoleAngularSum(4,2),quadrupoleAngularSum(2,2),quadrupoleAngularSum(0,2));
  // E2 channels of 3d, 4f and the s-state channels of 5s
  bool ok=false; const auto ch=ebkQuadrupoleChannels(3.0,2.5,ok); double t=0; for(auto& x: ch) t+=x.rate;
  std::printf("3d E2 channels (ok %d):",int(ok)); for(auto& x: ch) std::printf(" ->%d,%d %.4e",x.n2,x.l2,x.rate); std::printf("  total %.5e 1/s\n",t);
  const double e1=ebkE1Rate(3,2,2,1); std::printf("3d: A_E1(EBK) %.5e, E2/E1 %.4e (QM 9.164e-6)\n",e1,t/e1);
  const auto s5=sStateQuadrupoleChannelRates(5.0); std::printf("5s E2 channels by s:"); for(double x: s5) std::printf(" %.4e",x); std::printf("\n");
  const auto s4=sStateQuadrupoleChannelRates(4.0); std::printf("4s E2 channels by s:"); for(double x: s4) std::printf(" %.4e",x); std::printf("\n");
  // M1
  const double h=2*pi*hbar;
  const double a21=spinWaveM1Rate(h*1420.405751768e6);
  std::printf("M1 hydrogen 21 cm: %.5e 1/s (735 QM 2.8755e-15) rel %+.2e\n",a21,a21/2.8755e-15-1);
  applyPairFromOption("electron,positron");
  const double aPs=spinWaveM1Rate(h*203.39e9);
  std::printf("M1 Ps 1^3S1 -> 1^1S0 at 203.39 GHz: %.5e 1/s (735 QM 3.377e-8) rel %+.2e\n",aPs,aPs/3.377e-8-1);
  std::printf("Ps E2 rate (kappa_2 = 0): %.3e\n",ebkE2Rate(3,2,1,0));
  const bool okM=std::abs(a21/2.8755e-15-1)<1e-3&&std::abs(aPs/3.377e-8-1)<1e-3;
  std::printf("H736a: R2 %d, G %d, A_E2 %d, M1 %d -> %s\n",okR,okG,okA,okM,(okR&&okG&&okA&&okM)?"TAK":"NIE");
}
