#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 428: Ps n = 1, 2 S levels from the engine (after 427) + n = 2 P levels from 363/364; contact C_n test; 1S-2S.
int main(){
  gSingletChannelEligible=true;
  const double hP=2*pi*hbar, m=firstMass, k=pairCoulombStrength, mu=firstMagneticMoment, GHz=1e-9/hP;
  const double alpha=k/(hbar*c), a4=std::pow(alpha,4)*m*c*c, aPs=pairBohrRadius(activePair);
  auto Erel=[&](int n,double Lh){ return a4*(11.0/(64.0*n*n*n*n)-1.0/(4.0*n*n*n*Lh)); };
  const Vec3 ax{0,0,1}, P{1,0,0}, z{0,0,1};
  double US[3],UT[3];
  for(int n: {1,2}){ const double a=aPs*n*n, e=std::sqrt(1.0-0.25/(n*n));
    US[n]=spinCouplingEnergy(a,e,ax,P,z*(-mu),z*(-mu),0.5*hbar);   // e- moment -s1, e+ moment +s2; s1 = -s2: moments parallel
    UT[n]=spinCouplingEnergy(a,e,ax,P,z*(-mu),z*(mu),0.5*hbar);
    std::printf("n=%d: E_rel(S) %.4f  U_singlet %.4f  U_triplet %.4f  C_n %.4f (Darwin %.4f) GHz\n",n,Erel(n,0.5)*GHz,US[n]*GHz,UT[n]*GHz,a4/(4.0*n*n*n)*GHz,a4/(8.0*n*n*n)*GHz); }
  // n = 2 P levels (fs363.txt / singlet364.txt, reproduced byte for byte on the present code)
  const double EP[3]={-7.5976,-4.6187,-0.2868}, E1P1=-3.5357, meas[3]={18.49965,13.01242,8.62438}, measSP=11.180;
  for(int variant=0;variant<3;++variant){
    const double add=variant==0?0.0:(variant==1?a4/64.0:a4/32.0);
    const double ES1=(Erel(2,0.5)+UT[2]+add)*GHz;
    std::printf("%s: E(2^3S1) = %.4f GHz\n",variant==0?"engine after 427":(variant==1?"+ Darwin only":"+ full contact C_2"),ES1);
    for(int J=0;J<=2;++J) std::printf("   nu_%d model %+.3f measured %.3f residual %+.3f GHz\n",J,ES1-EP[J],meas[J],ES1-EP[J]-meas[J]);
    std::printf("   2^3S1-2^1P1 model %+.3f measured %.3f residual %+.3f GHz\n",ES1-E1P1,measSP,ES1-E1P1-measSP);
  }
  // 1S-2S triplet
  const double B1=pairBindingEnergy(activePair), lead=0.75*B1/hP*1e-6;   // MHz
  const double measured=1233607216.4;
  for(int withC=0;withC<=1;++withC){
    const double E41=Erel(1,0.5)+UT[1]+(withC?a4/4.0:0.0), E42=Erel(2,0.5)+UT[2]+(withC?a4/32.0:0.0);
    const double model=lead+(E42-E41)/hP*1e-6;
    std::printf("1S-2S %s: leading %.1f MHz, alpha^4 %+.1f MHz, model %.1f MHz, measured %.1f, residual %+.3f GHz\n",withC?"with C_n":"without C_n",lead,(E42-E41)/hP*1e-6,model,measured,(model-measured)/1000);
  }
}
