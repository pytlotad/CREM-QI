#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 514: Sommerfeld rule for Ps -- E(2^3P_J) = F(|L + lambda S|) + b <S12>_J, lambda = (2g - 1)/4; run with CREM_NO_S_CONTACT_ENERGY=1.
int main(){
  applyPairFromOption("electron,positron");
  const double aB=pairBohrRadius(activePair), a=4*aB, hP=2*pi*hbar, GHz=1e-9/hP, g=std::abs(firstGFactor), lam=(2*g-1)/4;
  const int l=1; const double r3=1.0/(aB*aB*aB*8*l*(l+0.5)*(l+1.0)); const double b=1e-7*firstMagneticMoment*secondMagneticMoment*r3;
  const double S12[3]={-4.0,2.0,-0.4}, ES1=7.4220, meas[3]={18.49965,13.01242,8.62438}, measS1P1=11.180;
  std::printf("Ps n = 2: g %.8f, lambda (2g-1)/4 = %.6f, b = %.5f GHz\n",g,lam,b*GHz);
  for(int mode=0;mode<3;++mode){ // 0 rule Q lengths, 1 QM lengths, 2 rule Q lengths with lambda = 1 (hydrogen form)
    const double lm=mode==2?1.0:lam; double E[3];
    const auto Lq=[&](double ll){ return mode==1?std::sqrt(ll*(ll+1)):ll+0.5; };
    for(int J=0;J<3;++J){ const double L=Lq(1), S=Lq(1), Jl=Lq(J); const double LS=(Jl*Jl-L*L-S*S)/2;
      const double Jeff=std::sqrt(std::max(0.0,L*L+lm*lm*S*S+2*lm*LS)); E[J]=(levelFineStructureEnergy(a,Jeff*hbar)+b*S12[J])*GHz; }
    const double L=Lq(1), S0=mode==1?0.0:0.5, J1=Lq(1); const double LS0=(J1*J1-L*L-S0*S0)/2;
    const double Jeff0=std::sqrt(std::max(0.0,L*L+lm*lm*S0*S0+2*lm*LS0)); const double E1P1=levelFineStructureEnergy(a,Jeff0*hbar)*GHz;
    std::printf("%s: P1-P0 %.4f (meas 5.4872) P2-P1 %.4f (meas 4.3880); residuals nu_J",mode==0?"rule Q lengths":(mode==1?"QM lengths    ":"rule Q, lambda 1"),E[1]-E[0],E[2]-E[1]);
    for(int J=0;J<3;++J) std::printf(" %+.4f",ES1-E[J]-meas[J]); std::printf("; 2^3S1-2^1P1 residual %+.4f GHz\n",ES1-E1P1-measS1P1);
  }
}
