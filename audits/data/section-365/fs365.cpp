#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 365: finite-size charge of the model's current loop (R = 2 mu/(e c)) in S states, n = 2 deficit and the 1S-2S line.
int main(){
  const double hP=2*pi*hbar, m=firstMass, k=pairCoulombStrength, mu=firstMagneticMoment, GHz=1e-9/hP;
  const double alpha=k/(hbar*c), a4=std::pow(alpha,4)*m*c*c, aPs=pairBohrRadius(activePair), psi0=1.0/(pi*aPs*aPs*aPs);
  const double e=std::sqrt(k*4*pi*epsilon0), R=2*mu/(e*c), B=pairBindingEnergy(activePair);
  std::printf("R = 2 mu/(e c) = %.4f fm (lambda-bar_C = %.4f fm); g/2 = %.6f\n",R*1e15,hbar/(m*c)*1e15,R/(hbar/(m*c)));
  auto fs=[&](int n){ return (2*pi/3)*k*(2*R*R)*psi0/(n*n*n); };
  auto fsA=[&](int n){ return a4*std::pow(R/(hbar/(m*c)),2)/(6.0*n*n*n); };
  auto Erel=[&](int n,double Lh){ return a4*(11.0/(64.0*n*n*n*n)-1.0/(4.0*n*n*n*Lh)); };
  auto cont=[&](int n){ return (2*mu0/3)*mu*mu*psi0/(n*n*n); };
  std::printf("dE_fs: 1S %.4f GHz (analytic %.4f), 2S %.4f GHz (analytic %.4f); quantum Darwin alpha^4mc^2/(8n^3) for reference: 1S %.4f, 2S %.4f\n",
    fs(1)*GHz,fsA(1)*GHz,fs(2)*GHz,fsA(2)*GHz,a4/8*GHz,a4/64*GHz);
  const double deficit364=22.1208;
  std::printf("2^3S1 deficit (364) %.3f GHz -> after the loop %.3f GHz (loop explains %.1f %%)\n",deficit364,deficit364-fs(2)*GHz,100*fs(2)*GHz/deficit364);
  // 1S-2S, triplet, level difference
  const double LO=0.75*B, dRel=Erel(2,0.5)-Erel(1,0.5), dCont=cont(2)-cont(1), dFs=fs(2)-fs(1);
  const double meas=1233607216.4e6*hP;  // J
  std::printf("1S-2S: LO (3/4)B = %.6f GHz; d E_rel = %+.3f; d contact = %+.3f; d loop = %+.3f GHz\n",LO*GHz,dRel*GHz,dCont*GHz,dFs*GHz);
  std::printf("  measured %.6f GHz\n",meas*GHz);
  std::printf("  model - measured: LO only %+.3f; + E_rel %+.3f; + contact %+.3f; + loop %+.3f GHz\n",(LO-meas)*GHz,(LO+dRel-meas)*GHz,(LO+dRel+dCont-meas)*GHz,(LO+dRel+dCont+dFs-meas)*GHz);
  std::printf("  in eV: measured %.6f, model (all) %.6f\n",meas/1.602176634e-19,(LO+dRel+dCont+dFs)/1.602176634e-19);
}
