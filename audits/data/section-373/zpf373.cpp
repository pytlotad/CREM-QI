#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 373: Welton jitter in the model's zero-point field and the S-state shifts it gives (positronium).
int main(){
  const double hP=2*pi*hbar, m=firstMass, k=pairCoulombStrength, q=std::sqrt(k*4*pi*epsilon0), GHz=1e-9/hP;
  const double alpha=k/(hbar*c), lbar=hbar/(m*c), a5=std::pow(alpha,5)*m*c*c, B=pairBindingEnergy(activePair);
  const double aPs=pairBohrRadius(activePair), psi0=1.0/(pi*aPs*aPs*aPs);
  auto modelVariance=[&](double lo,double hi,int N,double omega){ // free charge driven by the model's modes
    const ZeroPointField f=makeZeroPointField(lo,hi,N,1.0,373); const double A=f.amplitudeCoefficient*omega*omega; double v=0;
    for(double fac: f.frequencyFactor){ const double w=fac*omega; v+=q*q*A*A/(2*m*m*w*w*w*w); } return v; };
  auto welton=[&](double lo,double hi){ return 2*alpha/pi*lbar*lbar*std::log(hi/lo); };
  for(int n=1;n<=2;++n){
    const double wn=2*B/(hbar*n*n*n);   // orbital angular frequency of level n
    const double vb=modelVariance(0.3,3.0,20000,wn), wb=welton(0.3,3.0);
    const double hiF=m*c*c/(hbar*wn); const double vf=modelVariance(1.0,hiF,200000,wn), wf=welton(1.0,hiF);
    std::printf("n=%d: hbar omega_n = %.3f eV; model band [0.3,3]: <dr^2> modes %.4e / Welton %.4e (ratio %.4f); physical [1, %.0f]: %.4e / %.4e (ratio %.4f)\n",
      n,hbar*wn/1.602176634e-19,vb,wb,vb/wb,hiF,vf,wf,vf/wf);
    const double nc=psi0/(n*n*n);
    const double dEb=(2*pi/3)*k*4*wb*nc, dEf=(2*pi/3)*k*4*wf*nc;
    std::printf("      S shift: model band %.4f GHz, physical band %.4f GHz (analytic (2/3pi) a^5 mc^2 ln/n^3 = %.4f GHz)\n",dEb*GHz,dEf*GHz,(2/(3*pi))*a5*std::log(hiF)/(n*n*n)*GHz);
  }
}
