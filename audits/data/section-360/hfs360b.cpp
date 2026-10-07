#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
// Audit 360 follow-up: variant (a) -- exact isotropic average (trace over x, y, z) vs Monte Carlo, and the contact prediction.
int main(){
  const double a1=pairBohrRadius(activePair), mu=firstMagneticMoment, psi0=1.0/(pi*a1*a1*a1), h=2*pi*hbar;
  const double e=std::sqrt(1.0-0.25); const Vec3 Lh{0,0,1}, P{1,0,0};
  auto dE=[&](const Vec3& s){ return orbitAveragedDipoleEnergy(a1,e,Lh,P,s*mu,s*(-mu))-orbitAveragedDipoleEnergy(a1,e,Lh,P,s*mu,s*mu); };
  const Vec3 ax[3]={{1,0,0},{0,1,0},{0,0,1}};
  double exact=0; for(auto& a: ax) exact+=dE(a)/3.0;
  const double nP=orbitAveragedContactDensity(a1,e);
  std::printf("exact isotropic (trace) Delta E = %.6e GHz\n",exact/h*1e-9);
  std::printf("contact prediction (4 mu0/3) mu^2 <n_P>_orb = %.6e GHz  (<n_P> = %.4e |psi0|^2)\n",(4*mu0/3)*mu*mu*nP/h*1e-9,nP/psi0);
  std::printf("per-axis Delta E [GHz]: x %.4f  y %.4f  z %.4f\n",dE(ax[0])/h*1e-9,dE(ax[1])/h*1e-9,dE(ax[2])/h*1e-9);
  for(int M: {4000,16000,64000}){ std::mt19937_64 rng(360); std::normal_distribution<double> N(0,1); double sum=0,sq=0;
    for(int k=0;k<M;++k){ Vec3 s{N(rng),N(rng),N(rng)}; s=s/s.norm(); const double v=dE(s)/h*1e-9; sum+=v; sq+=v*v; }
    const double m=sum/M, sd=std::sqrt(sq/M-m*m); std::printf("Monte Carlo M=%6d: mean %.4f GHz +- %.4f (sd/sqrt M), sd %.2f GHz\n",M,m,sd/std::sqrt(M),sd); }
}
