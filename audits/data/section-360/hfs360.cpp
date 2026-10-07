#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
// Audit 360: 1s hyperfine splitting at n = 1, L = hbar/2 under the current rules.
// (a) engine orbit-averaged dipole energy (Plummer contact + tensor), averaged over spin orientations (isotropy);
// (b) contact term with the Quigg-Rosner density |psi_1s(0)|^2.
int main(){
  const double a1=pairBohrRadius(activePair), mu=firstMagneticMoment, psi0=1.0/(pi*a1*a1*a1), h=2*pi*hbar;
  const double L=0.5, e=std::sqrt(1.0-L*L);
  std::mt19937_64 rng(360); std::normal_distribution<double> N(0,1);
  const Vec3 Lh{0,0,1}, P{1,0,0};
  double sumO=0,sumP=0; const int M=4000;
  for(int k=0;k<M;++k){ Vec3 s{N(rng),N(rng),N(rng)}; s=s/s.norm();
    sumP+=orbitAveragedDipoleEnergy(a1,e,Lh,P,s*mu,s*mu);      // p-Ps: parallel moments
    sumO+=orbitAveragedDipoleEnergy(a1,e,Lh,P,s*mu,s*(-mu)); } // o-Ps: antiparallel moments
  const double dA=(sumO-sumP)/M;
  const double contactOrb=orbitAveragedContactDensity(a1,e)/psi0;
  std::printf("(a) engine <U>_o - <U>_p (isotropic average of %d orientations) = %.4e J = %.4e GHz ; orbit contact density = %.4e |psi0|^2\n",M,dA,dA/h*1e-9,contactOrb);
  const double dB=(4.0*mu0/3.0)*mu*mu*psi0;
  std::printf("(b) contact with |psi_1s(0)|^2: (4 mu0/3) mu^2 |psi0|^2 = %.4e J = %.4f GHz = %.4f x measured 203.3942 GHz\n",dB,dB/h*1e-9,dB/h*1e-9/203.3942);
  std::printf("    QM magnetic part (8 mu0/3) mu^2 |psi0|^2 = %.4f GHz; /(4/7) = %.4f GHz (LO total with annihilation)\n",2*dB/h*1e-9,2*dB/h*1e-9*7/4);
}
