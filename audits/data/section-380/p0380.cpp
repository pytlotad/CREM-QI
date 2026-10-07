#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <complex>
// Audit 380: internal clock of the moment's current loop (Omega = m c^2/hbar) on Kepler orbits of positronium.
int main(){
  const double mu=firstMass*secondMass/(firstMass+secondMass), K=pairCoulombStrength, a1=pairBohrRadius(activePair);
  const double m1=firstMass, m2=secondMass, Om1=m1*c*c/hbar, Om2=m2*c*c/hbar, lam=hbar/(m1*c);
  std::printf("Omega = %.6e rad/s, lambdabar_C = %.6e m, a_Ps/lambdabar = %.3f\n",Om1,lam,a1/lam);
  std::printf("loop R = 2 r* = %.6e m; c/R = %.6e rad/s (should equal Omega)\n",2*comptonBarrierRadius,c/(2*comptonBarrierRadius));
  std::printf("\n n  e      L/hbar   max|d(phi1-phi2)/dt|/Omega  dilation phase/orbit per clock /(pi J/2hbar)   <cos psi>   <sin psi>\n");
  for(int n=1;n<=4;++n) for(double e: {0.0,0.3,0.6,0.866,0.95}){
    const double a=n*n*a1, w=std::sqrt(K/(mu*a*a*a)), T=2*pi/w, L=std::sqrt(mu*K*a*(1-e*e));
    const int N=200000; double dil1=0,maxBeat=0; std::complex<double> avg=0;
    for(int i=0;i<N;++i){
      const double E=2*pi*(i+0.5)/N, tw=1-e*std::cos(E), dt=tw*(T/N);           // dt = (1 - e cos E) dE / w
      const double r=a*tw, v2=K/mu*(2/r-1/a);                                   // relative speed^2 (vis-viva)
      const double v1sq=v2*(m2/(m1+m2))*(m2/(m1+m2)), v2sq=v2*(m1/(m1+m2))*(m1/(m1+m2));
      const double rate1=Om1*(std::sqrt(1-v1sq/(c*c))), rate2=Om2*(std::sqrt(1-v2sq/(c*c)));
      maxBeat=std::max(maxBeat,std::abs(rate1-rate2)/Om1);
      dil1+=Om1*(1-std::sqrt(1-v1sq/(c*c)))*dt;
      avg+=std::polar(1.0,-Om1*r/c)*(dt/T);                                      // psi = -Omega r/c (+const)
    }
    const double J=std::sqrt(mu*K*a);
    std::printf("%2d %.3f  %7.4f   %.3e                    %.8f                         %+.4f     %+.4f\n",n,e,L/hbar,maxBeat,dil1/(pi*J/(2*hbar)),avg.real(),avg.imag());
  }
  // H380c scan: at n = 2, <cos psi> vs L/hbar on a fine grid -- look for structure at integer L
  std::printf("\nscan n = 2: L/hbar  |<e^{i psi}>|\n");
  for(double Lh=0.05;Lh<=2.0001;Lh+=0.05){
    const double a=4*a1, w=std::sqrt(K/(mu*a*a*a)), T=2*pi/w, e=std::sqrt(std::max(0.0,1-(Lh*hbar)*(Lh*hbar)/(mu*K*a)));
    const int N=400000; std::complex<double> avg=0;
    for(int i=0;i<N;++i){ const double E=2*pi*(i+0.5)/N, tw=1-e*std::cos(E); avg+=std::polar(1.0,-Om1*a*tw/c)*(tw/N); }
    std::printf("  %.2f  %.5f   (a e/lambdabar = %.1f)\n",Lh,std::abs(avg),a*e/lam);
  }
  // H380d: locking to a fixed-frequency background (per-orbit slip = dilation phase, both clocks) vs de Broglie wave control
  std::printf("\nlocking per orbit (circular): n  slip per clock/2pi = J/(4 hbar)   de Broglie wave phase/2pi = J/hbar\n");
  for(int n=1;n<=4;++n){ const double a=n*n*a1, J=std::sqrt(mu*K*a), v=std::sqrt(K/(mu*a)), T=2*pi*a/v;
    const double slip=Om1*(1-std::sqrt(1-(v/2)*(v/2)/(c*c)))*T, wave=mu*v*2*pi*a/hbar;   // relative orbit: k = mu v/hbar
    std::printf("   %d  %.6f   %.6f\n",n,slip/(2*pi),wave/(2*pi)); }
}
