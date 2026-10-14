// Audit 436 diagnostic: secular dE/dt = eps <cos theta dg/dt> on a FROZEN Kepler orbit (hbar = mu = k = 1), theta' = v^2.
#include <cmath>
#include <cstdio>
int main(){ const double e=0.6, r0=2.0;
  for(double J=1.80; J<=2.2001; J+=0.05){ const double a=J*J, w=1/(J*J*J), T=2*M_PI/w;
    double amp=0; for(double psi=0; psi<2*M_PI-1e-9; psi+=M_PI/4){
      const int orbits=200, N=4000; const double dt=T/N; double th=psi, s=0, gPrev=std::exp(-a*(1-e)/r0);
      for(int i=1;i<=orbits*N;++i){ const double t=i*dt, M=w*t; double E=M; for(int k=0;k<10;++k) E-=(E-e*std::sin(E)-M)/(1-e*std::cos(E));
        const double r=a*(1-e*std::cos(E)), v2=2/r-1/a, g=std::exp(-r/r0); th+=v2*dt; s+=std::cos(th)*(g-gPrev); gPrev=g; }
      const double Edot=s/(orbits*T); if(std::abs(Edot)>amp) amp=std::abs(Edot); }
    std::printf("J %.2f: max_psi |<cos theta dg/dt>| = %.3e  (x J^3 = dJ/dt per eps: %.3e)\n",J,amp,amp*J*J*J); }
}
