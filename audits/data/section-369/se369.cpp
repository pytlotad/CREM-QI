#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 369: engine S(e) (DFT harmonic table, interpolated in e^2) vs the classical ratio f(e) = omega <dL/dt> / <P>
// = (1-e^2)^(3/2) / (1 + e^2/2) for Larmor radiation of a Kepler orbit (<P> ~ (1+e^2/2)(1-e^2)^(-5/2), <dL/dt> ~ (1-e^2)^(-1)).
int main(){ double worst=0;
  for(double e: {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.866,0.9,0.95,0.98}){
    const double S=eccentricOrbitHazardSuppression(e), f=std::pow(1-e*e,1.5)/(1+e*e/2);
    worst=std::max(worst,std::abs(S/f-1)); std::printf("e=%.3f  S=%.6f  f=%.6f  S/f=%.5f\n",e,S,f,S/f); }
  std::printf("max |S/f - 1| = %.4f\n",worst); }
