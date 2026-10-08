#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
int main(){ std::mt19937_64 r(7); std::uniform_real_distribution<double> U(0,1);
  for(double e: {0.5528,0.6614,0.8660,0.9682,0.9860}){ int N=400000,c1=0; double mk=0;
    for(int i=0;i<N;++i){ const double k=eccentricOrbitHarmonicNumber(e,U(r)); if(k<=1.0) ++c1; mk+=k; }
    std::printf("e=%.4f  w1 sampled %.5f   <k> %.3f   S(e)<k> %.4f\n",e,c1/(double)N,mk/N,eccentricOrbitHazardSuppression(e)*mk/N); } }
