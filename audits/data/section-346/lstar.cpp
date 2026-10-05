#include "modules/crem_collapse.hpp"
#include <cstdio>
// Orbit-averaged contact density at n = 1 versus L/hbar (e = sqrt(1 - L^2)), against |psi_1s(0)|^2.
int main(){
  const double a1=pairBohrRadius(activePair), psi0=1.0/(pi*a1*a1*a1);
  for(double L: {1e-4,3e-4,1e-3,2e-3,3e-3,5e-3,1e-2,2e-2,5e-2,0.1,0.2,0.5,0.9}){
    const double e=std::sqrt(1.0-L*L);
    std::printf("L/hbar=%.4g e=%.8f rp/eps=%.4e n/psi0=%.4e\n",L,e,a1*(1-e)/magneticDipoleRadius(),
      orbitAveragedContactDensity(a1,std::min(e,0.999999))/psi0);
  }
}
