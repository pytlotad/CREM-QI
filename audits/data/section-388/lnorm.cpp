#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 388: does |L| stay bit-constant across s-state substeps (Ps 3s)?
int main(){
  const double mred=firstMass*secondMass/(firstMass+secondMass), a=9.0*pairBohrRadius(activePair);
  SecularSpinOrbitState st{Vec3{0,0,0.5*hbar},Vec3{0.3,0.4,0.866}*(firstMagneticMoment/1.0412),Vec3{-0.5,0.1,0.86}*(secondMagneticMoment/1.0005),0.0,Vec3{1,0,0}};
  for(int i=0;i<8;++i){ const auto x=advanceCoupledSecularSpinOrbit(st,a,mred,1e-12,0.05,1); st=x.state;
    std::printf("step %d substeps %d |L|/hbar = %.17e\n",i,x.substeps,st.orbitalAngularMomentum.norm()/hbar); }
}
