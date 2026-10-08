#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 385 unit test: advanceCoupledSecularSpinOrbit at n = 1, L = hbar/2, free spins, 100 ps; argv[1] = substep, argv[2] = B [T], optional argv[3] = single trial.
static Vec3 iso(std::mt19937_64& r){ std::uniform_real_distribution<double> u(0,1); const double c=2*u(r)-1,p=2*pi*u(r),s=std::sqrt(1-c*c); return {s*std::cos(p),s*std::sin(p),c}; }
int main(int argc,char**argv){
  const double sub=std::atof(argv[1]), B=std::atof(argv[2]);
  const double mred=firstMass*secondMass/(firstMass+secondMass), a=pairBohrRadius(activePair);
  const int t0=argc>3?std::atoi(argv[3]):0, t1=argc>3?t0+1:3;
  for(int t=t0;t<t1;++t){ std::mt19937_64 r(500+t);
    const Vec3 m1=iso(r)*firstMagneticMoment, m2=iso(r)*secondMagneticMoment; gExternalMagneticField=iso(r)*B;
    const SecularSpinOrbitState st{Vec3{0,0,0.5*hbar},m1,m2,0.0,Vec3{1,0,0}};
    const auto x=advanceCoupledSecularSpinOrbit(st,a,mred,1e-10,sub,1<<26);
    std::printf("trial %d completed %d substeps %d  m1 %.12e %.12e %.12e  m2 %.12e %.12e %.12e  L %.12e %.12e %.12e  cos %.15f\n",t,x.completed,x.substeps,
      x.state.firstDipole.x/firstMagneticMoment,x.state.firstDipole.y/firstMagneticMoment,x.state.firstDipole.z/firstMagneticMoment,
      x.state.secondDipole.x/secondMagneticMoment,x.state.secondDipole.y/secondMagneticMoment,x.state.secondDipole.z/secondMagneticMoment,
      x.state.orbitalAngularMomentum.x/hbar,x.state.orbitalAngularMomentum.y/hbar,x.state.orbitalAngularMomentum.z/hbar,
      dot(x.state.firstDipole,x.state.secondDipole)/(x.state.firstDipole.norm()*x.state.secondDipole.norm()));
    std::printf("        cos0 %.15f\n",dot(m1,m2)/(m1.norm()*m2.norm())); }
}
