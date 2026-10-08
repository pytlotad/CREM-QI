#include "modules/crem_collapse.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 388: cost of the s-state spin transport at Ps 3s (a = 9 a_pair, L = hbar/2); argv[1] = elapsed [s], argv[2] = seed.
static Vec3 iso(std::mt19937_64& r){ std::uniform_real_distribution<double> u(0,1); const double c=2*u(r)-1,p=2*pi*u(r),s=std::sqrt(1-c*c); return {s*std::cos(p),s*std::sin(p),c}; }
int main(int argc,char**argv){
  const double T=argc>1?std::atof(argv[1]):1e-10; std::mt19937_64 r(argc>2?std::atoi(argv[2]):388);
  const double mred=firstMass*secondMass/(firstMass+secondMass), a=9.0*pairBohrRadius(activePair);
  const Vec3 m1=iso(r)*firstMagneticMoment, m2=iso(r)*secondMagneticMoment;
  const SecularSpinOrbitState st{Vec3{0,0,0.5*hbar},m1,m2,0.0,Vec3{1,0,0}};
  const auto t0=std::chrono::steady_clock::now();
  const auto one=orbitAveragedBmtAngularVelocities(a,st.orbitalAngularMomentum,m1,m2,mred,0.0,st.periapsisDirection);
  const auto t1=std::chrono::steady_clock::now();
  const auto warm=advanceCoupledSecularSpinOrbit(st,a,mred,T,0.05,1<<26);
  const auto t1b=std::chrono::steady_clock::now();
  const auto x=advanceCoupledSecularSpinOrbit(st,a,mred,T,0.05,1<<26);
  const auto t2=std::chrono::steady_clock::now();
  std::printf("first advance (with any one-time setup) %.3e s, %d substeps\n",std::chrono::duration<double>(t1b-t1).count(),warm.substeps);
  const double c1=std::chrono::duration<double>(t1-t0).count(), c2=std::chrono::duration<double>(t2-t1b).count();
  std::printf("nodes %d  one average %.3e s  advance T=%.3e s: completed %d substeps %d  (%.3e s total, %.3e s/substep)\n",one.phaseNodes,c1,T,x.completed,x.substeps,c2,c2/std::max(1,x.substeps));
  std::printf("m1 %.15e %.15e %.15e\nm2 %.15e %.15e %.15e\nL %.15e %.15e %.15e\nperi %.15e %.15e %.15e\n",
    x.state.firstDipole.x,x.state.firstDipole.y,x.state.firstDipole.z,x.state.secondDipole.x,x.state.secondDipole.y,x.state.secondDipole.z,
    x.state.orbitalAngularMomentum.x,x.state.orbitalAngularMomentum.y,x.state.orbitalAngularMomentum.z,
    x.state.periapsisDirection.x,x.state.periapsisDirection.y,x.state.periapsisDirection.z);
}
