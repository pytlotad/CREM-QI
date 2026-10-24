// Audit 630: is a Kustaanheimo-Stiefel integrator needed?  MODE=nodes: production E1 loss quadrature (secularElectricDipoleOrbitAveragedEmission
// on the analytic Kepler orbit) against the closed form (larmorOrbitAveragedPower) for e -> 1.  MODE=drift: energy conservation of the mechanical
// path (conservative stochastic mode, no photons, moments zero) on Ps n = 2 orbits with e = ECC over PERIODS orbits.  env PAIR, MODE, ECC, PERIODS.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
static double envd(const char* k,double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
int main(){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  const std::string mode=std::getenv("MODE")?std::getenv("MODE"):"nodes";
  const double mu=pairReducedMass, k=pairCoulombStrength, a0=pairBohrRadius(activePair);
  if(mode=="nodes"){
    const double a=4*a0;
    for(double e: {0.0,0.5,0.9,0.97,0.99,0.997,0.999,0.9999}){
      const double L=std::sqrt((1-e*e)*mu*k*a);
      const SecularElectricDipoleEmission s=secularElectricDipoleOrbitAveragedEmission(a,Vec3{0,0,L},Vec3{0,0,0},Vec3{0,0,0},mu,Vec3{1,0,0});
      const double closed=larmorOrbitAveragedPower(a,e);
      std::printf("NODES e %.4f nodes %d kinematic/closed %.12f power/closed %.12f force/kinematic %.12f\n",e,s.phaseNodes,
                  s.kinematicPower/closed,s.power/closed,s.forcePower/s.kinematicPower);
    }
    return 0;
  }
  // drift
  firstMagneticMoment=0.0; secondMagneticMoment=0.0;
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
  const double e=envd("ECC",0.9), n=2;
  gInitialSeparationOverride=n*n*a0; gInitialAngularMomentumFraction=std::sqrt(1-e*e);
  const double a=n*n*a0, T=2*pi*std::sqrt(a*a*a*mu/k);
  SimulationOptions o; o.frameCount=2; o.observationTime=envd("PERIODS",5)*T; o.radiatedEnergyBookkeeping=false;
  double E0=0, maxDev=0, rmin=1e300; long steps=0; std::vector<std::pair<double,double>> perOrbit;
  int lastOrbit=-1;
  o.stepReady=[&](const State& s){
    State z=s; z.firstDipole={0,0,0}; z.secondDipole={0,0,0}; z.dipoleConstraintEnergy=0;
    const double E=conservativeParticleEnergy(z);
    if(steps==0) E0=E;
    maxDev=std::max(maxDev,std::abs(E-E0)/std::abs(E0)); rmin=std::min(rmin,(s.firstPosition-s.secondPosition).norm());
    const int orbit=static_cast<int>(s.time/T);
    if(orbit!=lastOrbit){ perOrbit.push_back({s.time,E}); lastOrbit=orbit; }
    ++steps; };
  simulate(splitMix64(630ULL),1,o);
  const double drift=perOrbit.size()>1?(perOrbit.back().second-perOrbit.front().second)/std::abs(E0)/((perOrbit.back().first-perOrbit.front().first)/T):0.0;
  std::printf("DRIFT e %.4f steps %ld orbits %.2f r_min/a %.3e max|dE/E| %.3e drift/orbit %.3e (dn/orbit %.3e)\n",e,steps,o.observationTime/T,rmin/a,maxDev,drift,0.5*n*drift);
}
