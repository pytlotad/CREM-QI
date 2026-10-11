#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <string>
// Audit 406: Ps Kepler orbit (n, e0) from apoapsis with both moments tilted by theta from the orbit normal
// (antiparallel "ortho" or parallel "para"), mechanical path, no photons; run with CREM_DUMP_FAIL + CREM_DUMP_BELOW_R
// to dump a state for h0406.  argv: orbits n e0 theta_deg para|ortho.
int main(int argc,char**argv){
  const double orbits=std::atof(argv[1]), level=std::atof(argv[2]), e0=std::atof(argv[3]), th=std::atof(argv[4])*pi/180.0;
  const bool para=std::string(argv[5])=="para";
  const double mu=pairReducedMass, a=level*level*pairBohrRadius(activePair);
  const double ra=a*(1.0+e0), va=std::sqrt(pairCoulombStrength/(mu*a)*(1.0-e0)/(1.0+e0));
  const double m1=firstMass, m2=secondMass;
  State s; s.firstPosition={ra*m2/(m1+m2),0,0}; s.secondPosition={-ra*m1/(m1+m2),0,0};
  s.firstVelocity={0,va*m2/(m1+m2),0}; s.secondVelocity={0,-va*m1/(m1+m2),0};
  const Vec3 axis{std::sin(th),0,std::cos(th)};
  s.firstProperDipole=axis*firstMagneticMoment; s.secondProperDipole=axis*(para?secondMagneticMoment:-secondMagneticMoment);
  const double T=2*pi*std::sqrt(mu*a*a*a/pairCoulombStrength);
  SimulationOptions o;
  const MechanicalTrajectoryResult r=runMechanicalTrajectory(s,orbits*T,comptonBarrierRadius,o,
      ChargeRadiationReactionModel::stochasticElectricDipole);
  std::printf("outcome %d elapsed %.3f orbits\n",(int)r.outcome,r.elapsedTime/T);
}
