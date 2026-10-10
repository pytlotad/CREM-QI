#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 399: bound Ps, circular n = 3, mechanical path, stochastic photons (hazard scaled by CREM_MECHANICAL_HAZARD_SCALE);
// argv[2] = a/a_Ps (default 9, n = 3); the per-photon energies come from CREM_DEBUG_LADDER on stderr; argv[1] = orbits to integrate.
int main(int argc,char**argv){
  const double orbits=argc>1?std::atof(argv[1]):40.0;
  const double mu=pairReducedMass, a=(argc>2?std::atof(argv[2]):9.0)*pairBohrRadius(activePair);
  const double v=std::sqrt(pairCoulombStrength/(mu*a)), m1=firstMass, m2=secondMass;
  State s; s.firstPosition={a*m2/(m1+m2),0,0}; s.secondPosition={-a*m1/(m1+m2),0,0};
  s.firstVelocity={0,v*m2/(m1+m2),0}; s.secondVelocity={0,-v*m1/(m1+m2),0};
  s.firstProperDipole=Vec3{0,0,1}*firstMagneticMoment; s.secondProperDipole=Vec3{0,0,-1}*secondMagneticMoment;
  const double T=2*pi*std::sqrt(mu*a*a*a/pairCoulombStrength);
  SimulationOptions o; o.collectFrames=true; o.frameCount=200;
  const double n0=std::sqrt(pairBindingEnergy(activePair)/-conservativeParticleEnergy(s));
  const MechanicalTrajectoryResult r=runMechanicalTrajectory(s,orbits*T,comptonBarrierRadius,o,
      ChargeRadiationReactionModel::stochasticElectricDipole);
  const double n1=std::sqrt(pairBindingEnergy(activePair)/-conservativeParticleEnergy(r.finalState));
  std::printf("start n = %.6f, end n = %.6f, outcome %d, elapsed %.3f orbits\n",n0,n1,(int)r.outcome,r.elapsedTime/T);
  for(size_t i=0;i<r.frames.size();i+=10){ const auto& f=r.frames[i]; (void)f; }
}
