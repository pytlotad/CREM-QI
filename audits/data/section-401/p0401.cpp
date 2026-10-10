#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <string>
// Audit 401: bound pair on a Kepler orbit of level n and eccentricity e0, started at apoapsis; mechanical path,
// stochastic photons (hazard scaled by CREM_MECHANICAL_HAZARD_SCALE); per-photon and per-window lines from
// CREM_DEBUG_LADDER / CREM_DEBUG_LADDER_SPLIT on stderr.
// argv: orbits (of the start orbit) [40], n [3], e0 [0], pair [default, e.g. proton,electron; "-" = default],
//       reaction model [stochastic | disabled], "nodipole" (magnetic moments set to zero; audit 401 diagnosis).
int main(int argc,char**argv){
  const double orbits=argc>1?std::atof(argv[1]):40.0;
  const double level=argc>2?std::atof(argv[2]):3.0;
  const double e0=argc>3?std::atof(argv[3]):0.0;
  if(argc>4&&std::string(argv[4])!="-") applyPairFromOption(argv[4]);
  const bool disabled=argc>5&&std::string(argv[5])=="disabled";
  const bool noDipole=argc>6&&std::string(argv[6])=="nodipole";
  if(noDipole){ firstMagneticMoment=0.0; secondMagneticMoment=0.0; }
  const double mu=pairReducedMass, a=level*level*pairBohrRadius(activePair);
  const double ra=a*(1.0+e0), va=std::sqrt(pairCoulombStrength/(mu*a)*(1.0-e0)/(1.0+e0));
  const double m1=firstMass, m2=secondMass;
  State s; s.firstPosition={ra*m2/(m1+m2),0,0}; s.secondPosition={-ra*m1/(m1+m2),0,0};
  s.firstVelocity={0,va*m2/(m1+m2),0}; s.secondVelocity={0,-va*m1/(m1+m2),0};
  s.firstProperDipole=Vec3{0,0,1}*firstMagneticMoment; s.secondProperDipole=Vec3{0,0,-1}*secondMagneticMoment;
  const double T=2*pi*std::sqrt(mu*a*a*a/pairCoulombStrength);
  SimulationOptions o; o.collectFrames=true; o.frameCount=200;
  const auto report=[](const char* tag,const State& x){
    const double E=conservativeParticleEnergy(x), n=std::sqrt(pairBindingEnergy(activePair)/-E);
    const double l=relativeOrbitalAngularMomentum(x)/hbar;
    const double nc=std::sqrt(pairBindingEnergy(activePair)/-conservedParticleEnergy(x));
    std::printf("%s n = %.6f, L/hbar = %.6f, e = %.6f, J_r = n - L/hbar = %.6f; n(conserved, no charge-dipole) = %.6f\n",tag,n,l,
                std::sqrt(std::max(0.0,1.0-l*l/(n*n))),n-l,nc);
  };
  report("start",s);
  const MechanicalTrajectoryResult r=runMechanicalTrajectory(s,orbits*T,comptonBarrierRadius,o,
      disabled?ChargeRadiationReactionModel::disabled:ChargeRadiationReactionModel::stochasticElectricDipole);
  report("end  ",r.finalState);
  std::printf("outcome %d, elapsed %.3f start orbits\n",(int)r.outcome,r.elapsedTime/T);
}
