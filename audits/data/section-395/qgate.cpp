#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 395: reproduce the validation's quantized-radiation-gating trajectory (stochastic, retarded) step by step.
int main(){
  const double r=pairBohrRadius(activePair), mu=pairReducedMass;
  const double m1=firstMass,m2=secondMass, fs=m2/(m1+m2), ss=m1/(m1+m2);
  const double v=std::sqrt(pairCoulombStrength/(mu*r));
  State s; s.firstPosition={fs*r,0,0}; s.secondPosition={-ss*r,0,0};
  s.firstVelocity={0,fs*v,0}; s.secondVelocity={0,-ss*v,0};
  s.firstProperDipole=Vec3{0.3,0.5,0.81}*(firstMagneticMoment/0.99929975);
  s.secondProperDipole=Vec3{-0.4,0.62,0.67}*(secondMagneticMoment/1.00344407);
  synchronizeCovariantDipoles(s);
  ClassicalTrajectoryEngine e(s,{.relativeTolerance=1.0e-8,.maximumDepth=12,
     .reactionModel=ChargeRadiationReactionModel::stochasticElectricDipole,.computeOutwardFlux=true,.useRetardedExternalForces=true});
  const double T=2*pi*std::sqrt(mu*r*r*r/pairCoulombStrength);
  const double E0=conservativeParticleEnergy(s);
  for(int k=0;k<64;++k){ const bool ok=e.advance(s,T/256.0);
    const MutualForces f=retardedExternalForces(s,e.history());
    const auto rad=particleMultipoleRadiation(s,f,e.history(),true,ChargeRadiationReactionModel::stochasticElectricDipole,true);
    std::printf("step %2d ok=%d finite=%d |r|/a=%.9f |F|=%.3e |react|=%.3e hist=%zu\n",k,ok,isFinite(s),separation(s)/r,
      f.first.norm(),rad.chargeReaction.first.norm(),e.history().size());
    if(!ok||!isFinite(s)) break; }
  const double larmor=std::pow(eCharge,2)*std::pow(pairCoulombStrength/(mu*r*r),2)/(6*pi*epsilon0*c*c*c);
  std::printf("dE_conservative over 64 steps (T/4) = %.6e J; coherent Larmor x T/4 = %.6e J; ratio %.4f\n",
    conservativeParticleEnergy(s)-E0,-larmor*T/4,(conservativeParticleEnergy(s)-E0)/(-larmor*T/4));
}
