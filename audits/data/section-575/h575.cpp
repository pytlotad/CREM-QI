#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 575: 1s hyperfine splitting of the model's level energy (spinCouplingEnergy) at n = 1, L = hbar/2,
// triplet (s1 = s2 along z) minus singlet (s1 = -s2), with a channel; env CREM_LO_QED=1 for the leading order.
int main(){
  gSingletChannelEligible=true;
  const double a=pairBohrRadius(activePair), m=firstMagneticMoment;
  const Vec3 z{0,0,1}, x{1,0,0}, y{0,1,0};
  // s_i = sign(q_i) m_i/|m_i|
  const double sg1=firstCharge<0?-1:1, sg2=secondCharge<0?-1:1;
  const Vec3 mT1=z*(sg1*m), mT2=z*(sg2*m), mS1=z*(sg1*m), mS2=z*(-sg2*m);
  const double L=0.5*hbar, e=std::sqrt(1-0.25);
  const double ET=spinCouplingEnergy(a,e,z,x,mT1,mT2,L), ES=spinCouplingEnergy(a,e,z,x,mS1,mS2,L);
  std::printf("QED %s: E_T %.6e J  E_S %.6e J  HFS %.4f GHz  (measured 203.3942, Ishida et al. 2014)\n",
    qedHigherOrder()?"alpha^2 / one loop":"leading order",ET,ES,(ET-ES)/(2*pi*hbar)/1e9);
  std::printf("rate factors: p-Ps %.6f  o-Ps %.6f  3g/2g ratio %.6e (LO %.6e)\n",paraRateQedFactor(),orthoRateQedFactor(),threePhotonRateRatio(),orePowellSuppression());
}
