// Audit 468 part A: the engine's M1 (spin-flip) emission rate for hydrogen 1s vs the classical formula and vs A_QM(21 cm).
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("proton,electron"); gSingletChannelEligible=false;
  const double m1=firstMass,m2=secondMass,mu=m1*m2/(m1+m2),k=pairCoulombStrength,aB=pairBohrRadius(activePair);
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, M1=firstMagneticMoment, M2=secondMagneticMoment;
  const double w21=2*pi*1420.405751768e6, muB=9.2740100783e-24;
  const double Aqm=mu0*w21*w21*w21*muB*muB/(3*pi*hbar*c*c*c);
  std::printf("# hydrogen 1s, a = %.4e m, L = hbar/2; |mu_p| %.4e |mu_e| %.4e J/T; A_QM(21 cm, formula) = %.4e 1/s\n",aB,M1,M2,Aqm);
  std::printf("# th12[deg]  |W_p| |W_e| [rad/s]  W_p.z W_e.z  P_M1 [W]  P_classical(avg rates) [W]  rate=P/(hbar w21)  rate/A_QM   U-flip(e) freq [Hz]\n");
  const Vec3 z{0,0,1}, x{1,0,0};
  for(double th: {10.,30.,60.,90.,120.,150.,170.}){
    const double t=th*pi/180; const Vec3 up{0,0,1}, u2{std::sin(t),0,std::cos(t)};  // proton spin along z, electron at th12
    const Vec3 mp=up*(s1*M1), me=u2*(s2*M2);
    const auto em=coherentMagneticDipoleOrbitAveragedEmission(aB,z*(0.5*hbar),mp,me,mu,0.0,x);
    const auto r=orbitAveragedBmtAngularVelocities(aB,z*(0.5*hbar),mp,me,mu,0.0,x);
    const Vec3 dd=cross(r.first,cross(r.first,mp))+cross(r.second,cross(r.second,me));
    const double Pcl=mu0*dd.squaredNorm()/(6*pi*c*c*c);
    // model's M1 photon energy for flipping the electron: U(before) - U(after) with me -> -me
    const double U0=spinCouplingEnergy(aB,std::sqrt(1-0.25*hbar*hbar/(mu*k*aB)),z,x,mp,me,0.5*hbar);
    const double U1=spinCouplingEnergy(aB,std::sqrt(1-0.25*hbar*hbar/(mu*k*aB)),z,x,mp,me*-1.0,0.5*hbar);
    std::printf("%6.1f  %.4e %.4e  %+.3e %+.3e  %.4e  %.4e  %.4e  %.4f  %.6e\n",th,r.first.norm(),r.second.norm(),r.first.z,r.second.z,
      em.power,Pcl,em.power/(hbar*w21),em.power/(hbar*w21)/Aqm,(U0-U1)/(2*pi*hbar));
  }
}
