// Audit 469: what the s-state isotropy rule does to the spin rates (transport) and what the M1 power reads.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  for(const char* pair: {"electron,positron","proton,electron"}){
    applyPairFromOption(pair);
    const double m1=firstMass,m2=secondMass,mu=m1*m2/(m1+m2),aB=pairBohrRadius(activePair);
    const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, M1=firstMagneticMoment, M2=secondMagneticMoment;
    const Vec3 z{0,0,1}, x{1,0,0}; const double th=M_PI/3; const Vec3 sh{std::sin(th),0,std::cos(th)};
    const double ph=70.5*M_PI/180; const Vec3 sp{std::sin(th+ph),0,std::cos(th+ph)};
    for(int para=0;para<3;++para){
      // spins parallel (0) / antiparallel (1) / at 70.5 deg to each other (2)
      const Vec3 mm1=sh*(s1*M1), mm2=(para==2?sp:sh)*((para==1?-1:1)*s2*M2);
      const Vec3 Lv=z*(0.5*hbar);
      const auto raw=orbitAveragedBmtAngularVelocities(aB,Lv,mm1,mm2,mu,0.0,x);
      const auto rule=secularRatesWithSStateRule(aB,Lv,mm1,mm2,mu,0.0,x);
      const auto orb=orbitAveragedBmtAngularVelocities(aB,Lv,mm1*1e-12,mm2*1e-12,mu,0.0,x);
      const auto em=coherentMagneticDipoleOrbitAveragedEmission(aB,Lv,mm1,mm2,mu,0.0,x);
      // independent contact-flow M1 (same formula written out again for H469a)
      const double f=(2.0/3.0)*mu0*quiggaRosnerContactDensity(aB,0.5*hbar), c1=-firstGyromagneticRatioOf()*f, c2=-secondGyromagneticRatioOf()*f;
      const Vec3 d1=cross(mm2,mm1)*c1, d2=cross(mm1,mm2)*c2, dd1=(cross(d2,mm1)+cross(mm2,d1))*c1, dd2=(cross(d1,mm2)+cross(mm1,d2))*c2;
      const double Pind=mu0*(dd1+dd2).squaredNorm()/(6*pi*c*c*c), w21=2*pi*1420.405751768e6, muB=9.2740100783e-24;
      const double Aqm=mu0*w21*w21*w21*muB*muB/(3*pi*hbar*c*c*c);
      std::printf("   contact-flow M1 independent %.6e W (engine/independent %.9f); contact precession |c1 m2| %.4e rad/s; H: rate/A_QM %.4f\n",
        Pind,Pind>0?em.power/Pind:0.0,(cross(mm2,mm1)*c1).norm()/mm1.norm(),em.power/(hbar*w21)/Aqm);
      std::printf("%s %s 1s: |w1| raw %.4e rule %.4e orbital-only %.4e; |w2| raw %.4e rule %.4e orbital-only %.4e; M1 power %.4e W\n",
        pair,para==0?"spins parallel":(para==1?"spins antiparallel":"spins at 70.5 deg"),raw.first.norm(),rule.first.norm(),orb.first.norm(),raw.second.norm(),rule.second.norm(),orb.second.norm(),em.power);
    }
  }
}
