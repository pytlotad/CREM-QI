#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 525: (B) the engine's numerical partner-dipole part of the BMT rates (H 2p) vs the analytic classical tensor; (C) proton spin precession on
// 2P1/2 averaged over the fine-structure precession of L, s about J, component along J, vs Delta E(F)/h of import 22 and the measurement.
// Run with and without CREM_RULE_DIPOLE_TENSOR=1.
int main(){
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), a=4*aB, L=1.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), e=std::sqrt(1-0.75*0.75);
  const Vec3 z{0,0,1}, x{1,0,0};
  std::printf("CREM_RULE_DIPOLE_TENSOR %s\n",std::getenv("CREM_RULE_DIPOLE_TENSOR")?"on":"off");
  { // (B)
    const double muk=mu*pairCoulombStrength, r3=muk*std::sqrt(muk)/(a*std::sqrt(a)*L*L*L);
    double worst=0;
    for(double th:{10.0,50.0,90.0,140.0}){ const double t=th*pi/180; const Vec3 d{std::sin(t),0,std::cos(t)};
      const Vec3 mp=d*firstMagneticMoment, me=Vec3{std::cos(t),std::sin(t),0}*(-secondMagneticMoment);
      const auto full=orbitAveragedBmtAngularVelocities(a,z*L,mp,me,mu,0.0,x), tiny=orbitAveragedBmtAngularVelocities(a,z*L,mp*1e-12,me*1e-12,mu,0.0,x);
      const double fermi=(2.0/3.0)*mu0*orbitAveragedContactDensity(a,e);
      const Vec3 num=full.first-tiny.first-me*(-firstGyromagneticRatioOf()*fermi);   // contact part removed as applyContactAboveHbar does
      const Vec3 ana=(z*dot(z,me)-me*(1.0/3.0))*(-1.5e-7*r3)*(-firstGyromagneticRatioOf());
      worst=std::max(worst,(num-ana).norm()/ana.norm()); }
    std::printf("(B) proton rate from the electron dipole: max |numerical - analytic|/|analytic| = %.2e\n",worst);
  }
  { // (C) J along z, |J| = 1 hbar from rule-Q lengths |L| = 1.5, |s| = 1: 1.5 sin aL = sin as, 1.5 cos aL + cos as = 1
    const double aL=std::acos((2.25+1.0-1.0)/(2*1.5*1.0)), as=std::acos((1.0+1.0-2.25)/(2*1.0*1.0));   // triangle |L| 1.5, |s| 1, |J| 1
    std::printf("(C) cone angles: L^ %.3f deg, s^ %.3f deg from J; cos(L^, s^) = %+.4f; |1.5 L^ + s^| = %.4f\n",aL*180/pi,as*180/pi,std::cos(aL+as),
      std::sqrt(std::pow(1.5*std::sin(aL)-std::sin(as),2)+std::pow(1.5*std::cos(aL)+std::cos(as),2)));
    const int N=72; Vec3 sum; Vec3 sumOrb;
    for(int k=0;k<N;++k){ const double ph=2*pi*k/N;
      const Vec3 Lh{std::sin(aL)*std::cos(ph),std::sin(aL)*std::sin(ph),std::cos(aL)}, sh{-std::sin(as)*std::cos(ph),-std::sin(as)*std::sin(ph),std::cos(as)};
      const Vec3 me=sh*(-secondMagneticMoment), mp=Vec3{1,0,0}*firstMagneticMoment;
      Vec3 per=cross(Lh,Vec3{0,1,0}); if(per.norm()<1e-6) per=cross(Lh,x); per=per*(1.0/per.norm());
      OrbitAveragedBmtAngularVelocities r=secularRatesWithSStateRule(a,Lh*L,mp,me,mu,0.0,per);
      applyContactAboveHbar(r,a,Lh*L,mp,me,mu); applyBornInverseCube(r,Lh*L); applyRuleDipoleTensor(r,a,Lh*L,mp,me);
      sum=sum+r.first*(1.0/N);
      OrbitAveragedBmtAngularVelocities o=secularRatesWithSStateRule(a,Lh*L,mp,me*1e-12,mu,0.0,per);   // electron moment off: orbital field only
      applyContactAboveHbar(o,a,Lh*L,mp,me*1e-12,mu); applyBornInverseCube(o,Lh*L); applyRuleDipoleTensor(o,a,Lh*L,mp,me*1e-12);
      sumOrb=sumOrb+o.first*(1.0/N); }
    { const double fo=std::abs(sumOrb.z)/(2*pi)*1e-6, ft=std::abs(sum.z)/(2*pi)*1e-6;
      // QM split (Born squares, g = 2): orbital L.J/J^2 = 4/3, spin (3/4 - S.J)/J^2 = 4/3 of the total 8/3 -> half each of 59.1687
      std::printf("(C) split: orbital (electron moment 1e-12) %.4f MHz, spin dipole %.4f MHz; QM halves 29.58 / 29.58 MHz\n",fo,ft-fo); }
    const double f=std::abs(sum.z)/(2*pi)*1e-6, fperp=std::sqrt(sum.x*sum.x+sum.y*sum.y)/(2*pi)*1e-6;
    std::printf("(C) proton precession averaged over the fine-structure cone: along J %.4f MHz (perp %.2e); import 22 Delta E(F)/h 59.1687 MHz; measured 59.22(14):"
                " ratio f/59.1687 = %.4f\n",f,fperp,f/59.1687);
  }
}
