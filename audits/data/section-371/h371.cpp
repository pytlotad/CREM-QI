#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 371: hydrogen 2p / 3p lifetime factors -- engine P_E1 vs analytic Kepler-Larmor, and the semiclassical tau for L choices.
int main(){
  applyPairFromOption("proton,electron");
  const double mred=firstMass*secondMass/(firstMass+secondMass), k=pairCoulombStrength, a0=pairBohrRadius(activePair);
  const double q=std::sqrt(k*4*pi*epsilon0);
  auto f=[](double e){ return std::pow(1-e*e,1.5)/(1+e*e/2); };
  auto Larmor=[&](double a,double e){ return q*q*k*k*(1+e*e/2)/(std::pow(a,4)*std::pow(1-e*e,2.5))/(6*pi*epsilon0*c*c*c*mred*mred); };
  struct S{const char* name; int n,l; double meas;} st[]={{"2p",2,1,1.600},{"3p",3,1,5.58}};
  for(auto& s: st){
    const double a=s.n*s.n*a0, w=std::sqrt(k/(mred*a*a*a));
    std::printf("%s (n=%d): measured %.3f ns\n",s.name,s.n,s.meas);
    for(int choice=0;choice<2;++choice){
      const double Lh=choice==0?s.l+0.5:std::sqrt(s.l*(s.l+1.0)); const double e=std::sqrt(std::max(0.0,1-std::pow(Lh/s.n,2)));
      const double PL=Larmor(a,e);
      const auto em=secularElectricDipoleOrbitAveragedEmission(a,Vec3{0,0,1}*(Lh*hbar),Vec3{1,0,0}*(firstMagneticMoment*1e-12),Vec3{0,1,0}*(secondMagneticMoment*1e-12),mred,Vec3{1,0,0});
      const double tauL=hbar*w/(f(e)*PL), tauE=em.valid?hbar*w/(f(e)*em.power):0;
      std::printf("  L=%-14s e=%.4f f=%.4f  P_L=%.4e W  P_engine/P_L=%.5f (force %.5f, spectral %.5f)  tau(Larmor)=%.4f ns (%+.1f%%)  tau(engine P)=%.4f ns (%+.1f%%)\n",
        choice==0?"(l+1/2)hbar":"sqrt(l(l+1))hbar",e,f(e),PL,em.power/PL,em.forcePower/PL,em.spectralFactor,tauL*1e9,100*(tauL*1e9/s.meas-1),tauE*1e9,100*(tauE*1e9/s.meas-1));
    }
  }
}
