#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 503: |L| under the secular transport alone in a uniform field, random spin / field directions (20 configurations).
// argv: B [T], time [s], configurations.
int main(int argc,char**argv){
  applyPairFromOption("proton,electron");
  const double Bm=std::atof(argv[1]), T=std::atof(argv[2]); const int K=std::atoi(argv[3]);
  const double aB=pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass), a=4*aB;
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  std::mt19937_64 r(5031); std::normal_distribution<double> G(0,1);
  auto rnd=[&](){ Vec3 v{G(r),G(r),G(r)}; return v*(1/v.norm()); };
  for(int k=0;k<K;++k){ gExternalMagneticField=rnd()*Bm;
    SecularSpinOrbitState st{}; st.orbitalAngularMomentum=Vec3{0,0,1.5*hbar}; Vec3 p=rnd(); p=p-Vec3{0,0,p.z}; st.periapsisDirection=p*(1/p.norm());
    st.firstDipole=rnd()*(g1*0.5*hbar); st.secondDipole=rnd()*(g2*0.5*hbar);
    double Lmin=1e9,Lmax=-1e9; int sub=0; bool ok=true;
    for(int i=0;i<200;++i){ const auto ad=advanceCoupledSecularSpinOrbit(st,a,mu,T/200); if(!ad.completed){ ok=false; break; } st=ad.state; sub+=ad.substeps;
      const double L=st.orbitalAngularMomentum.norm()/hbar; Lmin=std::min(Lmin,L); Lmax=std::max(Lmax,L); }
    std::printf("cfg %2d: %s substeps %6d |L|/hbar min %.5f max %.5f final %.5f\n",k,ok?"ok":"FAILED",sub,Lmin,Lmax,st.orbitalAngularMomentum.norm()/hbar); }
}
