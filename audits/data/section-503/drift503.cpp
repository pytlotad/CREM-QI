#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 503: |L| under the secular transport alone (no radiation) in a uniform field -- H, n = 2, L = 1.5 hbar, spins along random
// directions, B along z.  argv: B [T], total time [s], chunks.
int main(int argc,char**argv){
  applyPairFromOption("proton,electron");
  const double B=std::atof(argv[1]), T=std::atof(argv[2]); const int N=std::atoi(argv[3]);
  gExternalMagneticField=Vec3{0,0,B};
  const double aB=pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass), a=4*aB;
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  SecularSpinOrbitState st{}; const Vec3 Lh{std::sin(1.0),0,std::cos(1.0)}; st.orbitalAngularMomentum=Lh*(1.5*hbar); st.periapsisDirection=Vec3{0,1,0};
  const Vec3 s1{0.3,0.8,0.52}, s2{-0.6,0.2,0.77};
  st.firstDipole=s1*(g1*0.5*hbar/s1.norm()); st.secondDipole=s2*(g2*0.5*hbar/s2.norm());
  double Lmin=1e9,Lmax=-1e9; int sub=0;
  for(int i=0;i<N;++i){ const auto ad=advanceCoupledSecularSpinOrbit(st,a,mu,T/N); if(!ad.completed){ std::printf("FAILED at chunk %d\n",i); return 1; }
    st=ad.state; sub+=ad.substeps; const double L=st.orbitalAngularMomentum.norm()/hbar; Lmin=std::min(Lmin,L); Lmax=std::max(Lmax,L); }
  std::printf("B %.1f T, t %.3e s, substeps %d: |L|/hbar final %.6f, min %.6f, max %.6f, L.z/hbar final %.6f\n",B,T,sub,
    st.orbitalAngularMomentum.norm()/hbar,Lmin,Lmax,st.orbitalAngularMomentum.z/hbar);
}
