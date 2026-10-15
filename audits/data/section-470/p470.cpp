// Audit 470: time-averaged mu_z of a hydrogen 1s atom in B_z on the model's secular path (Stern-Gerlach observable).
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
int main(int argc,char**argv){
  applyPairFromOption("proton,electron");
  const double B=std::atof(argv[1]); const int i0=std::atoi(argv[2]), i1=std::atoi(argv[3]);
  gExternalMagneticField=Vec3{0,0,B};
  const double m1=firstMass,m2=secondMass,M=m1+m2,mu=m1*m2/M,aB=pairBohrRadius(activePair);
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, M1=firstMagneticMoment, M2=secondMagneticMoment, muB=9.2740100783e-24;
  // orbital moment of the relative orbit: particle i carries L_i = L (m_j/M); mu_orb = sum q_i L_i/(2 m_i)
  const double orbFactor=firstCharge*(m2/M)/(2*m1)+secondCharge*(m1/M)/(2*m2);
  const double T=200e-12; const int K=200;
  for(int i=i0;i<i1;++i){
    std::mt19937_64 rng(470000+i); std::normal_distribution<double> g(0,1);
    auto unit=[&]{ Vec3 v{g(rng),g(rng),g(rng)}; return v*(1.0/v.norm()); };
    const Vec3 Lh=unit(), u1=unit(), u2=unit();
    Vec3 per=cross(Lh,Vec3{1,0,0}); if(per.norm()<1e-6) per=cross(Lh,Vec3{0,1,0}); per=per*(1.0/per.norm());
    SecularSpinOrbitState st; st.orbitalAngularMomentum=Lh*(0.5*hbar); st.firstDipole=u1*(s1*M1); st.secondDipole=u2*(s2*M2); st.periapsisDirection=per;
    double acc=0, accSpin=0, accOrb=0; bool ok=true;
    for(int k=0;k<K;++k){
      const auto adv=advanceCoupledSecularSpinOrbit(st,aB,mu,T/K); if(!adv.completed){ ok=false; break; } st=adv.state;
      const double ms=(st.firstDipole+st.secondDipole).z, mo=orbFactor*st.orbitalAngularMomentum.z; acc+=ms+mo; accSpin+=ms; accOrb+=mo; }
    std::printf("%d %d %.6f %.6f %.6f %.6f %.6f\n",i,int(ok),acc/K/muB,accSpin/K/muB,accOrb/K/muB,dot(u2,Lh),Lh.z);
    std::fflush(stdout);
  }
}
