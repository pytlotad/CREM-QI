// Audit 462: J populations that import 18 assigns (engine's spinCouplingEnergy), production-like spins, Ps 2p/3d/4f.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <map>
#include <random>
int main(){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true;
  const double aB=pairBohrRadius(activePair), s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, M1=firstMagneticMoment, M2=secondMagneticMoment;
  const Vec3 z{0,0,1}, x{1,0,0}; std::mt19937_64 rng(462); std::normal_distribution<double> g(0,1);
  const auto unit=[&]{ Vec3 v{g(rng),g(rng),g(rng)}; return v*(1.0/v.norm()); };
  struct St{int n,l;}; for(St st: {St{2,1},St{3,2},St{4,3}}){
    const double a=aB*st.n*st.n, L=(st.l+0.5)*hbar;
    const auto E=[&](Vec3 u1,Vec3 u2){ return spinCouplingEnergy(a,0.0,z,x,u1*(s1*M1),u2*(s2*M2),L); };
    const double Eup=E(z,z), Edn=E(z*-1.0,z*-1.0);           // J = l+1 (S along L), J = l-1 (S against L)
    std::map<long long,long> cnt; std::map<long long,double> val; const int M=200000;
    for(int i=0;i<M;){ Vec3 u1=unit(), u2=unit(); if(dot(u1,u2)<0) continue; ++i; const double e=E(u1,u2);
      const long long k=std::llround(e/std::abs(Eup)*1e6); ++cnt[k]; val[k]=e; }
    std::printf("%d%c  E(J=l+1) %.9e J  E(J=l-1) %.9e J  levels found %zu\n",st.n,"spdf"[st.l],Eup,Edn,cnt.size());
    for(auto& [k,c]: cnt){ const double e=val[k]; int J=std::abs(e-Eup)<1e-9*std::abs(Eup)?st.l+1:(std::abs(e-Edn)<1e-9*std::abs(Eup)?st.l-1:st.l);
      const double p=double(c)/M, qm=(2*J+1)/(3.0*(2*st.l+1)), sg=std::sqrt(qm*(1-qm)/M);
      std::printf("   J=%d  E %.9e J  share %.4f  QM %.4f  (diff %+.1f sigma)\n",J,e,p,qm,(p-qm)/sg); }
  }
}
