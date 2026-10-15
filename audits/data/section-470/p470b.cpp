// Audit 470b: projection of spin 1 on the direction of the PARTNER's field acting on it, and cos(theta12), over time.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
int main(int argc,char**argv){
  applyPairFromOption(argv[1]); const int nlev=std::atoi(argv[2]); const double Lh_hbar=std::atof(argv[3]);
  const int N=std::atoi(argv[4]); const double T=std::atof(argv[5]); const int K=100;
  const double m1=firstMass,m2=secondMass,mu=m1*m2/(m1+m2),aB=pairBohrRadius(activePair), a=aB*nlev*nlev;
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, M1=firstMagneticMoment, M2=secondMagneticMoment;
  gSingletChannelEligible=false;
  auto partnerDir=[&](const SecularSpinOrbitState& st){
    const auto full=secularRatesWithSStateRule(a,st.orbitalAngularMomentum,st.firstDipole,st.secondDipole,mu,0.0,st.periapsisDirection);
    const auto own=secularRatesWithSStateRule(a,st.orbitalAngularMomentum,st.firstDipole*1e-12,st.secondDipole*1e-12,mu,0.0,st.periapsisDirection);
    Vec3 w=full.first-own.first;   // = -gamma_1 B_partner
    const double gam=firstGyromagneticRatioOf(); Vec3 Bp=w*(-1.0/gam); return Bp.norm()>0?Bp*(1.0/Bp.norm()):Bp; };
  int hi0=0,hiT=0,hiAvg=0; double maxdc=0;
  for(int i=0;i<N;++i){
    std::mt19937_64 rng(4700000+i); std::normal_distribution<double> g(0,1);
    auto unit=[&]{ Vec3 v{g(rng),g(rng),g(rng)}; return v*(1.0/v.norm()); };
    const Vec3 Lh=unit(), u1=unit(), u2=unit();
    Vec3 per=cross(Lh,Vec3{1,0,0}); if(per.norm()<1e-6) per=cross(Lh,Vec3{0,1,0}); per=per*(1.0/per.norm());
    SecularSpinOrbitState st; st.orbitalAngularMomentum=Lh*(Lh_hbar*hbar); st.firstDipole=u1*(s1*M1); st.secondDipole=u2*(s2*M2); st.periapsisDirection=per;
    auto proj=[&](const SecularSpinOrbitState& s){ const Vec3 b=partnerDir(s); return dot(s.firstDipole*(1.0/s.firstDipole.norm()),b); };
    auto c12=[&](const SecularSpinOrbitState& s){ return dot(s.firstDipole,s.secondDipole)*s1*s2/(s.firstDipole.norm()*s.secondDipole.norm()); };
    const double p0=proj(st), c0=c12(st); double pav=0, cmin=c0, cmax=c0;
    for(int k=0;k<K;++k){ const auto adv=advanceCoupledSecularSpinOrbit(st,a,mu,T/K); st=adv.state; const double p=proj(st); pav+=p/K; const double c=c12(st); cmin=std::min(cmin,c); cmax=std::max(cmax,c); }
    const double pT=proj(st); hi0+=std::abs(p0)>0.9; hiT+=std::abs(pT)>0.9; hiAvg+=std::abs(pav)>0.9; maxdc=std::max(maxdc,cmax-cmin);
    if(i<5) std::printf("  sample %d: p0 %+.4f pT %+.4f <p> %+.4f  cos12 %+.4f range %.2e\n",i,p0,pT,pav,c0,cmax-cmin);
  }
  std::printf("%s n=%d L=%.1f hbar N=%d T=%.2e s: |p|>0.9 at t=0 %d (%.3f), at T %d (%.3f), time-avg %d (%.3f); isotropic 0.100; max range of cos12 over T %.3e\n",
    argv[1],nlev,Lh_hbar,N,T,hi0,double(hi0)/N,hiT,double(hiT)/N,hiAvg,double(hiAvg)/N,maxdc);
}
