#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 500: how many of each motion at n = 2, L = 1.5 hbar during a dwell time tau (argv[2], s); argv[1] = pair.
// Orbits (mean motion), periapsis (apsidal rate, integrated by the transport), spin precession rates, L and spins about J, relative azimuth
// of the two spins about J (rotation vs libration) -- from the engine's own secular transport over tau.
static double azimuth(const Vec3& v,const Vec3& J,const Vec3& ref){ const Vec3 p=v-J*dot(v,J); const Vec3 r2=cross(J,ref); return std::atan2(dot(p,r2),dot(p,ref)); }
int main(int argc,char**argv){
  applyPairFromOption(argv[1]); const double tau=std::atof(argv[2]);
  const double aB=pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass), a=4*aB, Lm=1.5*hbar;
  const double nmean=std::sqrt(pairCoulombStrength/(mu*a*a*a)), ecc=std::sqrt(1-std::pow(Lm/(2*hbar),2));
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  std::printf("pair %s: a = 4 a_pair, L = 1.5 hbar, e = %.4f, tau = %.4e s; orbital period %.4e s -> orbits %.4e\n",argv[1],ecc,tau,2*pi/nmean,tau*nmean/(2*pi));
  const Vec3 Lhat{0,0,1}, x{1,0,0};
  // spin configurations: directions of spin 1 and spin 2 (unit), labels
  struct Cfg{const char* name; Vec3 s1,s2;};
  const double c=std::cos(pi/3), s=std::sin(pi/3);
  std::vector<Cfg> cf={{"parallel, along L",{0,0,1},{0,0,1}},{"parallel, 60 deg to L",{s,0,c},{s,0,c}},{"parallel, 90 deg to L",{1,0,0},{1,0,0}},
                       {"parallel, against L",{0,0,-1},{0,0,-1}},{"antiparallel, 1 along L",{0,0,1},{0,0,-1}},{"antiparallel, 90 deg",{1,0,0},{-1,0,0}}};
  for(const Cfg& k: cf){
    SecularSpinOrbitState st{}; st.orbitalAngularMomentum=Lhat*Lm; st.periapsisDirection=x;
    st.firstDipole=k.s1*(g1*0.5*hbar); st.secondDipole=k.s2*(g2*0.5*hbar);      // mu = gamma S, |S| = hbar/2
    const auto rates=secularRatesWithSStateRule(a,st.orbitalAngularMomentum,st.firstDipole,st.secondDipole,mu,0.0,st.periapsisDirection);
    const auto aps=orbitAveragedApsidalRate(a,st.orbitalAngularMomentum,st.firstDipole,st.secondDipole,mu,st.periapsisDirection);
    Vec3 J=st.orbitalAngularMomentum+st.firstDipole/g1+st.secondDipole/g2; const double Jn=J.norm(); J=J*(1/Jn);
    Vec3 ref=std::abs(J.z)<0.9?cross(J,Vec3{0,0,1}):cross(J,Vec3{1,0,0}); ref=ref*(1/ref.norm());
    const int N=4000; double accL=0,acc1=0,acc2=0,accRel=0,aps2=0; double relMin=1e9,relMax=-1e9;
    double pL=azimuth(st.orbitalAngularMomentum,J,ref), p1=azimuth(st.firstDipole/g1,J,ref), p2=azimuth(st.secondDipole/g2,J,ref);
    bool ok=true; double relUnw=0;
    for(int i=0;i<N;++i){ const SecularSpinOrbitAdvance ad=advanceCoupledSecularSpinOrbit(st,a,mu,tau/N);
      if(!ad.completed){ ok=false; break; } st=ad.state; aps2+=ad.apsidalAngle;
      auto wrap=[](double d){ while(d>pi) d-=2*pi; while(d<-pi) d+=2*pi; return d; };
      const double qL=azimuth(st.orbitalAngularMomentum,J,ref), q1=azimuth(st.firstDipole/g1,J,ref), q2=azimuth(st.secondDipole/g2,J,ref);
      const double dL=wrap(qL-pL), d1=wrap(q1-p1), d2=wrap(q2-p2);
      accL+=dL; acc1+=d1; acc2+=d2; relUnw+=d1-d2; relMin=std::min(relMin,relUnw); relMax=std::max(relMax,relUnw); pL=qL; p1=q1; p2=q2; }
    std::printf("  [%s] |J|/hbar %.3f: spin rates |w1| %.4e |w2| %.4e rad/s; apsidal rate %.4e rad/s (rel %.3e, mag %.3e)\n",k.name,Jn/hbar,
      rates.first.norm(),rates.second.norm(),aps.total(),aps.relativistic,aps.magnetic);
    std::printf("      over tau (transport %s): periapsis turns %.3f; L about J %.3f turns; spin 1 about J %.3f; spin 2 about J %.3f; relative azimuth span %.3f turns (%s)\n",
      ok?"ok":"FAILED",aps2/(2*pi),accL/(2*pi),acc1/(2*pi),acc2/(2*pi),(relMax-relMin)/(2*pi),(relMax-relMin)>2*pi?"rotation":"libration");
    std::printf("      naive rate x tau: spin 1 %.3f, spin 2 %.3f, periapsis %.3f turns\n",rates.first.norm()*tau/(2*pi),rates.second.norm()*tau/(2*pi),aps.total()*tau/(2*pi));
  }
}
