// Audit 486: inventory of every motion frequency in the model against the window 1e13-1e16 rad/s.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  const double W0=1e13, W1=1e16; auto inw=[&](double w){ return (w>=W0&&w<=W1)?"IN ":"   "; };
  struct St{const char* pair; const char* name; int n; double Lh;};
  St st[]={{"electron,positron","Ps 1s",1,0.5},{"electron,positron","Ps 2p",2,1.5},{"electron,positron","Ps 3d",3,2.5},{"electron,positron","Ps n10 l9",10,9.5},
           {"electron,positron","Ps n20 l19",20,19.5},{"proton,electron","H 1s",1,0.5},{"proton,electron","H 2p",2,1.5}};
  for(const auto& s: st){
    applyPairFromOption(s.pair);
    const double m1=firstMass,m2=secondMass,mu=m1*m2/(m1+m2),k=pairCoulombStrength,a=pairBohrRadius(activePair)*s.n*s.n;
    const double L=s.Lh*hbar, Lc=std::sqrt(mu*k*a), e=std::sqrt(std::max(0.0,1-(L/Lc)*(L/Lc))), worb=std::sqrt(k/(mu*a*a*a));
    const double kpk=0.504*std::pow(std::max(1-e,1e-6),-1.5);
    const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1; const Vec3 z{0,0,1}, x{1,0,0};
    const double th=M_PI/3; const Vec3 sh{std::sin(th),0,std::cos(th)};
    const Vec3 mm1=sh*(s1*firstMagneticMoment), mm2=sh*(s2*secondMagneticMoment);
    const auto r=orbitAveragedBmtAngularVelocities(a,z*L,mm1,mm2,mu,0.0,x);
    const auto ap=orbitAveragedApsidalRate(a,z*L,mm1,mm2,mu,x);
    const double wc=((r.first+r.second)*0.5).norm(), wd=(r.first-r.second).norm(), wmax=std::max(r.first.norm(),r.second.norm());
    const double wap=std::abs(ap.total()), larmor=std::max(std::abs(firstGyromagneticRatioOf()),std::abs(secondGyromagneticRatioOf()))*1.0;
    std::printf("== %s (a %.3e m, e %.3f)\n",s.name,a,e);
    std::printf("   %s orbit omega_orb                %.3e  rad/s\n",inw(worb),worb);
    std::printf("   %s orbit harmonics up to peak k~%.0f  %.3e\n",inw(kpk*worb),kpk,kpk*worb);
    std::printf("   %s spin common precession (= plane L about J) %.3e\n",inw(wc),wc);
    std::printf("   %s spin differential (mutual nutation)  %.3e\n",inw(wd),wd);
    std::printf("   %s fastest single spin rate          %.3e\n",inw(wmax),wmax);
    std::printf("   %s periapsis precession (rel %.2e, mag %.2e) %.3e\n",inw(wap),ap.relativistic,ap.magnetic,wap);
    std::printf("   %s within-orbit wobble of spin/L at omega_orb, amplitude Omega/omega = %.2e rad\n",inw(worb),wmax/worb);
    std::printf("   %s light crossing c/a                %.3e;  ping-pong ~ pi c/a %.3e\n",inw(c/a),c/a,pi*c/a);
    std::printf("   %s Larmor in 1 T (electron)          %.3e\n",inw(larmor),larmor);
  }
  const double rstar=comptonBarrierRadius;
  for(double R: {2*rstar,3*rstar,5*rstar}) std::printf("== current loop / sphere / ball R = %.0f r*: omega = c/R = %.3e rad/s %s\n",R/rstar,c/R,inw(c/R));
  std::printf("== Zitterbewegung 2 m_e c^2/hbar = %.3e rad/s\n",2*9.1093837015e-31*c*c/hbar);
}
