// Audit 467: engine spin precession (orbit-averaged Thomas-BMT, by sector) vs classical Breit-Pauli 1/c^2 terms.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  const char* sec=std::getenv("CREM_SPIN_RATE_SECTOR"); std::printf("# sector %s\n",sec?sec:"full");
  for(const char* pair: {"electron,positron","proton,electron"}){
    applyPairFromOption(pair);
    const double m1=firstMass,m2=secondMass,mu=m1*m2/(m1+m2),k=pairCoulombStrength,aB=pairBohrRadius(activePair);
    const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1, M1=firstMagneticMoment, M2=secondMagneticMoment;
    const double g1=std::abs(firstGFactor), g2=std::abs(secondGFactor);
    const double n=2, a=aB*n*n, L=1.5*hbar, Lc=std::sqrt(mu*k*a), e=std::sqrt(1-(L/Lc)*(L/Lc));
    const Vec3 z{0,0,1}, x{1,0,0}, y{0,1,0};
    const double r3=1/(a*a*a*std::pow(1-e*e,1.5));
    const double A1=(k/(c*c))*((g1-1)/(2*m1*m1)+g1/(2*m1*m2)), A2=(k/(c*c))*((g2-1)/(2*m2*m2)+g2/(2*m1*m2));
    std::printf("# %s: g %.6f %.6f, a %.4e m, e %.4f, <r^-3> %.4e; Breit SO: Omega_1 = %.6e, Omega_2 = %.6e rad/s along L\n",pair,g1,g2,a,e,r3,A1*r3*L,A2*r3*L);
    struct Cfg{const char* name; Vec3 u1,u2;};
    const double th=M_PI/3; const Vec3 sh{std::sin(th),0,std::cos(th)}, sh2{-std::sin(th),0,std::cos(th)};
    Cfg cfgs[]={{"along L",z,z},{"60deg parallel",sh,sh},{"perp (periapsis)",x,x},{"perp (y)",y,y},{"60deg antiparallel",sh,sh*-1.0},{"mixed 60/120",sh,sh2}};
    for(auto& cf: cfgs){
      const Vec3 mm1=cf.u1*(s1*M1), mm2=cf.u2*(s2*M2);
      const auto r=orbitAveragedBmtAngularVelocities(a,z*L,mm1,mm2,mu,0.0,x);
      // Breit-Pauli spin-spin, orbit average by quadrature in eccentric anomaly (time weight 1 - e cos E)
      Vec3 B1,B2; const int N=4096; double W=0;
      for(int i=0;i<N;++i){ const double E=2*M_PI*(i+0.5)/N, w=1-e*std::cos(E);
        const Vec3 rel=x*(a*(std::cos(E)-e))+y*(a*std::sqrt(1-e*e)*std::sin(E)); // r_1 - r_2
        const double rr=rel.norm(); const Vec3 nn=rel*(1/rr); const double c3=1e-7/(rr*rr*rr);
        B1+=(nn*(3*dot(nn,mm2))-mm2)*(c3*w); B2+=(nn*(3*dot(nn,mm1))-mm1)*(c3*w); W+=w; }
      B1=B1*(1/W); B2=B2*(1/W);
      const double gam1=g1*firstCharge/(2*m1), gam2=g2*secondCharge/(2*m2);
      const Vec3 SS1=B1*(-gam1), SS2=B2*(-gam2), SO1=z*(A1*r3*L), SO2=z*(A2*r3*L);
      auto pr=[](const char* t,const Vec3& v){ std::printf("    %-9s %+.6e %+.6e %+.6e\n",t,v.x,v.y,v.z); };
      std::printf("  [%s]\n",cf.name); pr("engine1",r.first); pr("engine2",r.second);
      pr("SO1",SO1); pr("SO2",SO2); pr("SS1",SS1); pr("SS2",SS2); pr("SO+SS1",SO1+SS1); pr("SO+SS2",SO2+SS2);
      std::printf("    common |(w1+w2)/2| engine %.6e\n",((r.first+r.second)*0.5).norm());
    }
  }
}
