#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  const double a1=pairBohrRadius(activePair), mu=firstMagneticMoment, mred=firstMass*secondMass/(firstMass+secondMass);
  const Vec3 Lh{0,0,1}, P{1,0,0}, Q{0,1,0};
  struct D{const char* n; Vec3 s;} dirs[]={{"L^",Lh},{"P^",P},{"Q^",Q},{"45(P,Q)",(P+Q)/std::sqrt(2.0)},{"45(L,P)",(Lh+P)/std::sqrt(2.0)}};
  for(double L: {0.5,1.0}) for(auto& d: dirs){
    const Vec3 s=d.s/d.s.norm(), m=s*mu;
    const auto F=orbitAveragedBmtAngularVelocities(a1,Lh*(L*hbar),m,m,mred,0.0,P);
    const Vec3 dw=F.first-F.second; const double perp=(dw-s*dot(dw,s)).norm();
    std::printf("p-Ps L=%.1f mu||%-8s |dw_perp| = %.4e rad/s  |dw| = %.4e\n",L,d.n,perp,dw.norm());
  }
}
