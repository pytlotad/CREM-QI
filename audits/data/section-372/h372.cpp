#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 372: analytic semiclassical np lifetimes (Langer and sqrt(l(l+1))) for hydrogen n = 2..6.
int main(){
  applyPairFromOption("proton,electron");
  const double mred=firstMass*secondMass/(firstMass+secondMass), k=pairCoulombStrength, a0=pairBohrRadius(activePair), q=std::sqrt(k*4*pi*epsilon0);
  auto f=[](double e){ return std::pow(1-e*e,1.5)/(1+e*e/2); };
  auto tau=[&](int n,double Lh){ const double a=n*n*a0, w=std::sqrt(k/(mred*a*a*a)), e=std::sqrt(1-std::pow(Lh/n,2));
    const double P=q*q*k*k*(1+e*e/2)/(std::pow(a,4)*std::pow(1-e*e,2.5))/(6*pi*epsilon0*c*c*c*mred*mred); return hbar*w/(f(e)*P); };
  struct M{int n; double v,s; const char* src;} meas[]={{2,1.600,0.004,"BG66"},{3,5.58,0.13,"BG66"},{3,5.41,0.18,"E70"},{4,11.25,0.78,"E70"},{5,21.9,3.0,"E70"}};
  for(auto& m: meas){ const double tl=tau(m.n,1.5)*1e9, ts=tau(m.n,std::sqrt(2.0))*1e9;
    std::printf("%dp (%s) measured %.3f +- %.3f ns | Langer %.4f ns (%+.1f%%, %+.2f sigma) | sqrt(l(l+1)) %.4f ns (%+.1f%%, %+.2f sigma)\n",
      m.n,m.src,m.v,m.s,tl,100*(tl/m.v-1),(tl-m.v)/m.s,ts,100*(ts/m.v-1),(ts-m.v)/m.s); }
}
