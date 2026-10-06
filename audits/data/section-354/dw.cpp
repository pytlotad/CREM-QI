#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
#include <algorithm>
// Audit 354: why w oscillates in p-Ps.  Delta omega_perp = |(omega1 - omega2) x mu_hat|, the rate at which the two moments
// lose their initial (anti)parallel alignment; split into orbital part (partner moment ~0) and partner-dipole part.
int main(){
  const double a1=pairBohrRadius(activePair), mu=firstMagneticMoment, mred=firstMass*secondMass/(firstMass+secondMass);
  std::mt19937_64 rng(354); std::normal_distribution<double> N(0,1);
  const Vec3 Lh{0,0,1}, P{1,0,0};
  std::printf("n = 1, a = %.4e m; Delta omega in rad/s; period 2 pi / Delta omega_perp in ps\n",a1);
  for(double L: {0.5,1.0}) for(int ch=1;ch<=2;++ch){
    std::vector<double> full,orb,part,par;
    for(int k=0;k<200;++k){
      Vec3 s{N(rng),N(rng),N(rng)}; s=s/s.norm();
      const Vec3 m1=s*mu, m2=s*(ch==1?mu:-mu);                 // p-Ps parallel moments, o-Ps antiparallel
      const Vec3 Lv=Lh*(L*hbar);
      auto perp=[&](const Vec3& d){ return (d-s*dot(d,s)).norm(); };
      const auto F=orbitAveragedBmtAngularVelocities(a1,Lv,m1,m2,mred,0.0,P);
      const auto O1=orbitAveragedBmtAngularVelocities(a1,Lv,m1,m2*1e-12,mred,0.0,P);   // orbital part of omega1
      const auto O2=orbitAveragedBmtAngularVelocities(a1,Lv,m1*1e-12,m2,mred,0.0,P);   // orbital part of omega2
      if(!F.valid||!O1.valid||!O2.valid) continue;
      const Vec3 dF=F.first-F.second, dO=O1.first-O2.second, dP=dF-dO;
      full.push_back(perp(dF)); orb.push_back(perp(dO)); part.push_back(perp(dP));
      par.push_back(std::abs(dot(F.first,s))+std::abs(dot(F.second,s)));
    }
    auto med=[](std::vector<double> v){ std::sort(v.begin(),v.end()); return v[v.size()/2]; };
    auto mx=[](std::vector<double> v){ return *std::max_element(v.begin(),v.end()); };
    std::printf("L=%.1f %s N=%zu  |dw_perp| med %.3e max %.3e (period med %.2f ps) | orbital part med %.3e | partner part med %.3e\n",
      L,ch==1?"p-Ps":"o-Ps",full.size(),med(full),mx(full),2*pi/med(full)*1e12,med(orb),med(part));
  }
}
