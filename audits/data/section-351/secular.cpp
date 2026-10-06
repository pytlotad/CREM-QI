#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
#include <algorithm>
// Audit 351: secular |L| mixing at n = 1.  Delaunay: dG/dt = -dU/dg, U = orbit-averaged dipole energy (engine function),
// G = |L|, g = argument of periapsis.  Rate and reachable range of G on the level U(G,g) = U(G0,g0), spins frozen.
static double a1;
static Vec3 rot(const Vec3& v,const Vec3& k,double t){ return v*std::cos(t)+cross(k,v)*std::sin(t)+k*(dot(k,v)*(1-std::cos(t))); }
static double U(double G,double g,const Vec3& Lh,const Vec3& P0,const Vec3& m1,const Vec3& m2){
  const double e=std::sqrt(std::max(0.0,1.0-G*G)); return orbitAveragedDipoleEnergy(a1,std::min(e,1.0-1e-9),Lh,rot(P0,Lh,g),m1,m2); }
int main(){
  a1=pairBohrRadius(activePair); const double mu=firstMagneticMoment;
  std::mt19937_64 rng(351); std::normal_distribution<double> N(0,1); std::uniform_real_distribution<double> Ug(0,2*pi);
  const Vec3 Lh{0,0,1}, P0{1,0,0};
  std::printf("n=1, a=%.4e m, |m|=%.4e J/T; tau_p=125.14 ps, tau_o=142.04 ns\n",a1,mu);
  std::printf("ch    L0     |dG/dt| med (hbar/ns)  t(0.5hbar)/ns med   reach dG med [min,max] (hbar)   |U| med (eV)\n");
  for(int ch=1;ch<=2;++ch) for(double L0: {0.05,0.1,0.2,0.3,0.5,0.7,0.9}){
    std::vector<double> rate,reach,uu;
    for(int k=0;k<24;++k){
      Vec3 s{N(rng),N(rng),N(rng)}; s=s/s.norm();
      const Vec3 m1=s*mu, m2=s*(ch==1?mu:-mu);       // p-Ps parallel (w=1), o-Ps antiparallel (w=0)
      const double g0=Ug(rng), h=1e-4;
      const double dUdg=(U(L0,g0+h,Lh,P0,m1,m2)-U(L0,g0-h,Lh,P0,m1,m2))/(2*h);
      rate.push_back(std::abs(dUdg)/hbar*1e-9);        // dG/dt in hbar per ns
      const double U0=U(L0,g0,Lh,P0,m1,m2); uu.push_back(std::abs(U0)/1.602176634e-19);
      // reachable G: walk outward from L0 while some g gives U = U0
      auto ok=[&](double G){ double lo=1e300,hi=-1e300; for(int j=0;j<48;++j){ double v=U(G,2*pi*j/48,Lh,P0,m1,m2); lo=std::min(lo,v); hi=std::max(hi,v);} return lo<=U0&&U0<=hi; };
      double up=L0,dn=L0; const double st=0.01;
      while(up+st<0.999&&ok(up+st)) up+=st; while(dn-st>0.005&&ok(dn-st)) dn-=st;
      reach.push_back(up-dn);
    }
    auto med=[](std::vector<double> v){ std::sort(v.begin(),v.end()); return v[v.size()/2]; };
    const double r=med(rate);
    std::printf("%s  %.2f   %.3e              %.3e           %.3f [%.3f,%.3f]              %.3e\n",ch==1?"p-Ps":"o-Ps",L0,r,0.5/r,
      med(reach),*std::min_element(reach.begin(),reach.end()),*std::max_element(reach.begin(),reach.end()),med(uu));
  }
}
