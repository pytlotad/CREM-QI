#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 538: the engine's s-band rate vectors for Ps 1S in a field vs the two-spin model: Omega_1 = (2/hbar)(K f n2 + Z z), Omega_2 = (2/hbar)(K f n1 - Z z)
// (dn/dt = Omega x n up to the engine's sign convention), K = magnetic contact, f = 1.748, Z = (g/2) mu_B B.
int main(){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true; gSingletSpinTransport=false;
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  const double h=2*pi*hbar, muB=9.2740100783e-24, g=std::abs(firstGFactor), K=204.6577e9*h/1.748264/4, f=sStateContactTorqueFactor();
  const Vec3 z{0,0,1}, x{1,0,0};
  for(double B:{1.0,3.0}){ gExternalMagneticField=z*B; const double Z=(g/2)*muB*B;
    for(double th:{0.0,40.0,90.0}){ const double t=th*pi/180;
      const Vec3 n1{std::cos(t/2),0,std::sin(t/2)}, n2{std::cos(t/2),0,-std::sin(t/2)};   // spins split symmetrically about x in the xz plane
      const Vec3 m1=n1*(-firstMagneticMoment), m2=n2*secondMagneticMoment;
      OrbitAveragedBmtAngularVelocities r=secularRatesWithSStateRule(a,z*L,m1,m2,mu,0.0,x);
      const Vec3 i1=r.first-r.firstExternal, i2=r.second-r.secondExternal;
      const Vec3 mi1=n2*(2*K*f/hbar), me1=z*(2*Z/hbar);   // model magnitudes (sign convention free)
      std::printf("B %.0f T split %4.0f deg: engine internal |O1| %.4e |O2| %.4e (model %.4e), cos(O1_int, n2) %+.4f, cos(O2_int, n1) %+.4f;"
                  " external |O1| %.4e |O2| %.4e (model %.4e), cos(O1_ext,z) %+.4f cos(O2_ext,z) %+.4f\n",B,th,i1.norm(),i2.norm(),mi1.norm(),
                  dot(i1,n2)/i1.norm(),dot(i2,n1)/i2.norm(),r.firstExternal.norm(),r.secondExternal.norm(),me1.norm(),dot(r.firstExternal,z)/r.firstExternal.norm(),dot(r.secondExternal,z)/r.secondExternal.norm()); }
  }
}
