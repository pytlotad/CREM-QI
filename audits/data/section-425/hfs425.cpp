#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 425: HFS 1s (n = 1, L = hbar/2, QR contact, s-state isotropy) for spin-length variants V0-V3.
int main(){
  const double a1=pairBohrRadius(activePair), mu=firstMagneticMoment, h=2*pi*hbar, L=0.5*hbar;
  const double e=std::sqrt(1.0-0.25); const Vec3 ax{0,0,1}, P{1,0,0};
  auto U=[&](const Vec3& m1,const Vec3& m2){ return spinCouplingEnergy(a1,e,ax,P,m1,m2,L); };
  const Vec3 z{0,0,1};
  const double v0=(U(z*mu,z*(-mu))-U(z*mu,z*mu))/h*1e-9;
  const double r3=std::sqrt(3.0);
  const double v1=(U(z*(r3*mu),z*(-r3*mu))-U(z*(r3*mu),z*(r3*mu)))/h*1e-9;
  const double half=0.5*std::acos(1.0/3.0);   // half-angle of the triplet cone
  const Vec3 t1{std::sin(half),0,std::cos(half)}, t2{-std::sin(half),0,std::cos(half)};
  const double v2=(U(z*(r3*mu),z*(-r3*mu))-U(t1*(r3*mu),t2*(r3*mu)))/h*1e-9;
  std::printf("mu = %.9e J/T, a1 = %.9e m\n",mu,a1);
  std::printf("V0 (hbar/2, parallel triplet):          %.6f GHz\n",v0);
  std::printf("V1 (sqrt3/2 hbar, parallel triplet):    %.6f GHz  (x %.6f)\n",v1,v1/v0);
  std::printf("V2 (vector model, cone cos = 1/3):      %.6f GHz  (x %.6f)\n",v2,v2/v0);
  std::printf("V3 (contact x 2):                       %.6f GHz\n",2*v0);
  std::printf("measured 203.3942 GHz; leading-order magnetic part 4/7 x (7/12) alpha^2 Ry_inf... see text\n");
  // spin-orbit: L.S spread parallel vs antiparallel at Langer L = 1.5 hbar (units hbar^2); QM 3/2
  for(double s: {0.5,r3/2}) std::printf("spin-orbit spread 2 L S at L = 1.5, |S| = %.4f: %.4f hbar^2 (QM 1.5)\n",s,2*1.5*s);
  std::printf("vector model: L = sqrt2, S = sqrt3/2, j = 3/2: L.S = %.4f, j = 1/2: %.4f, spread %.4f\n",
    (3.75-2-0.75)/2,(0.75-2-0.75)/2,(3.75-0.75)/2);
  // annihilation weight under V2 triplet cone with w = (1 - s1.s2)/2
  std::printf("w(V2 triplet, (1 - s1.s2)/2) = %.6f; w = 1 - |S1+S2|^2/(2 hbar^2): triplet %.6f singlet %.6f\n",
    (1-1.0/3)/2, 1-2.0/2, 1-0.0);
}
