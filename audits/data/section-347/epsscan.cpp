#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 347: L* (contact = measured p-Ps rate) versus the smearing radius, same quadrature as
// orbitAveragedContactDensity with the softening scaled; f=1 must reproduce the engine.
static double a1,psi0;
static double contactEps(double L,double soft){ const double e=std::min(std::sqrt(std::max(0.0,1-L*L)),1-1e-12);
  const int nodes=(int)std::clamp(512.0/std::sqrt(std::max(1.0-e,1e-12)),512.0,20001.0); double d=0,w=0;
  for(int k=0;k<nodes;++k){ const double u=(k+0.5)/nodes, an=pi*u*u*u*u, st=4*pi*u*u*u/nodes, tw=1-e*std::cos(an), r=a1*tw, rr=r*r+soft*soft;
    d+=tw*st*3*soft*soft/(4*pi*std::pow(rr,2.5)); w+=tw*st; } return d/w/psi0; }
int main(){ a1=pairBohrRadius(activePair); psi0=1/(pi*a1*a1*a1);
  const double target=0.99484, eps=magneticDipoleRadius();
  std::printf("check f=1, L=0.18: %.6e (engine %.6e)\n",contactEps(0.18,eps),orbitAveragedContactDensity(a1,std::sqrt(1-0.0324))/psi0);
  for(double f: {0.25,0.5,1.0,2.0,4.0}){ double lo=1e-4,hi=0.99; for(int k=0;k<100;++k){double m=.5*(lo+hi); (contactEps(m,f*eps)>target?lo:hi)=m;}
    const double L=.5*(lo+hi); std::printf("eps x %.2f: L* = %.5f hbar  rp/eps = %.3f  L*/f^(2/9) = %.5f\n",f,L,a1*(1-std::sqrt(1-L*L))/(f*eps),L/std::pow(f,2.0/9.0)); } }
