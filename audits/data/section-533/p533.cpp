#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 533: singlet contact energy (pair with a channel, antiparallel spins) -- new triangle form vs the old -3 x classical contact, Ps 1S and H 1S.
int main(){
  gSingletChannelEligible=true;
  for(const char* p:{"electron,positron","proton,electron"}){
    applyPairFromOption(p); const double aB=pairBohrRadius(activePair), MHz=1e-6/(2*pi*hbar); const Vec3 z{0,0,1}, x{1,0,0};
    const double a=aB, L=0.5*hbar, e=std::sqrt(1-0.25);
    for(double c:{-1.0,-0.999,-0.9}){   // cos between the two spins
      const Vec3 s1=z, s2=Vec3{std::sqrt(1-c*c),0,c};
      const Vec3 m1=s1*((firstCharge<0?-1:1)*firstMagneticMoment), m2=s2*((secondCharge<0?-1:1)*secondMagneticMoment);
      const double E=spinCouplingEnergy(a,e,z,x,m1,m2,L)*MHz;
      const double old=-3.0*(2.0/3.0)*mu0*dot(m1,m2)*quiggaRosnerContactDensity(a,L)*MHz;
      std::printf("%-18s 1S, cos(s1,s2) %+.3f: singlet energy new %.9f MHz, old %.9f MHz, rel. diff %.2e\n",p,c,E,old,(E-old)/old);
    }
  }
}
