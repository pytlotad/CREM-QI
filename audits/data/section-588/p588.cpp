// Audit 588: convention of contactTripletPolarization and the generator with an explicit triplet polarization.
#include "modules/contact_annihilation.hpp"
#include <cstdio>
int main(){
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  const auto dip=[&](const Vec3& s,double g){ return s*(g*hbar); };   // mu = g s (spin s in hbar)
  const Vec3 z{0,0,1};
  const auto report=[&](const char* name,const Vec3& s1,const Vec3& s2,double expectedP2){
    const ContactTripletPolarization t=contactTripletPolarization(dip(s1,g1),dip(s2,g2));
    double xi2=0; for(auto& x: t.xi) xi2+=std::norm(x);
    double dev=0;
    for(int i=0;i<=40;++i){ const double c=-1+i/20.0, sn=std::sqrt(std::max(0.0,1-c*c));
      for(int j=0;j<8;++j){ const double ph=2*pi*j/8; const Vec3 n{sn*std::cos(ph),sn*std::sin(ph),c};
        const double xn=std::norm(t.xi[0]*n.x+t.xi[1]*n.y+t.xi[2]*n.z)/xi2;
        const double ref=expectedP2<0?(1-c*c)/2:c*c; dev=std::max(dev,std::abs(xn-ref)); } }
    std::uint64_t stream=splitMix64(588);
    const int draws=200000; double p2=0,p22=0; int three=0;
    for(int k=0;k<draws;++k){ const ContactAnnihilationPhotons p=drawContactAnnihilationPhotons(0.0,z,stream,&t);
      if(p.count!=3) continue; ++three; const double c=p.normalSpinCosine; const double v=0.5*(3*c*c-1); p2+=v; p22+=v*v; }
    const double m=p2/three, sd=std::sqrt((p22/three-m*m)/three);
    std::printf("%-28s |xi|^2 %.15f  max ||xi.n|^2/|xi|^2 - ref| %.1e   <P2(n.z)> %.5f +- %.5f (expected %+.3f, %.1f sd)\n",
      name,xi2,dev,m,sd,expectedP2,(m-expectedP2)/sd);
  };
  report("s1 = s2 = z (m = +1)",z,z,-0.05);
  report("s1 = -s2 = z (T0)",z,z*-1.0,0.10);
}
