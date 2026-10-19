#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 535: Ps 1S in a field B along z, transport only: time-averaged contact weight w for o-Ps "m = 0" (parallel spins _|_ B), "m = +-1" (along B),
// and the p-Ps singlet with and without the singlet transport; vs QM mixing P_S = y^2/(1 + y^2) with the measured Delta E_HFS.
int main(){
  applyPairFromOption("electron,positron"); gSingletChannelEligible=true;
  const double aB=pairBohrRadius(activePair), a=aB, L=0.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass), e=std::sqrt(1-0.25);
  const double muB=9.2740100783e-24, g=std::abs(firstGFactor), dE=203.38910e9*2*pi*hbar;   // Ritter et al. 1984
  const Vec3 z{0,0,1}, x{1,0,0};
  { // model 1S HFS
    const Vec3 mt1=x*(-firstMagneticMoment), mt2=x*(-secondMagneticMoment)*(-1.0); // triplet parallel spins along x: s = sign(q) m^
    const Vec3 ms1=x*(-firstMagneticMoment), ms2=(x*(-1.0))*(secondMagneticMoment);
    const double Et=spinCouplingEnergy(a,e,z,x,mt1,mt2,L), Es=spinCouplingEnergy(a,e,z,x,ms1,ms2,L);
    std::printf("model 1S HFS (triplet - singlet) %.4f GHz; measured 203.38910(74) GHz (Ritter et al. 1984)\n",(Et-Es)/(2*pi*hbar)*1e-9);
  }
  const bool noST=std::getenv("CREM_NO_SINGLET_TRANSPORT")!=nullptr;
  for(double B:{0.1,0.3,1.0,3.0}){
    gExternalMagneticField=z*B;
    const double X=2*g*muB*B/dE, y=X/(1+std::sqrt(1+X*X)), PS=y*y/(1+y*y);
    struct C{const char* name; Vec3 s1,s2; bool singlet;};
    const C cs[]={{"o-Ps m=0 (s || x)",x,x,false},{"o-Ps m=+-1 (s || z)",z,z,false},{noST?"p-Ps singlet, no singlet transport":"p-Ps singlet, singlet transport",x,x*(-1.0),true}};
    std::printf("B %.1f T: x %.5f, QM P_S %.6e\n",B,X,PS);
    for(const C& q:cs){
      gSingletSpinTransport=q.singlet&&!noST;
      SecularSpinOrbitState s; s.orbitalAngularMomentum=z*L; s.periapsisDirection=x;
      s.firstDipole=q.s1*(-firstMagneticMoment); s.secondDipole=q.s2*(secondMagneticMoment);
      const double dt=2.5e-13, T=1e-9; double wsum=0, wmin=2, wmax=-1; long n=0;
      for(double t=0;t<T-1e-18;t+=dt){ s=advanceCoupledSecularSpinOrbit(s,a,mu,dt).state; const double w=contactTwoPhotonWeight(s.firstDipole,s.secondDipole); wsum+=w; ++n; wmin=std::min(wmin,w); wmax=std::max(wmax,w); }
      const double wm=wsum/n;
      std::printf("   %-38s <w> %.6e (min %.4e max %.4e); QM %s %.6e; ratio %.4f\n",q.name,wm,wmin,wmax,q.singlet?"1 - P_S":"P_S  ",q.singlet?1-PS:(q.s1.z>0.5?0.0:PS),
        q.singlet?(1-wm)/PS:(q.s1.z>0.5?0.0:wm/PS));
    }
  }
}
