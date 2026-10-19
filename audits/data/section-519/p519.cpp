#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 519: (A) Born <r^-3> vs the QM formula import 18 used; (B) electron-spin precession rate in hydrogen on the Langer orbit, as the
// transport takes it (secularRatesWithSStateRule + contact rule + applyBornInverseCube), vs the engine's fine-structure interval (import 21,
// lambda = g - 1) and the measured intervals.  Run with and without CREM_BORN_INVERSE_CUBE=1.
int main(){
  const double MHz=1e-6/(2*pi*hbar);
  std::printf("CREM_BORN_INVERSE_CUBE %s\n",std::getenv("CREM_BORN_INVERSE_CUBE")?"on":"off");
  for(const char* p:{"electron,positron","proton,electron"}){
    applyPairFromOption(p); const double aB=pairBohrRadius(activePair); double worst=0;
    for(int n=2;n<=12;++n) for(int l=1;l<n;++l){ const double a=n*n*aB, qm=1.0/(aB*aB*aB*n*n*n*l*(l+0.5)*(l+1.0));
      worst=std::max(worst,std::abs(bornInverseCubeRadius(a,(l+0.5)*hbar)/qm-1)); }
    std::printf("(A) %s: max |<r^-3>_Born/QM - 1| over n = 2..12, all l: %.2e\n",p,worst);
  }
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), mu=firstMass*secondMass/(firstMass+secondMass), lam=std::abs(secondGFactor)-1.0;
  struct C{int n,l; double meas,err; const char* src;};
  const C cs[]={{2,1,10969.045,0.015,"Hagley, Pipkin 1994 + Lundeen, Pipkin 1981"},{3,1,0,0,"-"},{3,2,0,0,"-"},{4,1,1371.104,0.018,"Berkeland et al. 1995"},
                {8,2,57.1342,0.0105,"de Beauvoir et al. 1997"},{12,2,16.9310,0.0117,"Schwob et al. 1999"}};
  for(const C& q:cs){
    const double a=q.n*q.n*aB, L=(q.l+0.5)*hbar; const Vec3 z{0,0,1}, x{1,0,0};
    const double interval=lam*(levelFineStructureEnergy(a,(q.l+1)*hbar)-levelFineStructureEnergy(a,q.l*hbar))*MHz;
    double wmin=1e300,wmax=0, along=0;
    for(double th:{20.0,60.0,100.0,150.0}){ const double t=th*pi/180; const Vec3 sh{std::sin(t),0,std::cos(t)};
      const Vec3 me=sh*(-secondMagneticMoment), mp=z*(firstMagneticMoment*1e-12);
      OrbitAveragedBmtAngularVelocities r=secularRatesWithSStateRule(a,z*L,mp,me,mu,0.0,x);
      applyContactAboveHbar(r,a,z*L,mp,me,mu); applyBornInverseCube(r,z*L);
      const double w=r.second.norm()*1e-6/(2*pi); wmin=std::min(wmin,w); wmax=std::max(wmax,w); along=dot(r.second,z)/r.second.norm(); }
    std::printf("(B) n %2d l %d: |Omega_e|/2pi %.4f .. %.4f MHz (cos(Omega, L) %+.6f); engine interval (import 21, lambda g-1) %.4f MHz; ratio %.6f;"
                " l(l+1)/(l+1/2)^2 = %.6f",q.n,q.l,wmin,wmax,along,interval,wmax/interval,q.l*(q.l+1.0)/((q.l+0.5)*(q.l+0.5)));
    if(q.meas>0) std::printf("; measured %.4f(%g) (%s): meas/Omega %.6f",q.meas,q.err,q.src,q.meas/wmax);
    std::printf("\n");
  }
}
