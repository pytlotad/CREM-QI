// Audit 811: (A) level energies beyond Bohr (levelEnergyBeyondBohr) and fine energies (sommerfeldLevelFineStructureEnergy) of all hydrogen
// levels n = 1..4 in J with %.17e, proton spin parallel (F = 1 label) / antiparallel (F = 0) to the electron spin -- same construction as p810;
// the 9 levels of 810 also in MHz with %.9f (control H811c).  (B) bornInverseCubeRadius of the engine for proton-electron and electron-positron,
// n = 2..4, l = 1..n-1, with the pair constants (masses, k, hbar, c) for the offline identity check.  No code change; reads only.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("proton,electron");
  const double h=2*pi*hbar, MHz=1.0/(h*1e6), aP=pairBohrRadius(activePair); const Vec3 z{0,0,1}, x{1,0,0};
  const char* ln="SPDF";
  for(int n=1;n<=4;++n) for(int l=0;l<n;++l) for(int es=-1;es<=1;es+=2){
    if(l==0&&es<0) continue;
    const int j2=l==0?1:(es>0?2*l+1:2*l-1);
    const double a=n*n*aP, L=(l+0.5)*hbar, e=std::sqrt(std::max(0.0,1.0-std::pow((l+0.5)/n,2)));
    const Vec3 esv=z*double(es);
    for(int F=1;F>=0;--F){
      const Vec3 ps=esv*(F?1.0:-1.0);
      const Vec3 mp=ps*firstMagneticMoment, me=esv*(-secondMagneticMoment);
      const double k=sommerfeldLevelAngularMomentum(z,mp,me,L)/hbar;
      const double E=levelEnergyBeyondBohr(a,e,z,x,mp,me,L);
      const double Ef=sommerfeldLevelFineStructureEnergy(a,z,mp,me,L);
      std::printf("levelJ %d%c%d/2 n %d l %d j2 %d F %d k %.1f E %.17e J fine %.17e J\n",n,ln[l],j2,n,l,j2,F,k,E,Ef);
      std::printf("levelMHz %d%c%d/2 F %d E %.9f MHz fine %.9f MHz\n",n,ln[l],j2,F,E*MHz,Ef*MHz);
    }
  }
  for(const char* pr: {"proton,electron","electron,positron"}){
    applyPairFromOption(pr); const double b=pairBohrRadius(activePair);
    std::printf("pair %s m1 %.17e m2 %.17e k %.17e hbar %.17e c %.17e aP %.17e\n",pr,firstMass,secondMass,pairCoulombStrength,hbar,c,b);
    for(int n=2;n<=4;++n) for(int l=1;l<n;++l)
      std::printf("born %s n %d l %d r3 %.17e\n",pr,n,l,bornInverseCubeRadius(n*n*b,(l+0.5)*hbar));
  }
}
