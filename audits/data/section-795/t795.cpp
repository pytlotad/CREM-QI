// Audit 795 unit check (S'): osculatingPeriapsisState for hydrogen and Ps on 8 tori, new start vs CREM_KEPLER_START=1.
// Prints the orbital energy (conserved - dipole-dipole - constraint) against the Sommerfeld energy E_S and the Kepler energy, the relativistic L,
// and the periapsis shift.  Dipoles: none, and Bohr-magneton moments parallel / antiparallel (spin energy must not move the orbit).
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  const double tori[8][2]={{3,1},{3,1.5},{3,2},{3,2.5},{3,3},{2,1},{2,1.5},{2,2}};
  const char* pairs[2]={"positron,electron","proton,electron"};
  for(const char* pp:pairs){ applyPairFromOption(pp);
    const double k=pairCoulombStrength/pairReducedMass, alpha=pairCoulombStrength/(hbar*c);
    for(int dip=0;dip<3;++dip) for(auto& t:tori){ const double n=t[0], L=t[1];
      OsculatingElements el; el.specificEnergy=-k/(2.0*n*n*pairBohrRadius(activePair)); el.specificAngularMomentum=L*hbar/pairReducedMass;
      const Vec3 z{0,0,1}, x{1,0,0}, m{0.6,0.0,0.8}; const double mu=9.2740100783e-24;
      const Vec3 d1=dip==0?Vec3{}:m*mu, d2=dip==0?Vec3{}:(dip==1?m*mu:m*(-mu));
      const double den=n-L+std::sqrt(L*L-alpha*alpha), ES=pairReducedMass*c*c*(1.0/std::sqrt(1.0+alpha*alpha/(den*den))-1.0);
      double eo[2], ln[2], rp[2];
      for(int kep=0;kep<2;++kep){ if(kep) setenv("CREM_KEPLER_START","1",1); else unsetenv("CREM_KEPLER_START");
        const State s=osculatingPeriapsisState(el,k,d1,d2,z,x,0.0);
        const Vec3 r=s.firstPosition-s.secondPosition;
        eo[kep]=conservedParticleEnergy(s)-pairDipoleInteractionEnergy(r,s.firstDipole,s.secondDipole)-s.dipoleConstraintEnergy;
        ln[kep]=(cross(s.firstPosition,momentum(s.firstVelocity,firstMass))+cross(s.secondPosition,momentum(s.secondVelocity,secondMass))).norm()/hbar;
        rp[kep]=r.norm(); }
      std::printf("%s dip %d n %g L %g | nowy: Eorb/E_S-1 %+.3e Lrel %.12f | Kepler: Eorb/E_S-1 %+.3e Lrel %.12f | r_new/r_Kep-1 %+.3e | E_S/E_nom-1 %+.3e\n",
        pp,dip,n,L,eo[0]/ES-1,ln[0],eo[1]/ES-1,ln[1],rp[0]/rp[1]-1,ES/(pairReducedMass*el.specificEnergy)-1); } }
}
