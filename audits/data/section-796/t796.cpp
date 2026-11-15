// Audit 796 unit check: osculatingPeriapsisState for hydrogen and Ps on 8 tori -- new two-body start, CREM_ONEBODY_START=1 (795), CREM_KEPLER_START=1.
// Orbital energy (conserved - dipole-dipole - constraint) against E_T = E_S - nu mu c^2 alpha^4/(8 n^4) and E_S; canonical L computed from its
// definition, sum r_i x (gamma m v_i + dL_D/dv_i) with L_D = q1 q2/(2 c^2 r)[v1.v2 + (v1.n)(v2.n)]; prediction for 795: L_kan/L - 1 = nu alpha^2 a/r_p.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  const double tori[8][2]={{3,1},{3,1.5},{3,2},{3,2.5},{3,3},{2,1},{2,1.5},{2,2}};
  const char* pairs[2]={"positron,electron","proton,electron"};
  const char* modes[3]={"nowy","795","Kepler"};
  for(const char* pp:pairs){ applyPairFromOption(pp);
    const double k=pairCoulombStrength/pairReducedMass, alpha=pairCoulombStrength/(hbar*c), M=firstMass+secondMass, nu=firstMass*secondMass/(M*M);
    for(int dip=0;dip<3;++dip) for(auto& t:tori){ const double n=t[0], L=t[1];
      OsculatingElements el; el.specificEnergy=-k/(2.0*n*n*pairBohrRadius(activePair)); el.specificAngularMomentum=L*hbar/pairReducedMass;
      const Vec3 z{0,0,1}, x{1,0,0}, m{0.6,0.0,0.8}; const double mub=9.2740100783e-24;
      const Vec3 d1=dip==0?Vec3{}:m*mub, d2=dip==0?Vec3{}:(dip==1?m*mub:m*(-mub));
      const double den=n-L+std::sqrt(L*L-alpha*alpha), ES=pairReducedMass*c*c*(1.0/std::sqrt(1.0+alpha*alpha/(den*den))-1.0);
      const double ET=ES-nu*pairReducedMass*c*c*std::pow(alpha/n,4)/8.0;
      std::printf("%s dip %d n %g L %g",pp,dip,n,L);
      for(int md=0;md<3;++md){ unsetenv("CREM_ONEBODY_START"); unsetenv("CREM_KEPLER_START");
        if(md==1) setenv("CREM_ONEBODY_START","1",1); if(md==2) setenv("CREM_KEPLER_START","1",1);
        const State s=osculatingPeriapsisState(el,k,d1,d2,z,x,0.0);
        const Vec3 r=s.firstPosition-s.secondPosition; const double rn=r.norm(); const Vec3 nh=r*(1.0/rn);
        const double eo=conservedParticleEnergy(s)-pairDipoleInteractionEnergy(r,s.firstDipole,s.secondDipole)-s.dipoleConstraintEnergy;
        const double pref=coulomb*pairChargeProduct/(2.0*c*c*rn);
        const Vec3 a1=(s.secondVelocity+nh*dot(nh,s.secondVelocity))*pref, a2=(s.firstVelocity+nh*dot(nh,s.firstVelocity))*pref;
        const Vec3 cm=(s.firstPosition*firstMass+s.secondPosition*secondMass)*(1.0/M);
        const double lmech=(cross(s.firstPosition-cm,momentum(s.firstVelocity,firstMass))+cross(s.secondPosition-cm,momentum(s.secondVelocity,secondMass))).norm()/hbar;
        const double lkan=(cross(s.firstPosition-cm,momentum(s.firstVelocity,firstMass)+a1)+cross(s.secondPosition-cm,momentum(s.secondVelocity,secondMass)+a2)).norm()/hbar;
        const double pred=nu*alpha*alpha*pairBohrRadius(activePair)/rn;
        std::printf(" | %s: Eorb/E_T-1 %+.3e Eorb/E_S-1 %+.3e Lkan/L-1 %+.3e Lmech/L-1 %+.3e pred %.4e",modes[md],eo/ET-1,eo/ES-1,lkan/L-1,lmech/L-1,pred); }
      std::printf("\n"); } }
}
