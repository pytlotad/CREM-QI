// Audit 810: energy probe (from p763) -- (1) j-splittings of the hydrogen level energy at n = 2..4 (as p763); (2) level energies beyond Bohr
// (levelEnergyBeyondBohr: fine + spin coupling) of 1S, 2S, 2P1/2, 2P3/2, 3S, 3P1/2, 3P3/2, 3D3/2, 3D5/2 with the proton spin parallel (F = 1 label)
// and antiparallel (F = 0) to the electron spin; (3) breitLevelFineStructureEnergy on a grid n = 1..4, L = 1..n step 0.25 (hbar) and the S
// branch (L = 0.5).  Compiles on the tree before 810 too: the remainder R is the new - new_off difference of (3), taken offline (ana810.py).
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  applyPairFromOption("proton,electron");
  const double h=2*pi*hbar, MHz=1.0/(h*1e6), aP=pairBohrRadius(activePair); const Vec3 z{0,0,1}, x{1,0,0};
  for(int n=2;n<=4;++n) for(int l=1;l<n;++l){
    const double a=n*n*aP, L=(l+0.5)*hbar; double E[2]; int K[2];
    for(int s=0;s<2;++s){ const Vec3 m=z*(s==0?1.0:-1.0)*9.2740100783e-24;
      K[s]=(int)std::lround(sommerfeldLevelAngularMomentum(z,m,m,L)/hbar); E[s]=sommerfeldLevelFineStructureEnergy(a,z,m,m,L); }
    const int hi=K[0]>K[1]?0:1;
    std::printf("split n %d l %d k %d %d  E(k=l+1) - E(k=l) = %.9f MHz\n",n,l,K[hi],K[1-hi],(E[hi]-E[1-hi])/h/1e6);
  }
  // levels: (name, n, l); electron spin +z / -z picks k = l + 1 / l (l >= 1); proton spin = +-electron spin
  struct Lv { const char* name; int n, l; int es; };
  const Lv lv[]={{"1S1/2",1,0,1},{"2S1/2",2,0,1},{"2P1/2",2,1,-1},{"2P3/2",2,1,1},{"3S1/2",3,0,1},{"3P1/2",3,1,-1},{"3P3/2",3,1,1},
                 {"3D3/2",3,2,-1},{"3D5/2",3,2,1}};
  for(const Lv& v: lv){
    const double a=v.n*v.n*aP, L=(v.l+0.5)*hbar, e=std::sqrt(std::max(0.0,1.0-std::pow((v.l+0.5)/v.n,2)));
    const Vec3 es=z*double(v.es);
    for(int F=1;F>=0;--F){
      const Vec3 ps=es*(F?1.0:-1.0);
      const Vec3 mp=ps*firstMagneticMoment, me=es*(-secondMagneticMoment);
      const double k=sommerfeldLevelAngularMomentum(z,mp,me,L)/hbar;
      const double E=levelEnergyBeyondBohr(a,e,z,x,mp,me,L);
      const double Ef=sommerfeldLevelFineStructureEnergy(a,z,mp,me,L);
      std::printf("level %s F %d k %.1f E %.9f MHz fine %.9f MHz\n",v.name,F,k,E*MHz,Ef*MHz);
    }
  }
  for(int n=1;n<=4;++n){
    std::printf("grid n %d L 0.50 breit %.17e J\n",n,breitLevelFineStructureEnergy(n*n*aP,0.5*hbar));
    for(int q=4;q<=4*n;++q) std::printf("grid n %d L %.2f breit %.17e J\n",n,q/4.0,breitLevelFineStructureEnergy(n*n*aP,q/4.0*hbar));
  }
}
