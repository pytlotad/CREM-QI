#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
#include <map>
// Audit 520: import 22 (hyperfine of l >= 1 levels, hydrogen) -- (a) vector model without the quantum rules, (b) j, F basins for isotropic
// electron and proton spins, (c) Delta E(F) per j from spinCouplingEnergy vs the measured 2P1/2 interval.
static double proj(double LL,double SS,double JJ,double sr2,bool born){ // [L.J + 3 (S.r)^2 - S.J ... ] -> vector-model projection factor (g = 2, no recoil)
  const double LJ=0.5*(JJ+LL-SS), SJ=0.5*(JJ+SS-LL), SL=0.5*(JJ-LL-SS);
  if(born) return (LJ+3*sr2-SJ)/JJ;                       // r^.J = r^.S, (S.r^)^2 given
  return (LJ+0.5*SJ-1.5*SL*LJ/LL)/JJ;                      // classical orbit average 3<(S.r^)r^> = (3/2)(S - (S.L^)L^)
}
int main(){
  std::printf("(a) projection factor l(l+1)/(j(j+1)) (QM) vs vector model\n");
  for(auto [l,j]: std::vector<std::pair<int,double>>{{1,0.5},{1,1.5},{2,1.5},{2,2.5},{3,2.5},{3,3.5}}){
    const double qm=l*(l+1.0)/(j*(j+1.0));
    std::printf("  l %d j %.1f: QM %.6f; classical, rule-Q lengths %.4f (%.3f QM); classical, QM lengths %.4f (%.3f QM);"
                " Born squares + (S.r)^2 = 1/4: %.6f (%.2e)\n",l,j,qm,proj((l+.5)*(l+.5),1,(j+.5)*(j+.5),0,false),proj((l+.5)*(l+.5),1,(j+.5)*(j+.5),0,false)/qm,
                proj(l*(l+1.0),0.75,j*(j+1),0,false),proj(l*(l+1.0),0.75,j*(j+1),0,false)/qm,proj(l*(l+1.0),0.75,j*(j+1),0.25,true),proj(l*(l+1.0),0.75,j*(j+1),0.25,true)/qm-1);
  }
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), MHz=1e-6/(2*pi*hbar); const Vec3 z{0,0,1}, x{1,0,0};
  std::mt19937_64 r(520); std::uniform_real_distribution<double> U(0,1);
  auto iso=[&](){ const double c=2*U(r)-1, p=2*pi*U(r), s=std::sqrt(1-c*c); return Vec3{s*std::cos(p),s*std::sin(p),c}; };
  for(auto [n,l]: std::vector<std::pair<int,int>>{{2,1},{3,1},{3,2},{4,1}}){
    const double a=n*n*aB, L=(l+0.5)*hbar, e=std::sqrt(1-std::pow((l+0.5)/n,2)); const long N=400000;
    std::map<std::pair<int,int>,std::pair<long,double>> B; double spread=0;
    std::map<std::pair<int,int>,std::pair<double,double>> mm;
    for(long i=0;i<N;++i){ const Vec3 se=iso(), sp=iso(); const Vec3 me=se*(-secondMagneticMoment), mp=sp*firstMagneticMoment;
      const double E=spinCouplingEnergy(a,e,z,x,mp,me,L)*MHz;
      const Vec3 J=z*(l+0.5)+se; const int k=std::clamp((int)std::lround(J.norm()),l,l+1); const double j=k-0.5;
      const int F=std::clamp((int)std::lround((J*(k/J.norm())+sp).norm()-0.5),(int)std::lround(std::abs(j-0.5)),(int)std::lround(j+0.5));
      auto key=std::make_pair(2*k-1,F); B[key].first++; B[key].second+=E;
      auto& q=mm[key]; if(B[key].first==1) q={E,E}; q.first=std::min(q.first,E); q.second=std::max(q.second,E); }
    std::printf("(b,c) n %d l %d (%ld isotropic pairs of spin directions):\n",n,l,N);
    std::map<int,std::map<int,double>> E;
    for(auto& [key,v]:B){ const int j2=key.first, F=key.second; const double share=v.first/(double)N, qm=(2*F+1.0)/(4*(2*l+1.0));
      E[j2][F]=v.second/v.first; spread=std::max(spread,mm[key].second-mm[key].first);
      std::printf("   j %d/2 F %d: share %.4f (QM (2F+1)/(4(2l+1)) = %.4f), E %+.4f MHz\n",j2,F,share,qm,v.second/v.first); }
    for(auto& [j2,m]:E) if(m.size()==2){ auto it=m.begin(); const double lo=it->second; ++it;
      std::printf("   Delta E(F = j + 1/2 - F = j - 1/2), j = %d/2: %.4f MHz",j2,it->second-lo);
      if(n==2&&l==1&&j2==1) std::printf("  measured 59.22(14) MHz (Lundeen, Jessop, Pipkin 1975): residual %+.3f MHz (%+.2f sigma)",59.22-(it->second-lo),(59.22-(it->second-lo))/0.14);
      std::printf("\n"); }
    std::printf("   max spread of E within a basin %.2e MHz\n",spread);
  }
}
