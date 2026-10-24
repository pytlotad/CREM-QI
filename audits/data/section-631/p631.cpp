// Audit 631, H631a: photonJWeights against exact 6j (printed for check631.py), label after the turn, and the drawn distribution.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
  for(int l=1;l<=12;++l) for(int dl: {-1,1}){ const int lp=l+dl; if(lp<1) continue;
    for(int J=l-1;J<=l+1;++J){ if(J<0) continue; const auto w=photonJWeights(l,J,lp); double t=0; for(double x: w) t+=x;
      std::printf("W l %d J %d lp %d",l,J,lp); for(double x: w) std::printf(" %.17e",x/t); std::printf("\n"); } }
  std::uint64_t st=splitMix64(631631ULL); long bad=0, n=0;
  const auto iso=[&](){ const double z=2*drawUniformUnit(st)-1, ph=2*pi*drawUniformUnit(st), s=std::sqrt(1-z*z); return Vec3{s*std::cos(ph),s*std::sin(ph),z}; };
  long counts[13][2][3][3]={};
  for(long i=0;i<300000;++i){ const int l=1+static_cast<int>(drawUniformUnit(st)*12); const int dl=drawUniformUnit(st)<0.5?-1:1; const int lp=l+dl; if(lp<1) continue;
    const Vec3 Lh=iso(), Sh=iso(); const int J=tripletJLabel(l,Lh,Sh);
    const int Jp=std::abs(lp-1)+formationDraw(photonJWeights(l,J,lp),drawUniformUnit(st));
    const double LQ=lp+0.5, target=std::acos(std::clamp(((Jp+0.5)*(Jp+0.5)-LQ*LQ-2.25)/(2*LQ*1.5),-1.0,1.0));
    const Vec3 S2=formationTurn(Sh,Lh,target,Sh); ++n; if(tripletJLabel(lp,Lh,S2)!=Jp) ++bad;
    if(l<=12) ++counts[l][dl>0][J-(l-1)][Jp-std::abs(lp-1)]; }
  std::printf("TURN draws %ld, label after turn != J' : %ld\n",n,bad);
  for(int l=1;l<=12;++l) for(int u=0;u<2;++u) for(int j=0;j<3;++j){ long t=0; for(int k=0;k<3;++k) t+=counts[l][u][j][k]; if(!t) continue;
    std::printf("C l %d up %d J %d n %ld",l,u,l-1+j,t); for(int k=0;k<3;++k) std::printf(" %ld",counts[l][u][j][k]); std::printf("\n"); }
}
