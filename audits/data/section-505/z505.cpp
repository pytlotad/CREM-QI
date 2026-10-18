#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
// Audit 505: the engine's Delta m choice (zeemanPhotonChoice) for l -> l - 1 at Langer lengths vs QM Clebsch-Gordan branching.
static double cg2(int l,int m,int q){ // |<l-1, m-q; 1, q | l, m>|^2, j1 = l - 1, J = j1 + 1
  const double j1=l-1, M=m; if(std::abs(m-q)>l-1) return 0;
  if(q==1) return (j1+M)*(j1+M+1)/((2*j1+1)*(2*j1+2));
  if(q==0) return (j1-M+1)*(j1+M+1)/((2*j1+1)*(j1+1));
  return (j1-M)*(j1-M+1)/((2*j1+1)*(2*j1+2)); }
int main(){
  std::mt19937_64 r(505); std::uniform_real_distribution<double> U(0,1); const Vec3 z{0,0,1}; const long N=1000000;
  for(int l: {1,2,3}){
    const double Lb=(l+0.5)*hbar, La=(l-0.5)*hbar;
    std::vector<std::array<double,3>> cnt(2*l+1,{0,0,0}); std::vector<double> nb(2*l+1,0), nfin(2*l-1,0); double dm[3]={0,0,0}; long bad=0, inval=0;
    for(long i=0;i<N;++i){ const double cz=2*U(r)-1, ph=2*pi*U(r), sz=std::sqrt(1-cz*cz); const Vec3 Lh{sz*std::cos(ph),sz*std::sin(ph),cz};
      const ZeemanPhotonChoice ch=zeemanPhotonChoice(Lh,Lb,La,z,U(r)); if(!ch.valid){ ++inval; continue; }
      const int m=(int)std::lround(Lh.z*Lb/hbar), q=(int)ch.deltaM, mf=(int)std::lround(ch.afterHat.z*La/hbar);
      nb[m+l]++; cnt[m+l][q+1]++; dm[q+1]++; if(mf!=m-q) ++bad; if(std::abs(mf)<=l-1) nfin[mf+l-1]++; else ++bad; }
    std::printf("l %d -> %d (|L| %.1f -> %.1f hbar), N %ld, invalid %ld, final-basin mismatches %ld\n",l,l-1,l+0.5,l-0.5,N,inval,bad);
    double maxd=0;
    for(int m=-l;m<=l;++m){ std::printf("  m_l %+d (share %.4f, QM %.4f): P(dm=-1,0,+1) engine %.4f %.4f %.4f | QM %.4f %.4f %.4f\n",m,nb[m+l]/N,1.0/(2*l+1),
        cnt[m+l][0]/nb[m+l],cnt[m+l][1]/nb[m+l],cnt[m+l][2]/nb[m+l],cg2(l,m,-1),cg2(l,m,0),cg2(l,m,1));
      for(int q=-1;q<=1;++q) maxd=std::max(maxd,std::abs(cnt[m+l][q+1]/nb[m+l]-cg2(l,m,q))); }
    std::printf("  total P(dm=-1,0,+1) engine %.4f %.4f %.4f (QM 1/3 each); final m_l':",dm[0]/N,dm[1]/N,dm[2]/N);
    for(int k=0;k<2*l-1;++k) std::printf(" %.4f",nfin[k]/N); std::printf(" (QM %.4f each); max |P_engine - P_QM| per basin %.4f\n",1.0/(2*l-1),maxd);
  }
}
