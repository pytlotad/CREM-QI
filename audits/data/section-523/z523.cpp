#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
// Audit 523: zeemanPhotonChoice with the mode sign.  (A) integer M = m hbar: weights vs CG^2, Delta l = -1 (sign +1) and +1 (sign -1), l = 0..6;
// (B) 1e6 isotropic L: 2p, 3d, 4f (sign +1) and 3s -> 2p (|L| 1/2 -> 3/2, sign -1), basin shares vs QM; (C) axis = L^ (no field) with the same rule.
static double cgd(int l,int m,int q){ return q==1?(l+m)*(l+m-1.0):(q==0?2.0*(l*l-m*m):(l-m)*(l-m-1.0)); }
static double cgu(int l,int m,int q){ return q==1?(l-m+1.0)*(l-m+2.0):(q==0?2.0*((l+1)*(l+1)-m*m):(l+m+1.0)*(l+m+2.0)); }
int main(){
  const Vec3 z{0,0,1};
  double worst[2]={0,0};
  for(int dl:{-1,1}) for(int l=(dl<0?1:0);l<=6;++l) for(int m=-l;m<=l;++m){
    const double Ln=(l+0.5)*hbar, La=(l+dl+0.5)*hbar, cz=m*hbar/Ln, sz=std::sqrt(std::max(0.0,1-cz*cz));
    const Vec3 Lh{sz,0,cz}; ZeemanPhotonChoice ch=zeemanPhotonChoice(Lh,Ln,La,z,0.5,dl<0?1.0:-1.0);
    double ref[3], s=0; for(int k=0;k<3;++k){ const int q=k-1; ref[k]=std::abs(m-q)<=l+dl?(dl<0?cgd(l,m,q):cgu(l,m,q)):0; s+=ref[k]; }
    for(int k=0;k<3;++k) worst[dl>0]=std::max(worst[dl>0],std::abs(ch.weight[k]-ref[k]/s)); }
  std::printf("(A) integer M: max |w - CG^2| Delta l = -1 (l = 1..6): %.2e; Delta l = +1 (l = 0..6): %.2e\n",worst[0],worst[1]);
  std::mt19937_64 r(523); std::uniform_real_distribution<double> U(0,1); const long N=1000000;
  struct T{double Lb,La,sg; const char* name;}; const T ts[]={{1.5,0.5,1,"2p -> 1s"},{2.5,1.5,1,"3d -> 2p"},{3.5,2.5,1,"4f -> 3d"},{0.5,1.5,-1,"3s -> 2p (Delta l = +1)"}};
  for(const T& tt:ts){ double cnt[3]={0,0,0}; long inval=0; double maxdev=0; const int l=(int)std::lround(tt.Lb-0.5);
    std::vector<std::array<double,4>> bas(2*l+1,{0,0,0,0});
    for(long i=0;i<N;++i){ const double cz=2*U(r)-1, ph=2*pi*U(r), sz=std::sqrt(1-cz*cz); const Vec3 Lh{sz*std::cos(ph),sz*std::sin(ph),cz};
      const ZeemanPhotonChoice ch=zeemanPhotonChoice(Lh,tt.Lb*hbar,tt.La*hbar,z,U(r),tt.sg); if(!ch.valid){++inval;continue;}
      const int q=(int)ch.deltaM; cnt[q+1]++; const int m=(int)std::lround(Lh.z*tt.Lb); bas[m+l][q+1]++; bas[m+l][3]++; }
    for(int m=-l;m<=l;++m){ double s=0, ref[3]; for(int k=0;k<3;++k){ const int q=k-1; ref[k]=std::abs(m-q)<=(int)std::lround(tt.La-0.5)?(tt.sg>0?cgd(l,m,q):cgu(l,m,q)):0; s+=ref[k]; }
      for(int k=0;k<3;++k) maxdev=std::max(maxdev,std::abs(bas[m+l][k]/bas[m+l][3]-ref[k]/s)); }
    std::printf("(B) %-24s P(dm = -1, 0, +1) = %.4f %.4f %.4f; max |P - CG^2| per basin %.4f; invalid %ld\n",tt.name,cnt[0]/N,cnt[1]/N,cnt[2]/N,maxdev,inval); }
  for(const T& tt:ts){ const ZeemanPhotonChoice ch=zeemanPhotonChoice(z,tt.Lb*hbar,tt.La*hbar,z,0.5,tt.sg);
    std::printf("(C) axis = L^, %-24s w(dm = -1, 0, +1) = %.4f %.4f %.4f  (import 2: m = %+d)\n",tt.name,ch.weight[0],ch.weight[1],ch.weight[2],(int)tt.sg); }
}
