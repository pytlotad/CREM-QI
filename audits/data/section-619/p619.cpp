// Audit 619 T2: formation helpers on random isotropic configurations -- weights, rounding after the turn, 1S field shares.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <random>
int main(){
  std::mt19937_64 rng(619); std::normal_distribution<double> N(0,1); std::uniform_real_distribution<double> U(0,1);
  auto iso=[&](){ Vec3 v{N(rng),N(rng),N(rng)}; return v*(1.0/v.norm()); };
  const int M=200000;
  for(int l=1;l<=2;++l){   // Ps triplet J
    const double LQ=l+0.5; std::vector<long> cnt(3,0); long bad=0;
    for(int i=0;i<M;++i){ Vec3 Lh=iso(), Sh=iso(); const double th=std::acos(std::clamp(dot(Lh,Sh),-1.0,1.0));
      const int J=std::abs(l-1)+formationDraw(formationBornWeights(l,1.0,th),U(rng));
      const double tg=std::acos(std::clamp(((J+0.5)*(J+0.5)-LQ*LQ-2.25)/(2*LQ*1.5),-1.0,1.0));
      const Vec3 S2=formationTurn(Sh,Lh,tg,Sh); const Vec3 tot=Lh*LQ+S2*1.5;
      const int Jr=std::clamp((int)std::lround(tot.norm()-0.5),std::abs(l-1),l+1); if(Jr!=J) ++bad; ++cnt[J-std::abs(l-1)]; }
    std::printf("Ps l=%d: J weights",l); for(int k=0;k<3;++k){ const int J=std::abs(l-1)+k; if(J>l+1) break; std::printf("  J=%d %.4f (exact %.4f)",J,cnt[k]/double(M),(2*J+1)/(3.0*(2*l+1))); }
    std::printf("   rounding != drawn: %ld\n",bad); }
  for(int l=1;l<=2;++l){   // hydrogen j, F
    const double LQ=l+0.5; long cj[2]={0,0}, cf[2][2]={{0,0},{0,0}}, bad=0;
    for(int i=0;i<M;++i){ Vec3 Lh=iso(), sh=iso(), ih=iso();
      const double th=std::acos(std::clamp(dot(Lh,sh),-1.0,1.0)); const int ji=formationDraw(formationBornWeights(l,0.5,th),U(rng));
      const int k=l+ji; const double j=k-0.5; const double tg=std::acos(std::clamp((double(k)*k-LQ*LQ-1.0)/(2*LQ),-1.0,1.0));
      const Vec3 s2=formationTurn(sh,Lh,tg,sh); const Vec3 Jv=Lh*LQ+s2; const Vec3 Jh=Jv*(1.0/Jv.norm());
      const int kr=std::clamp((int)std::lround(Jv.norm()),l,l+1);
      const double th2=std::acos(std::clamp(dot(Jh,ih),-1.0,1.0)); const int fi=formationDraw(formationBornWeights(j,0.5,th2),U(rng));
      const double F=std::abs(j-0.5)+fi; const double tg2=std::acos(std::clamp(((F+0.5)*(F+0.5)-double(k)*k-1.0)/(2.0*k),-1.0,1.0));
      const Vec3 i2=formationTurn(ih,Jh,tg2,ih);
      const int Fr=std::clamp((int)std::lround((Jh*double(k)+i2).norm()-0.5),(int)std::lround(std::abs(j-0.5)),(int)std::lround(j+0.5));
      if(kr!=k||Fr!=(int)std::lround(F)) ++bad; ++cj[ji]; ++cf[ji][fi]; }
    std::printf("H l=%d: j weights %.4f / %.4f (exact %.4f / %.4f); F|j: ",l,cj[0]/double(M),cj[1]/double(M),l/(2.0*l+1),(l+1)/(2.0*l+1));
    for(int a=0;a<2;++a){ const double j=l-0.5+a; std::printf(" j=%.1f: F=%.1f %.4f (exact %.4f)",j,std::abs(j-0.5),cf[a][0]/double(cj[a]),(2*std::abs(j-0.5)+1)/(2*(2*j+1))); }
    std::printf("   rounding != drawn: %ld\n",bad); }
  // 1S in a field: product |s s> isotropic, Born on eigenstates via psFieldEigenAtPhoton
  const double h=2*pi*hbar, dE=203.3811e9*h, mu=9.284765e-24;
  for(double B: {0.01,0.1,1.0,3.0}){ long c[4]={0,0,0,0}; const int K=100000;
    for(int i=0;i<K;++i){ PsFieldEigenState st; psFieldFrame(Vec3{0,0,B},st.e1,st.e2,st.b); st.a1=mu*B; st.a2=-mu*B;
      const Vec3 sh=iso(); std::complex<double> ch[2]; psCoherentHalf(sh,ch); std::complex<double> p4[4];
      for(int a=0;a<2;++a) for(int b=0;b<2;++b) p4[2*a+b]=ch[a]*ch[b];
      for(int a=0;a<4;++a) for(int cc=0;cc<4;++cc) st.rho[a][cc]=std::real(p4[a]*std::conj(p4[cc]));
      gPsFieldEigen=st; const PsFieldEigenImage im=psFieldEigenAtPhoton(dE,U(rng),2*pi*U(rng));
      ++c[im.m==1?0:im.m==-1?1:(im.singletWeight<0.5?2:3)]; }
    std::printf("1S B=%.2f T: T+1 %.4f T-1 %.4f T0' %.4f S' %.4f  (M=0 block %.4f, exact 1/3)\n",B,c[0]/double(K),c[1]/double(K),c[2]/double(K),c[3]/double(K),(c[2]+c[3])/double(K)); }
}
