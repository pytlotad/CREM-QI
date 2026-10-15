// Audit 461: |J| attractor on the spin-orbit precession loop (424/453 law moved to the precession angle), Ps 2p-4f.
// Omega_c, N_c from section-460/dephase460.txt; Ps orbit period T_n = n^3 x 0.303966 fs (453 regression).
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
static const double PI=3.14159265358979323846;
struct St{int n,l; double Oc,Nc;};
static double slope(double X,double lam,double q){              // dX/dN_p
  if(lam==0) return 0; const double th=2*PI*(X-0.5), g=lam*(1-q)*(1-q)/(2*PI*q);
  return -g*q*std::sin(th)/(1-2*q*std::cos(th)+q*q);
}
static double evolve(double X,double N,double lam,double q,double lo,double hi){
  if(lam==0||N<=0) return X; const int steps=std::max(4,(int)std::ceil(lam*N/0.05)); const double h=N/steps;
  for(int i=0;i<steps;++i){ const double k1=slope(X,lam,q),k2=slope(X+0.5*h*k1,lam,q),k3=slope(X+0.5*h*k2,lam,q),k4=slope(X+h*k3,lam,q);
    X=std::clamp(X+h*(k1+2*k2+2*k3+k4)/6,lo,hi);} return X;
}
static bool nearRung(double X){ return std::abs(X-(std::round(X-0.5)+0.5))<0.05; }
int main(){
  std::vector<St> st; FILE* f=std::fopen("../section-460/dephase460.txt","r"); char line[512];
  while(std::fgets(line,sizeof line,f)){ int n,l; double Oc,Od1,Od2,tau,Nc;
    if(std::sscanf(line,"%d %d %lf %lf %lf %lf %lf",&n,&l,&Oc,&Od1,&Od2,&tau,&Nc)==7&&l>=1) st.push_back({n,l,Oc,Nc}); }
  std::fclose(f);
  const double tau=3, q=std::exp(-1/tau); const double lams[]={0,1e-3,1e-2,0.1,0.3,1,3}; const int M=1000;
  std::printf("# H461a/b: model spin length |S| = 1 hbar, isotropic cos(L,S), tau = 3 precession loops\n");
  std::printf("# state  N_c  orbits/precession  tau_orbits | fraction within +-0.05 of a rung at N = N_c, and decay-averaged, per lambda\n");
  for(auto s: st){
    const double T=s.n*s.n*s.n*0.303966e-15, opp=(2*PI/s.Oc)/T, L=s.l+0.5, S=1.0, lo=std::abs(L-S), hi=L+S;
    std::printf("%d%c  %.2f  %.3e  %.3e |",s.n,"spdf"[s.l],s.Nc,opp,tau*opp);
    for(double lam: lams){ int a=0; double avg=0;
      for(int i=0;i<M;++i){ const double c=-1+2*(i+0.5)/M, X0=std::sqrt(L*L+S*S+2*L*S*c);
        if(nearRung(evolve(X0,s.Nc,lam,q,lo,hi))) ++a;
        // decay-averaged: N ~ exponential with mean N_c (midpoint quadrature in the survival fraction)
        // N sorted ascending (u descending); evolve incrementally
        const int K=40; double w=0, X=X0, Nprev=0; for(int k=K-1;k>=0;--k){ const double u=(k+0.5)/K, N=-s.Nc*std::log(u); X=evolve(X,N-Nprev,lam,q,lo,hi); Nprev=N; w+=nearRung(X); } avg+=w/K; }
      std::printf("  lam=%g: %.3f/%.3f",lam,double(a)/M,avg/M); }
    std::printf("\n");
  }
  std::printf("# H461c: multiplet populations from the basins (J = clamp(round(|J|/hbar - 1/2), |l-1|, l+1)), isotropic preparation, vs QM (2J+1)/(3(2l+1))\n");
  const double lens[]={1.0,std::sqrt(2.0),1.5,2.0}; const char* nm[]={"model 1","sqrt(S(S+1))","Langer 1.5","two EBK 2"};
  for(int l=1;l<=3;++l) for(int k=0;k<4;++k){ const double L=l+0.5,S=lens[k]; double p[3]={0,0,0};
    // analytic: density of X is X/(2LS) on [|L-S|, L+S]; integrate over each basin
    const double lo=std::abs(L-S),hi=L+S;
    for(int J=l-1;J<=l+1;++J){ double a=(J==l-1)?-1e9:J, b=(J==l+1)?1e9:J+1; a=std::max(a,lo); b=std::min(b,hi); if(b>a) p[J-l+1]=(b*b-a*a)/(4*L*S); }
    double dmax=0; std::printf("l=%d  |S|=%-13s",l,nm[k]);
    for(int J=l-1;J<=l+1;++J){ const double qm=(2*J+1)/(3.0*(2*l+1)); dmax=std::max(dmax,std::abs(p[J-l+1]-qm)); std::printf("  J=%d %.6f (QM %.6f)",J,p[J-l+1],qm); }
    std::printf("  max|diff| %.3e\n",dmax); }
  // Production-like ensemble (crem_trajectory.hpp:2035/2120): both spin directions isotropic, import 18 acts when spin cos >= 0
  // (not antiparallel); S = (s1 + s2) hbar/2, L = (l + 1/2) hbar along z; J = clamp(round(|L + S|/hbar - 1/2)).
  std::printf("# H461c (production): isotropic s1, s2 with s1.s2 >= 0, |S| = |s1 + s2|/2 hbar in [0.707, 1]; import 18 rounding\n");
  std::mt19937_64 rng(461); std::normal_distribution<double> g(0,1);
  for(int l=1;l<=3;++l){ const double L=l+0.5; long cnt[3]={0,0,0}, tot=0; double Ssum=0;
    while(tot<2000000){ double a[3],b[3]; double na=0,nb=0; for(int i=0;i<3;++i){a[i]=g(rng);b[i]=g(rng);na+=a[i]*a[i];nb+=b[i]*b[i];}
      na=std::sqrt(na); nb=std::sqrt(nb); double c=0; for(int i=0;i<3;++i){a[i]/=na;b[i]/=nb;c+=a[i]*b[i];} if(c<0) continue;
      const double Sx=0.5*(a[0]+b[0]),Sy=0.5*(a[1]+b[1]),Sz=0.5*(a[2]+b[2]); Ssum+=std::sqrt(Sx*Sx+Sy*Sy+Sz*Sz);
      const double X=std::sqrt(Sx*Sx+Sy*Sy+(L+Sz)*(L+Sz)); const int J=std::clamp((int)std::lround(X-0.5),l-1,l+1); ++cnt[J-l+1]; ++tot; }
    double dmax=0; std::printf("l=%d  <|S|> %.4f",l,Ssum/tot);
    for(int J=l-1;J<=l+1;++J){ const double p=double(cnt[J-l+1])/tot, qm=(2*J+1)/(3.0*(2*l+1)); dmax=std::max(dmax,std::abs(p-qm)); std::printf("  J=%d %.4f (QM %.4f)",J,p,qm);} 
    std::printf("  max|diff| %.3f  (MC 2e6, sigma ~3e-4)\n",dmax); }
}
