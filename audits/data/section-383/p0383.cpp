#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <random>
// Audit 383: photon-count share of Kepler harmonics vs the engine table; lifetimes of hydrogen-like states under four rules.
static double Jd(int k,double x){ return 0.5*(std::cyl_bessel_j(k-1,x)-std::cyl_bessel_j(k+1,x)); }
static double wk(int k,double e){ const double x=k*e, J=std::cyl_bessel_j(k,x), D=Jd(k,x); const double v=k*(D*D+(1-e*e)/(e*e)*J*J); return std::isfinite(v)?v:0.0; }
static std::vector<double> shares(double e,int K){ std::vector<double> w(K+1,0.0); double s=0; for(int k=1;k<=K;++k){ w[k]=wk(k,e); s+=w[k]; } for(auto& v:w) v/=s; return w; }
static double tableP1(double e){ std::mt19937_64 r(1); std::uniform_real_distribution<double> U(0,1); int c=0,N=400000;
  for(int i=0;i<N;++i) if(std::lround(eccentricOrbitHarmonicNumber(e,U(r)))<=1) ++c; return c/(double)N; }
static std::vector<double> tableShares(double e,int K){ std::mt19937_64 r(2); std::uniform_real_distribution<double> U(0,1); std::vector<double> w(K+1,0.0); int N=400000;
  for(int i=0;i<N;++i){ long k=std::lround(eccentricOrbitHarmonicNumber(e,U(r))); if(k<1)k=1; if(k>K)k=K; w[k]+=1.0/N; } return w; }
int main(){
  // pair: hydrogen (proton, electron) via reduced mass; Ps via the engine's active pair
  const double me=9.1093837015e-31, mp=1.67262192369e-27, K=pairCoulombStrength;   // e^2/(4 pi eps0)
  const double muPs=firstMass*secondMass/(firstMass+secondMass), muH=me*mp/(me+mp);
  const auto rate=[&](double mu,double n,double Lh){ const double a=n*n*hbar*hbar/(mu*K), e2=1-(Lh/n)*(Lh/n);
    const double r3=1/(a*a*a*std::pow(e2<1?1-e2:1e-300,1.5)); return 2*K*K*(Lh*hbar)*r3/(3*c*c*c*mu*mu*hbar); };
  // H383a/b
  std::printf("photon-count share of harmonic 1 (exact vs engine table):\n");
  for(double e: {0.5528,0.6614,0.8660,0.9682,0.9860}){ const int Kmax=e>0.98?60000:(e>0.9?8000:800); auto w=shares(e,Kmax); auto w2=shares(e,Kmax/2); std::printf("  (convergence: w1 at Kmax/2 %.6f)",w2[1]);
    std::printf("  e=%.4f  exact w1=%.5f w2=%.5f   table P(k=1)=%.5f\n",e,w[1],w[2],tableP1(e)); }
  { const double n=3,Lh=0.5,a=n*n*hbar*hbar/(muPs*K), T=2*pi*std::sqrt(muPs*a*a*a/K);
    std::printf("Ps 3s: R = %.4e /s, per orbit %.4e (engine CREM_L_BALANCE: 3.289e-06), T = %.4e s\n",rate(muPs,n,Lh),rate(muPs,n,Lh)*T,T); }
  // lifetimes: allowed k for (n,l): l'=l-1 (l>=1) or 1 (l=0); n'=n-k must have n' >= l'+1
  struct St{int n,l; const char* qm;}; St sts[]={{2,1,"1.596"},{3,0,"158"},{3,1,"5.40"},{3,2,"15.6"},{4,0,"226"},{4,1,"12.4"},{4,2,"36.5"},{4,3,"72.6"},{5,1,"23.6"}};
  std::printf("\nhydrogen lifetimes [ns]  (QM theory reference, not measurement; measured: 2p 1.600(4), 3p 5.58(13), 4p 11.25(78), 5p 21.9(3.0))\n");
  std::printf(" state  e       V0 engine V1 exact  V2 trim   V3 Bohr     QM      | V1 branching by n'\n");
  for(auto s: sts){ const double Lh=s.l+0.5, e=std::sqrt(1-(Lh/s.n)*(Lh/s.n)); const int lp=s.l>=1?s.l-1:1;
    const int Kmax=e>0.98?60000:(e>0.9?8000:2000); auto w=shares(e,Kmax); auto wt=tableShares(e,Kmax);
    const double R=rate(muH,s.n,Lh); double a1=0,a0=0,ap=0; std::string br; std::vector<double> toN(s.n,0.0);
    for(int k=1;k<=Kmax;++k){ const int np=std::max(s.n-k,1); if(np>=lp+1){ a1+=w[k]; a0+=wt[k]; toN[np]+=w[k]; if(s.n-k>=1) ap+=w[k]; } }
    char buf[256]; for(int np=s.n-1;np>=1;--np) if(toN[np]>0){ std::snprintf(buf,sizeof buf," n'=%d %.3f",np,toN[np]/a1); br+=buf; }
    const double v3=ap; // pure correspondence: only k <= n-1, allowed
    // V4 (s states only): Kramers -- rate and harmonic share on the intermediate orbit n - k/2, L = (L + L')/2 = hbar
    double v4=0; if(s.l==0) for(int k=1;k<s.n;++k){ const int np=s.n-k; if(np<2) continue; const double nb=s.n-0.5*k, Lb=1.0;
      const double eb=std::sqrt(std::max(1e-12,1-(Lb/nb)*(Lb/nb))); auto wb=shares(eb,eb>0.98?60000:(eb>0.9?8000:2000)); v4+=rate(muH,nb,Lb)*wb[k]; }
    if(s.l==0) std::printf("   V4 Kramers (s): %.3f ns\n",1e9/v4);
    std::printf(" %d%c   %.4f  %8.3f  %8.3f  %8.3f  %8.3f    %-6s  |%s\n",s.n,"spdf"[s.l],e,1e9/(R*a0),1e9/(R*a1),1e9/R,v3>0?1e9/(R*v3):1e300,s.qm,br.c_str());
  }
}
