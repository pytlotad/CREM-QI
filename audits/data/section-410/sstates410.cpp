// Audit 410: hydrogen ns lifetimes from the classical mode picture vs QM (theory: 3s 158.4, 4s 226, 5s 352, 6s 541 ns;
// branching 4s -> 2p 0.58 / 3p 0.42).  Mode amplitudes in closed form (Bessel), ABSOLUTE power per mode (no
// normalization by a truncated sum -- the error of audit 408's s-state numbers):
//   zeta = x + i y = sum_s c_s e^{i s M},  c_{+-s} = (A_s +- B_s)/2,  A_s = (2/s) J_s'(s e),  B_s = (2 sqrt(1-e^2)/(s e)) J_s(s e)
//   a circular component of radius |c| a at s omega radiates P = q^2 (s omega)^4 |c|^2 a^2/(6 pi eps0 c^3); components add.
// Check: sum over modes against the closed-form Larmor power.
// The s-state photon: the counter-rotating mode m = -1 (Delta l = +1); Delta n = s, n' >= 2 (no floor), or floored to n' = 2.
// Also: the same for p states (m = +1, floor n' = max(n-s, l'+1)) to re-check audit 408/409's p, d numbers.
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <algorithm>
static double Jp(int s,double x){ return 0.5*(std::cyl_bessel_j(s-1.0,x)-std::cyl_bessel_j(s+1.0,x)); }
static double cplus(int s,double e){ const double A=2.0/s*Jp(s,s*e), B=2.0*std::sqrt(1-e*e)/(s*e)*std::cyl_bessel_j(double(s),s*e); return 0.5*(A+B); }
static double cminus(int s,double e){ const double A=2.0/s*Jp(s,s*e), B=2.0*std::sqrt(1-e*e)/(s*e)*std::cyl_bessel_j(double(s),s*e); return 0.5*(A-B); }
int main(){
  const double hbar=1.054571817e-34, eps0=8.8541878128e-12, cl=299792458.0, qe=1.602176634e-19, me=9.1093837015e-31, mp=1.67262192369e-27;
  const double mu=me*mp/(me+mp), k=qe*qe/(4*M_PI*eps0), a0=hbar*hbar/(mu*k), C=qe*qe/(6*M_PI*eps0*cl*cl*cl);
  const auto larmor=[&](double e,double a){ return C*(k/mu)*(k/mu)*(1+e*e/2)/(std::pow(a,4)*std::pow(1-e*e,2.5)); };
  for(double e: {0.5,0.9,0.986}){ const double a=9*a0, w=std::sqrt(k/(mu*a*a*a)); double sum=0; const int smax=int(40/std::pow(1-e*e,1.5))+200;
    for(int s=1;s<=smax;++s){ const double cp=cplus(s,e), cm=cminus(s,e); sum+=C*std::pow(s*w,4)*(cp*cp+cm*cm)*a*a; }
    std::printf("check e=%.3f: mode sum (s<=%d) / closed-form Larmor = %.6f\n",e,smax,sum/larmor(e,a)); }
  const double qmS[7]={0,0,0,158.4,226.0,352.0,541.0};
  // QM (theory) s-state branchings: 4s -> 2p 0.58, 3p 0.42; 5s -> 2p 0.47, 3p 0.21, 4p 0.32 (approx.)
  // s states
  const auto tauS=[&](int n,double Lh,bool floor,double* br){ const double e=std::sqrt(1-(Lh/n)*(Lh/n)), a=n*n*a0, w=std::sqrt(k/(mu*a*a*a));
    double r=0; if(br) for(int i=0;i<8;++i) br[i]=0; const int smax=floor?int(40/std::pow(1-e*e,1.5))+200:n;
    for(int s=1;s<=smax;++s){ int n2=n-s; if(floor) n2=std::max(n2,2); if(n2<2||n2>=n) continue;
      const double cm=cminus(s,e); const double P=C*std::pow(s*w,4)*cm*cm*a*a; const double rate=P/(s*hbar*w); r+=rate; if(br) br[n2]+=rate; }
    if(br) for(int i=0;i<8;++i) br[i]/=r; return 1e9/r; };
  std::printf("s states, Langer L = hbar/2, photon = m = -1 mode:\n");
  for(int n=3;n<=6;++n){ double b[8]; const double t0=tauS(n,0.5,false,b), t1=tauS(n,0.5,true,nullptr);
    std::printf("  %ds: QM %.1f  no floor %.2f  floor %.3f ns  branches(no floor):",n,qmS[n],t0,t1); for(int i=2;i<n;++i) std::printf(" ->%dp %.2f",i,b[i]); std::printf("\n"); }
  std::printf("L_eff/hbar reproducing the QM s lifetime (no floor):\n");
  for(int n=3;n<=6;++n){ double best=0, err=1e9; for(int i=1;i<=2000;++i){ const double L=0.001*i; if(L>=n) break; const double t=tauS(n,L,false,nullptr);
      const double d=std::abs(std::log(t/qmS[n])); if(d<err){err=d;best=L;} }
    double b[8]; const double t=tauS(n,best,false,b); std::printf("  %ds: L_eff = %.3f hbar (tau %.1f ns, |log err| %.3f) branches:",n,best,t,err);
    for(int i=2;i<n;++i) std::printf(" ->%dp %.2f",i,b[i]); std::printf("\n"); }
  // both modes from l = 0 (m_l = 0: a photon of either m about L-hat leaves l' = 1), no floor, Langer L = hbar/2
  std::printf("s states, BOTH modes m = +-1 (l = 0 -> l' = 1 either way), no floor:\n");
  for(int n=3;n<=6;++n){ const double Lh=0.5, e=std::sqrt(1-(Lh/n)*(Lh/n)), a=n*n*a0, w=std::sqrt(k/(mu*a*a*a)); double r=0, br[8]={};
    for(int s=1;s<=n-2;++s){ const double cp=cplus(s,e), cm=cminus(s,e); const double rate=C*std::pow(s*w,4)*(cp*cp+cm*cm)*a*a/(s*hbar*w); r+=rate; br[n-s]+=rate; }
    std::printf("  %ds: QM %.1f  both modes %.2f ns (ratio %.3f)  branches:",n,qmS[n],1e9/r,(1e9/r)/qmS[n]); for(int i=2;i<n;++i) std::printf(" ->%dp %.2f",i,br[i]/r); std::printf("\n"); }
  // p, d, f states, m = +1, floor n' = max(n-s, l'+1)
  struct Lv{int n,l; double ref; const char* nm;};
  std::printf("l >= 1 (m = +1, floor), recheck of 408/409:\n");
  for(const Lv v: {Lv{2,1,1.600,"2p meas"},Lv{3,1,5.58,"3p meas"},Lv{3,2,15.6,"3d QM"},Lv{4,1,12.3,"4p QM"},Lv{4,2,36.5,"4d QM"},Lv{4,3,72.7,"4f QM"}}){
    const double Lh=v.l+0.5, e=std::sqrt(1-(Lh/v.n)*(Lh/v.n)), a=v.n*v.n*a0, w=std::sqrt(k/(mu*a*a*a)); double r=0; const int smax=int(40/std::pow(1-e*e,1.5))+200;
    for(int s=1;s<=smax;++s){ const int l2=v.l-1, n2=std::max(v.n-s,l2+1); if(n2>=v.n) continue; const double cp=cplus(s,e);
      r+=C*std::pow(s*w,4)*cp*cp*a*a/(s*hbar*w); }
    std::printf("  %s %.2f: %.4f ns\n",v.nm,v.ref,1e9/r); }
}
