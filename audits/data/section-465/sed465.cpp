// Audit 465: Boyer-Puthoff power balance vs the SED Fokker-Planck of audit 387 -- direct simulation of a Kepler orbit
// (relative coordinate, units hbar = mu = e^2 = K = 1, c = 1/alpha) in a Gaussian isotropic zero-point field (dipole
// approximation, one-sided spectrum per component S(w) = (2/(3 pi)) alpha^3 w^3 in a band; first run had 4/(3 pi): osc-first-2x.txt), Landau-Lifshitz reaction
// tau d/dt(Coulomb), tau = (2/3) alpha^3, t0 = 1/tau.
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include <algorithm>
typedef std::complex<double> cd;
static void fft(std::vector<cd>& a){ const int n=int(a.size());
  for(int i=1,j=0;i<n;++i){ int bit=n>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit; if(i<j) std::swap(a[i],a[j]); }
  for(int len=2;len<=n;len<<=1){ const cd w=std::polar(1.0,-2*M_PI/len);
    for(int i=0;i<n;i+=len){ cd v=1; for(int j=0;j<len/2;++j){ cd x=a[i+j],y=a[i+j+len/2]*v; a[i+j]=x+y; a[i+j+len/2]=x-y; v*=w; } } } }
// ---- FP drift of audit 387 at (x, u) ------------------------------------------------------------------------
static void harmonics(double u,int KMAX,std::vector<double>& rp,std::vector<double>& rm){
  const int NF=1<<16; const double e=std::sqrt(std::max(0.0,1-u*u)); std::vector<cd> z(NF); double E=0;
  for(int j=0;j<NF;++j){ const double M=2*M_PI*j/NF; for(int it=0;it<60;++it){ const double f=E-e*std::sin(E)-M, d=1-e*std::cos(E), st=f/d; E-=st; if(std::abs(st)<1e-15) break; }
    z[j]={std::cos(E)-e,u*std::sin(E)}; }
  fft(z); rp.assign(KMAX+1,0); rm.assign(KMAX+1,0);
  for(int k=1;k<=KMAX;++k){ rp[k]=std::norm(z[k])/(double(NF)*NF); rm[k]=std::norm(z[NF-k])/(double(NF)*NF); } }
static void fpDrift(double x,double u,double& Ax,double& Ag,double& Rx,double& Rg){
  const int K=4000; const double du=1e-5; std::vector<double> p0,m0,p1,m1,p2,m2;
  harmonics(u,K,p0,m0); const double ulo=u-du, uhi=std::min(1.0,u+du); harmonics(ulo,K,p1,m1); harmonics(uhi,K,p2,m2);
  double V1x=0,V1g=0,V2x=0,V2g=0,V3x=0,V3g=0; const double g=u*x;
  for(int k=1;k<=K;++k) for(int s=-1;s<=1;s+=2){ const double rho=s>0?p0[k]:m0[k], drho=((s>0?p2[k]:m2[k])-(s>0?p1[k]:m1[k]))/(uhi-ulo);
    const double k3=double(k)*k*k, w=k3*rho, w2=k3*(-5.0*k*rho+drho*(s-k*u));
    V1x+=w*k; V1g+=w*s; V2x+=w2*k; V2g+=w2*s; V3x+=w*s*k; V3g+=w*s*s; }
  const double x5=std::pow(x,-5); Ax=x5*(-V1x+0.5*(V2x/x+V3x/g)); Ag=x5*(-V1g+0.5*(V2g/x+V3g/g)); Rx=-x5*V1x; Rg=-x5*V1g; }
// ---- noise ---------------------------------------------------------------------------------------------------
struct Noise{ std::vector<double> e[3]; double h; };
static void makeNoise(Noise& nz,int N,double h,double wlo,double whi,double alpha,double scale,std::mt19937_64& rng){
  std::normal_distribution<double> G(0,1); nz.h=h; const double dw=2*M_PI/(N*h);
  for(int c=0;c<3;++c){ std::vector<cd> A(N,0.0);
    for(int j=1;j<N/2;++j){ const double w=j*dw; if(w<wlo||w>whi) continue; const double S=2.0/(3*M_PI)*alpha*alpha*alpha*w*w*w*scale;
      const double sg=std::sqrt(S*dw/4); A[j]=cd(sg*G(rng),sg*G(rng)); }
    for(auto& v:A) v=std::conj(v); fft(A); nz.e[c].resize(N); for(int n=0;n<N;++n) nz.e[c][n]=2*std::conj(A[n]).real(); } }
static inline void field(const Noise& nz,double t,double* E){ if(nz.e[0].empty()){E[0]=E[1]=E[2]=0;return;}
  const double f=t/nz.h; int i=int(f); const int N=int(nz.e[0].size()); if(i>=N-1) i=N-2; const double w=f-i;
  for(int c=0;c<3;++c) E[c]=nz.e[c][i]*(1-w)+nz.e[c][i+1]*w; }
// ---- dynamics --------------------------------------------------------------------------------------------------
struct S6{ double r[3],v[3]; };
static inline void deriv(const S6& s,double t,const Noise& nz,double tau,bool osc,S6& d){
  double E[3]; field(nz,t,E); const double r2=s.r[0]*s.r[0]+s.r[1]*s.r[1]+s.r[2]*s.r[2], r=std::sqrt(r2);
  for(int c=0;c<3;++c) d.r[c]=s.v[c];
  if(osc){ for(int c=0;c<3;++c) d.v[c]=-s.r[c]+E[c]-tau*s.v[c]; return; }
  const double r3=r2*r, rv=s.r[0]*s.v[0]+s.r[1]*s.v[1]+s.r[2]*s.v[2];
  for(int c=0;c<3;++c) d.v[c]=-s.r[c]/r3+E[c]+tau*(-s.v[c]/r3+3*s.r[c]*rv/(r3*r2)); }
static void rk4(S6& s,double t,double h,const Noise& nz,double tau,bool osc){ S6 k1,k2,k3,k4,y;
  deriv(s,t,nz,tau,osc,k1); for(int c=0;c<3;++c){y.r[c]=s.r[c]+0.5*h*k1.r[c]; y.v[c]=s.v[c]+0.5*h*k1.v[c];}
  deriv(y,t+0.5*h,nz,tau,osc,k2); for(int c=0;c<3;++c){y.r[c]=s.r[c]+0.5*h*k2.r[c]; y.v[c]=s.v[c]+0.5*h*k2.v[c];}
  deriv(y,t+0.5*h,nz,tau,osc,k3); for(int c=0;c<3;++c){y.r[c]=s.r[c]+h*k3.r[c]; y.v[c]=s.v[c]+h*k3.v[c];}
  deriv(y,t+h,nz,tau,osc,k4); for(int c=0;c<3;++c){ s.r[c]+=h*(k1.r[c]+2*k2.r[c]+2*k3.r[c]+k4.r[c])/6; s.v[c]+=h*(k1.v[c]+2*k2.v[c]+2*k3.v[c]+k4.v[c])/6; } }
static void xg(const S6& s,double& x,double& g){ const double v2=s.v[0]*s.v[0]+s.v[1]*s.v[1]+s.v[2]*s.v[2], r=std::sqrt(s.r[0]*s.r[0]+s.r[1]*s.r[1]+s.r[2]*s.r[2]);
  const double E=0.5*v2-1/r; x=E<0?1/std::sqrt(-2*E):NAN; const double L[3]={s.r[1]*s.v[2]-s.r[2]*s.v[1],s.r[2]*s.v[0]-s.r[0]*s.v[2],s.r[0]*s.v[1]-s.r[1]*s.v[0]};
  g=std::sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]); }
int main(int argc,char**argv){
  const std::string mode=argv[1];
  if(mode=="fp"){ for(int a=2;a+1<argc;a+=2){ const double x=std::atof(argv[a]),u=std::atof(argv[a+1]); double Ax,Ag,Rx,Rg; fpDrift(x,u,Ax,Ag,Rx,Rg);
      std::printf("FP x=%.3f u=%.3f: A_x %.5f A_g %.5f (radiation only %.5f %.5f) [hbar/t0]; Puthoff circular -x^-5(1-1/x) = %.5f\n",x,u,Ax,Ag,Rx,Rg,-std::pow(x,-5)*(1-1/x)); } return 0; }
  const double alpha=std::atof(argv[2]); const double tau=2.0/3*alpha*alpha*alpha, t0=1/tau; const int Ns=std::atoi(argv[3]);
  if(mode=="osc"){ // 3D isotropic oscillator: stationary <E> (Boyer 3/2)
    const double h=2*M_PI/64, T=6*t0; int N=1; while(N*h<T+10) N<<=1; double sum=0,sum2=0;
    #pragma omp parallel for reduction(+:sum,sum2) schedule(dynamic)
    for(int i=0;i<Ns;++i){ std::mt19937_64 rng(465000+i); Noise nz; makeNoise(nz,N,h,0.2,5.0,alpha,1.0,rng);
      S6 s{{1,0,0},{0,1,0}}; double t=0, acc=0; long n=0; while(t<T){ rk4(s,t,h,nz,tau,true); t+=h; if(t>3*t0){ double E=0; for(int c=0;c<3;++c) E+=0.5*(s.v[c]*s.v[c]+s.r[c]*s.r[c]); acc+=E; ++n; } }
      const double m=acc/n; sum+=m; sum2+=m*m; }
    const double m=sum/Ns, se=std::sqrt((sum2/Ns-m*m)/Ns); std::printf("osc alpha=%.3f N=%d: <E> = %.4f +- %.4f (Boyer 1.5)\n",alpha,Ns,m,se); return 0; }
  // kepler x0 u0 tFinal[t0] zpf(0/1) wlo whi (in units of the orbital frequency) stepsPerOrbit
  const double x0=std::atof(argv[4]),u0=std::atof(argv[5]),tf=std::atof(argv[6])*t0; const bool zpf=std::atoi(argv[7]);
  const double wlo=std::atof(argv[8]),whi=std::atof(argv[9]); const int spo=std::atoi(argv[10]);
  const double a=x0*x0, e=std::sqrt(std::max(0.0,1-u0*u0)), nmean=std::pow(a,-1.5), Torb=2*M_PI/nmean, h=Torb/spo, w0=nmean;
  int N=1; while(N*h<tf+10*h) N<<=1;
  double sx1=0,sx2=0,sg1=0,sg2=0,qx=0,qg=0; int bad=0;
  #pragma omp parallel for reduction(+:sx1,sx2,sg1,sg2,qx,qg,bad) schedule(dynamic)
  for(int i=0;i<Ns;++i){ std::mt19937_64 rng(4650000+i); Noise nz; if(zpf) makeNoise(nz,N,h,wlo*w0,whi*w0,alpha,1.0,rng);
    std::uniform_real_distribution<double> U(0,2*M_PI); double M=U(rng),E=M; for(int it=0;it<60;++it) E-=(E-e*std::sin(E)-M)/(1-e*std::cos(E));
    const double fac=nmean*a/(1-e*std::cos(E)); S6 s{{a*(std::cos(E)-e),a*u0*std::sin(E),0},{-fac*std::sin(E),fac*u0*std::cos(E),0}};
    double t=0,xh=0,gh=0,xf=0,gf=0; const long nsteps=long(tf/h), half=nsteps/2;
    for(long k=0;k<nsteps;++k){ rk4(s,t,h,nz,tau,false); t+=h; if(k+1==half) xg(s,xh,gh); }
    xg(s,xf,gf); if(!std::isfinite(xf)||!std::isfinite(xh)){ ++bad; continue; }
    const double tt=(nsteps-half)*h/t0, dx=(xf-xh)/tt, dg=(gf-gh)/tt; sx1+=dx; sx2+=dx*dx; sg1+=dg; sg2+=dg*dg;
    double x0c,g0c; S6 s0{{a*(std::cos(E)-e),a*u0*std::sin(E),0},{-fac*std::sin(E),fac*u0*std::cos(E),0}}; xg(s0,x0c,g0c); qx+=xf-x0c; qg+=gf-g0c; }
  const int n=Ns-bad; const double mx=sx1/n, mg=sg1/n, ex=std::sqrt((sx2/n-mx*mx)/n), eg=std::sqrt((sg2/n-mg*mg)/n);
  std::printf("kepler alpha=%.3f x0=%.3f u0=%.3f zpf=%d band [%.2f,%.2f]w steps/orbit %d, t=%.3f t0 (%.0f orbits), N=%d (bad %d): A_x = %.4f +- %.4f, A_g = %.4f +- %.4f [hbar/t0]; total dx %.4f dg %.4f\n",
    alpha,x0,u0,int(zpf),wlo,whi,spo,tf/t0,tf/Torb,n,bad,mx,ex,mg,eg,qx/n,qg/n);
}
