// Audit 387: SED Fokker-Planck of a Kepler orbit in action variables (x = J/hbar, g = |L|/hbar), Ito Monte Carlo.
// Channels (k, s): R = x^-5 k^3 rho_ks(u), u = g/x; diffusion (hbar/2) R n n^T, n = (k, s).  Time unit t0 = 1/((2/3) alpha^3 mu K^2/hbar^3).
#include <cmath>
#include <cstdio>
#include <vector>
#include <random>
#include <algorithm>
constexpr int NFFT=1<<16, KMAX=16384, NU=3000; constexpr double ULO=0.08;
#include <complex>
struct Row{ double u,V1x,V1g,V2x,V2g,V3x,V3g,Mxx,Mxg,Mgg,P,T; };
static std::vector<Row> tab(NU);
static double uOf(int i){ return std::exp(std::log(ULO)*(1.0-double(i)/(NU-1))); }
static void fft(std::vector<std::complex<double>>& a){ const int n=int(a.size());
  for(int i=1,j=0;i<n;++i){ int bit=n>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit; if(i<j) std::swap(a[i],a[j]); }
  for(int len=2;len<=n;len<<=1){ const std::complex<double> w=std::polar(1.0,-2*M_PI/len);
    for(int i=0;i<n;i+=len){ std::complex<double> v=1; for(int j=0;j<len/2;++j){ auto x=a[i+j],y=a[i+j+len/2]*v; a[i+j]=x+y; a[i+j+len/2]=x-y; v*=w; } } } }
// Harmonics of zeta = xi + i eta (units a) sampled uniformly in mean anomaly: c_k = <zeta e^{-ikM}>;
// rho_{k,+} = |c_k|^2, rho_{k,-} = |c_{-k}|^2.  No special functions (libstdc++ Bessel fails near order 2000).
static void harmonics(double u,std::vector<double>& rp,std::vector<double>& rm){
  const double e=std::sqrt(std::max(0.0,1-u*u)); std::vector<std::complex<double>> z(NFFT);
  double E=0; for(int j=0;j<NFFT;++j){ const double M=2*M_PI*j/NFFT; if(j==0) E=0; else E=E; // Newton from previous
    double Ej=(j==0)?0.0:E; for(int it=0;it<60;++it){ const double f=Ej-e*std::sin(Ej)-M, d=1-e*std::cos(Ej); const double st=f/d; Ej-=st; if(std::abs(st)<1e-15) break; }
    E=Ej; z[j]={std::cos(E)-e,u*std::sin(E)}; }
  fft(z); rp.assign(KMAX+1,0); rm.assign(KMAX+1,0);
  for(int k=1;k<=KMAX;++k){ rp[k]=std::norm(z[k])/(double(NFFT)*NFFT); rm[k]=std::norm(z[NFFT-k])/(double(NFFT)*NFFT); } }
static void build(){
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<NU;++i){ const double u=uOf(i), du=1e-5*u; std::vector<double> p0,m0,p1,m1,p2,m2;
    harmonics(u,p0,m0); harmonics(std::max(ULO*0.5,u-du),p1,m1); harmonics(std::min(1.0,u+du),p2,m2);
    const double span=std::min(1.0,u+du)-std::max(ULO*0.5,u-du); Row r{}; r.u=u;
    for(int k=1;k<=KMAX;++k) for(int s=-1;s<=1;s+=2){
      const double rho=s>0?p0[k]:m0[k], drho=((s>0?p2[k]:m2[k])-(s>0?p1[k]:m1[k]))/span;
      const double k3=double(k)*k*k, w=k3*rho, w2=k3*(-5.0*k*rho+drho*(s-k*u));
      r.V1x+=w*k; r.V1g+=w*s; r.V2x+=w2*k; r.V2g+=w2*s; r.V3x+=w*s*k; r.V3g+=w*s*s;
      r.Mxx+=w*k*k; r.Mxg+=w*k*s; r.Mgg+=w*s*s; r.P+=k3*k*rho; r.T+=w*s; }
    tab[i]=r; } }
static Row at(double u){ const double t=(std::log(std::clamp(u,ULO,1.0))/std::log(ULO)); const double f=(1.0-t)*(NU-1);
  int i=std::min(NU-2,int(f)); const double w=f-i; const Row&a=tab[i],&b=tab[i+1]; Row r;
  #define L(m) r.m=a.m+(b.m-a.m)*w
  L(u);L(V1x);L(V1g);L(V2x);L(V2g);L(V3x);L(V3g);L(Mxx);L(Mxg);L(Mgg);L(P);L(T); return r; }
int main(int argc,char**argv){
  const double WALL=argc>1?std::atof(argv[1]):ULO; const bool skipChecks=argc>1;
  build();
  std::printf("# reflecting wall at u = %.3f\n",WALL);
  if(!skipChecks){
  std::printf("# H387a sum rules (truncated k <= %d)\n#  u      power/closed   torque/closed\n",KMAX);
  double worst=0;
  for(double u: {1.0,0.9,0.7,0.5,0.4,0.3,0.25,0.2,0.15,0.12,0.1,0.08}){ const Row r=at(u); const double e2=1-u*u;
    const double p=r.P/((1+e2/2)/std::pow(u,5)), q=r.T/(1/(u*u)); if(u>=0.25) worst=std::max({worst,std::abs(p-1),std::abs(q-1)});
    std::printf("  %.3f  %.8f     %.8f\n",u,p,q); }
  std::printf("# H387a worst for u >= 0.25: %.3e\n",worst); std::fflush(stdout);
  // H387b: 1D SED oscillator, R = J: dJ = (-J + 1/2) dt + sqrt(J dt) xi.
  { std::mt19937_64 rng(387); std::normal_distribution<double> N(0,1); double sum=0; const int P=4000;
    // Fixed-step ladder (a relative step rule stalls at the reflecting J = 0 boundary).
    for(double h: {1e-3,3e-4,1e-4}){ sum=0;
      for(int p=0;p<P;++p){ double J=1.0,t=0; while(t<20){ const double dt=std::min(h,20-t);
          J+=(-J+0.5)*dt+std::sqrt(J*dt)*N(rng); if(J<0) J=-J; t+=dt; } sum+=J; }
      std::printf("# H387b oscillator dt=%.0e <J> = %.4f (Boyer 0.5)\n",h,sum/P); std::fflush(stdout); } }
  }
  const double snaps[]={0,0.5,1,2,5,10,20,50,100,200};
  for(double g0: {0.5,1.0}){
    const int P=4000; std::vector<double> X(P,1.0),G(P,g0); std::vector<int> st(P,0); std::vector<double> umin(P,g0); // 0 alive, 1 collapse, 2 ionized
    std::mt19937_64 rng(1000+int(10*g0)); std::normal_distribution<double> N(0,1); double t=0;
    std::printf("\n# start (x,g) = (1, %.1f)\n#  t     alive  collapse  ionized   med x   med g   med u   <r>     <r^2>    <1/r>   <1/r^2>  frac g<0.75\n",g0);
    for(double ts: snaps){
      for(int p=0;p<P;++p){ double tp=t; double &x=X[p],&g=G[p]; long steps=0;
        while(st[p]==0&&tp<ts){ if(++steps>20000000L){ st[p]=3; break; } const Row r=at(g/x); const double x5=std::pow(x,-5);
          const double Ax=x5*(-r.V1x+0.5*(r.V2x/x+r.V3x/g)), Ag=x5*(-r.V1g+0.5*(r.V2g/x+r.V3g/g));
          const double Cxx=x5*r.Mxx,Cxg=x5*r.Mxg,Cgg=x5*r.Mgg;
          double dt=ts-tp; dt=std::min({dt,2e-3*x/std::abs(Ax),2e-3*g/std::abs(Ag),4e-4*x*x/Cxx,4e-4*g*g/Cgg});
          const double l11=std::sqrt(Cxx*dt), l21=Cxg*dt/l11, l22=std::sqrt(std::max(0.0,Cgg*dt-l21*l21));
          const double z1=N(rng),z2=N(rng); x+=Ax*dt+l11*z1; g+=Ag*dt+l21*z1+l22*z2; tp+=dt;
          if(x<=0.3){st[p]=1;break;} if(x>=20){st[p]=2;break;}
          if(g<WALL*x) g=2*WALL*x-g; umin[p]=std::min(umin[p],g/x); if(g>x) g=2*x-g; } }
      t=ts; std::vector<double> xs,gs,us; double m1=0,m2=0,m3=0,m4=0,low=0; int c1=0,c2=0,c3=0;
      for(int p=0;p<P;++p){ if(st[p]==1){++c1;continue;} if(st[p]==2){++c2;continue;} if(st[p]==3){++c3;continue;}
        const double x=X[p],g=G[p]; xs.push_back(x); gs.push_back(g); us.push_back(g/x);
        m1+=(3*x*x-g*g)/2; m2+=x*x*(5*x*x-3*g*g)/2; m3+=1/(x*x); m4+=1/(x*x*x*g); low+=(g<0.75); }
      const int n=int(xs.size()); auto med=[](std::vector<double> v)->double{ if(v.empty()) return NAN; std::nth_element(v.begin(),v.begin()+v.size()/2,v.end()); return v[v.size()/2]; };
      std::printf("  %-5g %5d  %5d    %5d     %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n",ts,n,c1,c2,med(xs),med(gs),med(us),
        n?m1/n:NAN,n?m2/n:NAN,n?m3/n:NAN,n?m4/n:NAN,n?low/n:NAN); if(c3) std::printf("    (stuck at step cap: %d)\n",c3); std::fflush(stdout); }
    { std::vector<double> ui,uc; for(int p=0;p<P;++p){ if(st[p]==2) ui.push_back(umin[p]); if(st[p]==1) uc.push_back(umin[p]); }
      std::sort(ui.begin(),ui.end()); std::sort(uc.begin(),uc.end());
      auto q=[](const std::vector<double>& v,double f){ return v.empty()?NAN:v[size_t(f*(v.size()-1))]; };
      std::printf("#  min u before ionization: 10/50/90 %% = %.3f / %.3f / %.3f;  before collapse: %.3f / %.3f / %.3f\n",
        q(ui,0.1),q(ui,0.5),q(ui,0.9),q(uc,0.1),q(uc,0.5),q(uc,0.9)); }
    // final histograms of survivors
    std::printf("#  survivors at t=200: histogram x (0.3..20, 20 bins log) and g\n#  x:");
    { std::vector<int> h(20,0); for(int p=0;p<P;++p) if(st[p]==0){ int b=int(20*std::log(X[p]/0.3)/std::log(20/0.3)); h[std::clamp(b,0,19)]++; } for(int v:h) std::printf(" %d",v); }
    std::printf("\n#  g (0..5, 20 bins):"); { std::vector<int> h(20,0); for(int p=0;p<P;++p) if(st[p]==0){ int b=int(G[p]/0.25); h[std::clamp(b,0,19)]++; } for(int v:h) std::printf(" %d",v); }
    std::printf("\n"); }
}
