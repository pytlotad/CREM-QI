// Audit 497 (from 488): two classical spins, contact A S1.S2 plus the orbit-averaged dipole tensor B[-S1.S2/2 + 3/2 (S1.n)(S2.n)]
// (orbit normal n at theta_n from the field axis z), in a slowly decreasing field -- separatrix crossing, final |S|.
// args: l, (A or B) T, w2/w1, N, A, B, theta_n [deg].  With A = 1, B = 0 identical to sep488.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include <algorithm>
struct V{double x,y,z;};
static V operator+(V a,V b){return {a.x+b.x,a.y+b.y,a.z+b.z};} static V operator*(V a,double k){return {a.x*k,a.y*k,a.z*k};}
static V cross(V a,V b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};} static double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static double nrm(V a){return std::sqrt(dot(a,a));}
int main(int argc,char**argv){
  const double l=std::atof(argv[1]), AT=std::atof(argv[2]), ratio=std::atof(argv[3]); const int N=std::atoi(argv[4]);   // ratio = w2/w1 (-1 Ps)
  const double A=std::atof(argv[5]), B=std::atof(argv[6]), thn=std::atof(argv[7])*M_PI/180; const V n{std::sin(thn),0,std::cos(thn)};
  const double K=std::max(A,B), w0=50*K, T=AT/K, dt=std::min(0.002, 0.02/(w0+2*l*(A+2*B))); const long steps=long(1.25*T/dt);
  std::mt19937_64 rng(488); std::normal_distribution<double> G(0,1);
  std::vector<double> Sf(N), Sa(N); const double checks[4]={40*K,20*K,5*K,1*K}; std::vector<int> rot(4,0);
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<N;++i){
    std::mt19937_64 r(488000+i); std::normal_distribution<double> g(0,1);
    V s1{g(r),g(r),g(r)}, s2{g(r),g(r),g(r)}; s1=s1*(l/nrm(s1)); s2=s2*(l/nrm(s2));
    double t=0; int ci=0; double Ssum=0; long Scount=0; double phiPrev=0, wind=0; bool tracking=false; double winStart=0;
    for(long k=0;k<steps;++k){
      const double w1=(t<T)?w0*(1-t/T):0.0, w2=ratio*w1;
      auto f=[&](V a,V b,V& da,V& db){ V O1=b*(A-0.5*B)+n*(1.5*B*dot(b,n))+V{0,0,w1}, O2=a*(A-0.5*B)+n*(1.5*B*dot(a,n))+V{0,0,w2}; da=cross(O1,a); db=cross(O2,b); };
      V a1,b1,a2,b2,a3,b3,a4,b4; f(s1,s2,a1,b1); f(s1+a1*(0.5*dt),s2+b1*(0.5*dt),a2,b2); f(s1+a2*(0.5*dt),s2+b2*(0.5*dt),a3,b3); f(s1+a3*dt,s2+b3*dt,a4,b4);
      s1=s1+(a1+a2*2+a3*2+a4)*(dt/6); s2=s2+(b1+b2*2+b3*2+b4)*(dt/6); s1=s1*(l/nrm(s1)); s2=s2*(l/nrm(s2)); t+=dt;
      if(t>T){ Ssum+=nrm(s1+s2); ++Scount; }
      // classify libration/rotation of the relative azimuth in a window around selected field values
      const double dphi=std::atan2(s2.y,s2.x)-std::atan2(s1.y,s1.x);
      if(ci<4){ const double target=T*(1-checks[ci]/w0);
        if(!tracking && t>=target){ tracking=true; wind=0; phiPrev=dphi; winStart=t; }
        else if(tracking){ double d=dphi-phiPrev; while(d>M_PI) d-=2*M_PI; while(d<-M_PI) d+=2*M_PI; wind+=d; phiPrev=dphi;
          if(t-winStart>20.0/K){ if(std::abs(wind)>2*M_PI){
              #pragma omp atomic
              rot[ci]++; }
            tracking=false; ++ci; } } }
    }
    Sf[i]=nrm(s1+s2);
    Sa[i]=Ssum/std::max(1L,Scount);
  }
  std::printf("l %.3f, A %.1f, B %.1f, theta_n %.0f, K T %.0f, w2/w1 %.5f, N %d: fraction rotating (relative azimuth winds > 2 pi in 20/K) at w1/K =",l,A,B,thn*180/M_PI,AT,ratio,N);
  for(int c=0;c<4;++c) std::printf(" %.0f: %.3f",checks[c]/K,double(rot[c])/N);
  std::printf("\n  final |S|/l histogram (20 bins 0..2):");
  std::vector<int> h(20,0); for(double v: Sf) h[std::min(19,int(v/l/0.1))]++; for(int v: h) std::printf(" %d",v);
  int below=0; for(double v: Sf) below+=(v<1.0);   // rule-Q basin boundary |S| = hbar
  const int mx=*std::max_element(h.begin(),h.end());
  std::printf("\n  largest bin %.3f of the ensemble; fraction with |S| < hbar %.4f\n",double(mx)/N,double(below)/N);
  std::vector<int> ha(20,0); for(double v: Sa) ha[std::min(19,int(v/l/0.1))]++; int ba=0; for(double v: Sa) ba+=(v<1.0);
  std::printf("  |S| averaged over the last T/4 (no field): histogram"); for(int v: ha) std::printf(" %d",v);
  std::printf("\n  largest bin %.3f; fraction with <|S|> < hbar %.4f\n",double(*std::max_element(ha.begin(),ha.end()))/N,double(ba)/N);
}
