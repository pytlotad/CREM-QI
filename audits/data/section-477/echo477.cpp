// Audit 477: as audit 473 (retarded / advanced / time-symmetric / consistent dipole echo between the two spins of Ps on a
// Kepler orbit), with the electric field of the moving moments (p = v x m/c^2) and any (n, L).
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
struct V{double x,y,z;};
static V operator+(V a,V b){return {a.x+b.x,a.y+b.y,a.z+b.z};} static V operator-(V a,V b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static V operator*(V a,double k){return {a.x*k,a.y*k,a.z*k};} static double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static V cross(V a,V b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};} static double nrm(V a){return std::sqrt(dot(a,a));}
const double hbar=1.054571817e-34,c=299792458.0,me=9.1093837015e-31,e0=1.602176634e-19,eps0=8.8541878128e-12,mu0=1.25663706212e-6,g=2.00231930436;
const double k=e0*e0/(4*M_PI*eps0), mu=me/2, aB=4*M_PI*eps0*hbar*hbar/(mu*e0*e0), muE=g*e0*hbar/(4*me);
double a, L, ecc, nmean, Torb, OmegaSO; bool useE=true;
const double aAnom=(g-2)/2;
static void relrv(double t,V& r,V& v){ double M=std::fmod(nmean*t,2*M_PI); if(M<0) M+=2*M_PI; double E=M+ecc*std::sin(M);
  for(int i=0;i<60;++i){ const double st=(E-ecc*std::sin(E)-M)/(1-ecc*std::cos(E)); E-=st; if(std::abs(st)<1e-15) break; }
  const double ce=std::cos(E), se=std::sin(E), q=std::sqrt(1-ecc*ecc), f=nmean*a/(1-ecc*ce);
  r={a*(ce-ecc),a*q*se,0}; v={-f*se,f*q*ce,0}; }
static V rel(double t){ V r,v; relrv(t,r,v); return r; }
// spin directions s1 (electron), s2 (positron); moments m1 = -muE s1, m2 = +muE s2; gamma1 = -g e/2m, gamma2 = +g e/2m
const double gam1=-g*e0/(2*me), gam2=g*e0/(2*me);
static V Bdip(V r,V m){ double rr=nrm(r); V n=r*(1/rr); return (n*(3*dot(n,m))-m)*(1e-7/(rr*rr*rr)); }
static V relv(double t){ V r,v; relrv(t,r,v); return v; }
static V Edip(V r,V p){ double rr=nrm(r); V n=r*(1/rr); return (n*(3*dot(n,p))-p)*(1/(4*M_PI*eps0*rr*rr*rr)); }
// particle 1 (electron) at +r/2, particle 2 (positron) at -r/2; r = r1 - r2
static void fields(double tpos,double tmom,V m1,V m2,V dm1,V dm2,double tau,V& B1,V& B2,V& E1,V& E2,V& v1,V& v2,bool consistent){
  V r=rel(tpos); V vr=relv(tpos);
  V m1s=m1+dm1*tau, m2s=m2+dm2*tau;  // moments at shifted time (first order)
  if(consistent){ m1s=m1; m2s=m2; }
  v1=vr*0.5; v2=vr*(-0.5);
  B1=Bdip(r,m2s); B2=Bdip(r*(-1.0),m1s);
  V p1=cross(v1,m1s)*(1/(c*c)), p2=cross(v2,m2s)*(1/(c*c));
  if(useE){ E1=Edip(r,p2); E2=Edip(r*(-1.0),p1); } else { E1={0,0,0}; E2={0,0,0}; } }
static V omega(double gam,V B,V E,V v){ const V zS{0,0,OmegaSO}; V Beff=B-cross(v,E)*((aAnom+0.5)/((aAnom+1)*c*c)); return Beff*(-gam)+zS; }
static void rates(double t,V s1,V s2,int mode,V& d1,V& d2){
  V m1=s1*(-muE), m2=s2*muE;
  // instantaneous rates for the moment derivatives
  V B1,B2,E1,E2,v1,v2; fields(t,t,m1,m2,{0,0,0},{0,0,0},0,B1,B2,E1,E2,v1,v2,true);
  V w1=omega(gam1,B1,E1,v1), w2=omega(gam2,B2,E2,v2); V dm1=cross(w1,m1), dm2=cross(w2,m2);
  if(mode!=0){
    const double tau=nrm(rel(t))/c;
    auto at=[&](double sh,V& b1,V& b2,V& e1,V& e2,V& u1,V& u2){ fields(t+sh*tau,t,m1,m2,dm1,dm2,sh*tau,b1,b2,e1,e2,u1,u2,false); };
    if(mode==1){ at(-1,B1,B2,E1,E2,v1,v2); }
    else if(mode==2){ at(+1,B1,B2,E1,E2,v1,v2); }
    else if(mode==3){ V a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6; at(-1,a1,a2,a3,a4,a5,a6); at(+1,b1,b2,b3,b4,b5,b6);
      B1=(a1+b1)*0.5; B2=(a2+b2)*0.5; E1=(a3+b3)*0.5; E2=(a4+b4)*0.5; v1=(a5+b5)*0.5; v2=(a6+b6)*0.5; }
    else if(mode==4){ // consistent to first order: retarded position + tau * retarded velocity, moment m(t-tau) + tau m' = m(t)
      V rr=rel(t-tau), vv=relv(t-tau); V rc=rr+vv*tau; V vr=relv(t);
      v1=vr*0.5; v2=vr*(-0.5); B1=Bdip(rc,m2); B2=Bdip(rc*(-1.0),m1);
      V p1=cross(v1,m1)*(1/(c*c)), p2=cross(v2,m2)*(1/(c*c));
      if(useE){ E1=Edip(rc,p2); E2=Edip(rc*(-1.0),p1); } else { E1={0,0,0}; E2={0,0,0}; } }
    w1=omega(gam1,B1,E1,v1); w2=omega(gam2,B2,E2,v2); }
  d1=cross(w1,s1); d2=cross(w2,s2); }
int main(int argc,char**argv){
  const int N=std::atoi(argv[1]); const double T=std::atof(argv[2]); const int nlev=std::atoi(argv[3]); const double Lh=std::atof(argv[4]); const int spo=std::atoi(argv[5]);
  useE=std::atoi(argv[6])!=0;
  a=nlev*nlev*aB; L=Lh*hbar; ecc=std::sqrt(1-L*L/(mu*k*a)); nmean=std::sqrt(k/(mu*a*a*a)); Torb=2*M_PI/nmean;
  { const double r3=1/(a*a*a*std::pow(1-ecc*ecc,1.5)); OmegaSO=(2*g-1)*k*L*r3/(2*me*me*c*c); } const double h=Torb/spo; const long steps=long(T/h);
  std::printf("# Ps n=%d L=%.2f hbar, E-part %d: a %.4e m, e %.4f, T_orb %.4e s, steps/orbit %d, steps %ld, Omega_SO %.4e rad/s\n",nlev,Lh,int(useE),a,ecc,Torb,spo,steps,OmegaSO);
  const char* nm[5]={"INST","RET","ADV","SYM","RETC"};
  std::vector<double> X0(N), XT[5];
  for(int v=0;v<5;++v) XT[v].resize(N);
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<N;++i){ std::mt19937_64 rng(477000+i); std::normal_distribution<double> G(0,1);
    auto unit=[&]{ V u{G(rng),G(rng),G(rng)}; return u*(1/nrm(u)); }; V s10=unit(), s20=unit(); X0[i]=dot(s10,s20);
    for(int v=0;v<5;++v){ V s1=s10,s2=s20; double t=0;
      for(long kk=0;kk<steps;++kk){ V a1,a2,b1,b2,c1,c2,d1,d2;
        rates(t,s1,s2,v,a1,a2); rates(t+0.5*h,s1+a1*(0.5*h),s2+a2*(0.5*h),v,b1,b2); rates(t+0.5*h,s1+b1*(0.5*h),s2+b2*(0.5*h),v,c1,c2); rates(t+h,s1+c1*h,s2+c2*h,v,d1,d2);
        s1=s1+(a1+b1*2+c1*2+d1)*(h/6); s2=s2+(a2+b2*2+c2*2+d2)*(h/6); s1=s1*(1/nrm(s1)); s2=s2*(1/nrm(s2)); t+=h; }
      XT[v][i]=dot(s1,s2); } }
  std::printf("# bin X0      n   <X_INST(T)-X0>   D_RET        D_ADV        D_SYM        D_RETC       (D = <X_v(T) - X_INST(T)>)\n");
  for(int b=0;b<10;++b){ double lo=-1+0.2*b, hi=lo+0.2; int n=0; double di=0,dr=0,da=0,ds=0,dc=0;
    for(int i=0;i<N;++i) if(X0[i]>=lo&&X0[i]<hi){ ++n; di+=XT[0][i]-X0[i]; dr+=XT[1][i]-XT[0][i]; da+=XT[2][i]-XT[0][i]; ds+=XT[3][i]-XT[0][i]; dc+=XT[4][i]-XT[0][i]; }
    if(n) std::printf("  [%+.1f,%+.1f) %4d  %+.4e  %+.4e  %+.4e  %+.4e  %+.4e\n",lo,hi,n,di/n,dr/n,da/n,ds/n,dc/n); }
  double mr=0,ma=0,ms=0,mc=0; for(int i=0;i<N;++i){ mr=std::max(mr,std::abs(XT[1][i]-XT[0][i])); ma=std::max(ma,std::abs(XT[2][i]-XT[0][i])); ms=std::max(ms,std::abs(XT[3][i]-XT[0][i])); mc=std::max(mc,std::abs(XT[4][i]-XT[0][i])); }
  std::printf("# max |X_v - X_INST|: RET %.3e ADV %.3e SYM %.3e RETC %.3e\n",mr,ma,ms,mc);
}
