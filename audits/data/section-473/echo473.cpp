// Audit 473: retarded / advanced / time-symmetric dipole echo between the two spins of Ps 2p on a Kepler orbit.
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
const double a=4*aB, L=1.5*hbar, ecc=std::sqrt(1-L*L/(mu*k*a)), nmean=std::sqrt(k/(mu*a*a*a)), Torb=2*M_PI/nmean;
const double OmegaSO=2.296758e10;            // common spin-orbit precession (audit 467), along z
static V rel(double t){ double M=std::fmod(nmean*t,2*M_PI); double E=M; for(int i=0;i<40;++i) E-=(E-ecc*std::sin(E)-M)/(1-ecc*std::cos(E));
  return {a*(std::cos(E)-ecc),a*std::sqrt(1-ecc*ecc)*std::sin(E),0}; }
// spin directions s1 (electron), s2 (positron); moments m1 = -muE s1, m2 = +muE s2; gamma1 = -g e/2m, gamma2 = +g e/2m
const double gam1=-g*e0/(2*me), gam2=g*e0/(2*me);
static V Bdip(V r,V m){ double rr=nrm(r); V n=r*(1/rr); return (n*(3*dot(n,m))-m)*(1e-7/(rr*rr*rr)); }
static void rates(double t,V s1,V s2,int mode,V& d1,V& d2){
  const V zS{0,0,OmegaSO}; V m1=s1*(-muE), m2=s2*muE;
  auto field=[&](double sh,V& B1,V& B2){ // partner evaluated at t + sh*r/c (sh = -1 retarded, +1 advanced, 0 instantaneous)
    V r=rel(t); double tau=nrm(r)/c*sh; V rs=rel(t+tau);
    // partner moment at shifted time, first order: m(t+tau) = m + tau * dm/dt (instantaneous rates)
    V B1i=Bdip(r,m2), B2i=Bdip(r*(-1.0),m1);
    V w1=B1i*(-gam1)+zS, w2=B2i*(-gam2)+zS; V dm1=cross(w1,m1), dm2=cross(w2,m2);
    B1=Bdip(rs,m2+dm2*tau); B2=Bdip(rs*(-1.0),m1+dm1*tau); };
  V B1,B2;
  if(mode==4){ // RETC: consistent to first order -- moment m(t-tau) + tau m' and position r(t-tau) + tau v(t-tau)
    V r=rel(t); double tau=nrm(r)/c; V rr=rel(t-tau); const double dt=1e-21; V vv=(rel(t-tau+dt)-rel(t-tau-dt))*(1/(2*dt));
    V rc=rr+vv*tau; V B1i=Bdip(r,m2), B2i=Bdip(r*(-1.0),m1);
    V w1=B1i*(-gam1)+zS, w2=B2i*(-gam2)+zS; V dm1=cross(w1,m1), dm2=cross(w2,m2);
    V m2r=m2-dm2*tau, m1r=m1-dm1*tau;          // moments at the retarded time (first order)
    B1=Bdip(rc,m2r+dm2*tau); B2=Bdip(rc*(-1.0),m1r+dm1*tau);  // + induction term restores m(t)
  } else if(mode==0) field(0,B1,B2); else if(mode==1) field(-1,B1,B2); else if(mode==2) field(1,B1,B2);
  else { V a1,a2,b1,b2; field(-1,a1,a2); field(1,b1,b2); B1=(a1+b1)*0.5; B2=(a2+b2)*0.5; }
  V w1=B1*(-gam1)+zS, w2=B2*(-gam2)+zS; d1=cross(w1,s1); d2=cross(w2,s2); }
int main(int argc,char**argv){
  const int N=std::atoi(argv[1]); const double T=std::atof(argv[2]); const int spo=64; const double h=Torb/spo; const long steps=long(T/h);
  std::printf("# Ps n=2 L=1.5 hbar: a %.4e m, e %.4f, T_orb %.4e s, steps %ld, tau/T_orb at a: %.3e\n",a,ecc,Torb,steps,a/c/Torb);
  const char* nm[5]={"INST","RET","ADV","SYM","RETC"};
  std::vector<double> X0(N), XT[5];
  for(int v=0;v<5;++v) XT[v].resize(N);
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<N;++i){ std::mt19937_64 rng(473000+i); std::normal_distribution<double> G(0,1);
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
