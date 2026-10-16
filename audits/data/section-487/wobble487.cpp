// Audit 487 (from 477): within-orbit wobble of the spins, instantaneous fields vs orbit-averaged rates.
// Audit 477 origin: as audit 473 (retarded / advanced / time-symmetric / consistent dipole echo between the two spins of Ps on a
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
double gSOinst=0; // instantaneous spin-orbit rate along z, set per time in rates()
static V omega(double gam,V B,V E,V v){ const V zS{0,0,gSOinst}; V Beff=B-cross(v,E)*((aAnom+0.5)/((aAnom+1)*c*c)); return Beff*(-gam)+zS; }
static void rates(double t,V s1,V s2,int mode,V& d1,V& d2){
  { V rr,vv; relrv(t,rr,vv); const double r=nrm(rr); gSOinst=(2*g-1)*k*(mu*nrm(cross(rr,vv)))/(2*me*me*c*c*r*r*r); }
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

static void rk4(double t,double h,V& s1,V& s2){ V a1,a2,b1,b2,c1,c2,d1,d2;
  rates(t,s1,s2,0,a1,a2); rates(t+0.5*h,s1+a1*(0.5*h),s2+a2*(0.5*h),0,b1,b2); rates(t+0.5*h,s1+b1*(0.5*h),s2+b2*(0.5*h),0,c1,c2); rates(t+h,s1+c1*h,s2+c2*h,0,d1,d2);
  s1=s1+(a1+b1*2+c1*2+d1)*(h/6); s2=s2+(a2+b2*2+c2*2+d2)*(h/6); s1=s1*(1/nrm(s1)); s2=s2*(1/nrm(s2)); }
static V rot(V v,V w,double dt){ double th=nrm(w)*dt; if(th<1e-300) return v; V ax=w*(1/nrm(w)); return v*std::cos(th)+cross(ax,v)*std::sin(th)+ax*(dot(ax,v)*(1-std::cos(th))); }
int main(int argc,char**argv){
  const int nlev=std::atoi(argv[1]); const double Lh=std::atof(argv[2]); const int spo=std::atoi(argv[3]); const int pair=std::atoi(argv[4]); // pair 0 Ps, 1 H (electron spin only; proton fixed)
  a=nlev*nlev*aB; L=Lh*hbar; ecc=std::sqrt(1-L*L/(mu*k*a)); nmean=std::sqrt(k/(mu*a*a*a)); Torb=2*M_PI/nmean; useE=true;
  if(pair==1){ std::printf("# hydrogen variant not implemented in this reduced program\n"); return 0; }
  const double h=Torb/spo; const int K=3*spo;
  std::printf("# Ps n=%d L=%.2f hbar: a %.4e m, e %.4f, omega_orb %.4e rad/s, (3/4) r_e/a = %.3e\n",nlev,Lh,a,ecc,nmean,0.75*2.8179403262e-15/a);
  std::mt19937_64 rng(487); std::normal_distribution<double> G(0,1);
  for(int cfg=0;cfg<7;++cfg){
    V s10,s20; if(cfg==0){ double th=M_PI/3; s10={std::sin(th),0,std::cos(th)}; s20=s10; }
    else { s10={G(rng),G(rng),G(rng)}; s10=s10*(1/nrm(s10)); s20={G(rng),G(rng),G(rng)}; s20=s20*(1/nrm(s20)); }
    // orbit-averaged rates at the start configuration (secular): average d/dt over one orbit with frozen spins -> angular velocities
    V w1{0,0,0},w2{0,0,0};
    for(int i=0;i<spo;++i){ const double t=(i+0.5)*h; V d1,d2; rates(t,s10,s20,0,d1,d2);
      // recover omega from d = omega x s (component perpendicular to s): use the instantaneous fields directly
      V m1=s10*(-muE), m2=s20*muE, B1,B2,E1,E2,v1,v2; fields(t,t,m1,m2,{0,0,0},{0,0,0},0,B1,B2,E1,E2,v1,v2,true);
      w1=w1+omega(gam1,B1,E1,v1)*(1.0/spo); w2=w2+omega(gam2,B2,E2,v2)*(1.0/spo); }
    V s1=s10,s2=s20, q1=s10,q2=s20; double maxang=0, rms=0, maxtilt=0, maxaz=0; double S12f=0,S12s=0;
    std::vector<double> dev(K); 
    for(int kk=0;kk<K;++kk){ const double t=kk*h; rk4(t,h,s1,s2); q1=rot(q1,w1,h); q2=rot(q2,w2,h);
      const double ang=std::acos(std::min(1.0,dot(s1,q1))); maxang=std::max(maxang,ang); rms+=ang*ang/K;
      // decomposition about L (z): tilt = change of the polar angle, azimuth = rest
      const double tilt=std::abs(std::acos(std::max(-1.0,std::min(1.0,s1.z)))-std::acos(std::max(-1.0,std::min(1.0,q1.z))));
      maxtilt=std::max(maxtilt,tilt); maxaz=std::max(maxaz,std::sqrt(std::max(0.0,ang*ang-tilt*tilt)));
      dev[kk]=std::atan2(cross(q1,s1).z,dot(q1,s1));
      V r=rel(t+h); V n=r*(1/nrm(r)); S12f+=(3*dot(s1,n)*dot(s2,n)-dot(s1,s2))/K; S12s+=(3*dot(q1,n)*dot(q2,n)-dot(q1,q2))/K; }
    // harmonic content of the azimuthal deviation (signed rotation about z)
    std::printf("  cfg %d: max angle %.3e rad, rms %.3e; max tilt %.3e, max azimuthal %.3e; <S12> full %+.6f secular %+.6f diff %+.3e\n     harmonics |c_k| k=1..8:",
      cfg,maxang,std::sqrt(rms),maxtilt,maxaz,S12f,S12s,S12f-S12s);
    for(int kh=1;kh<=8;++kh){ double re=0,im=0; for(int kk=0;kk<K;++kk){ double ph=2*M_PI*kh*kk/double(spo); re+=dev[kk]*std::cos(ph)/K; im+=dev[kk]*std::sin(ph)/K; } std::printf(" %.2e",2*std::sqrt(re*re+im*im)); }
    std::printf("\n");
  }
}
