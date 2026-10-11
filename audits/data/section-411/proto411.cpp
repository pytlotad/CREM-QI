// Audit 411: prototype -- does a self-generated wave of wave number kappa p/hbar (path memory) make the action
// J = (1/2pi) oint p.dx quantized?  2D, units hbar = mu = k = 1 (Bohr: r_n = n^2, J = n, T_n = 2 pi n^3).
// argv: eps kappa n0 [orbits=60] [gamma=2e-4] [tauOrbits=3]  (memory in LOCAL orbits, at most 4096 samples)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
struct V2 { double x=0,y=0; };
static V2 operator+(V2 a,V2 b){return {a.x+b.x,a.y+b.y};} static V2 operator-(V2 a,V2 b){return {a.x-b.x,a.y-b.y};}
static V2 operator*(V2 a,double s){return {a.x*s,a.y*s};} static double dot(V2 a,V2 b){return a.x*b.x+a.y*b.y;}
static double norm(V2 a){return std::sqrt(dot(a,a));}
struct Sample { double t, S; V2 x, p; };
static double eps, kappa, gam, tauM, tauOrbitsG;
static std::deque<Sample> mem;
// wave force on the particle at (x, t) with its own phase S
static V2 waveForce(V2 x,double t,double S,double Tloc){
  V2 F; if(eps==0.0) return F;
  const double tmin=0.5*Tloc; const double tau=tauOrbitsG*Tloc;
  for(const Sample& m: mem){ const double dt=t-m.t; if(dt<tmin) continue;
    const double w=std::exp(-dt/tau); const V2 d=x-m.x; const double pn=norm(m.p); if(pn<=0) continue;
    const double sigma=0.25*2*M_PI/(kappa*pn); const double d2=dot(d,d); if(d2>16*sigma*sigma) continue;
    const double G=std::exp(-0.5*d2/(sigma*sigma));
    const double phi=kappa*(m.S+dot(m.p,d)-S);
    // V = -eps w G cos(phi);  F = -grad V = eps w [ grad G cos phi - G sin phi kappa p ]
    const V2 gradG=d*(-G/(sigma*sigma));
    F=F+(gradG*std::cos(phi)-m.p*(kappa*G*std::sin(phi)))*(eps*w); }
  return F; }
int main(int argc,char**argv){
  eps=std::atof(argv[1]); kappa=std::atof(argv[2]); const double n0=std::atof(argv[3]);
  const double orbits=argc>4?std::atof(argv[4]):60; gam=argc>5?std::atof(argv[5]):2e-4; const double tauOrbits=argc>6?std::atof(argv[6]):3;
  const double r0=n0*n0; V2 x{r0,0}, v{0,1.0/n0}; double S=0, t=0;
  const double T0=2*M_PI*std::pow(n0,3); tauM=tauOrbits*T0; tauOrbitsG=tauOrbits; const double tEnd=orbits*T0;
  double nextSample=0;
  const auto deriv=[&](V2 xx,V2 vv,double tt,double SS,double Tloc,V2& a){ const double r=norm(xx);
    a=xx*(-1.0/(r*r*r))+vv*(-gam)+waveForce(xx,tt,SS,Tloc); };
  long step=0; double lastPrint=0;
  while(t<tEnd){
    const double r=norm(x), Tloc=2*M_PI*std::pow(r,1.5), dt=Tloc/2000;
    if(t>=nextSample){ mem.push_back({t,S,x,v}); nextSample=t+Tloc/256; while(!mem.empty()&&(t-mem.front().t>6*tauOrbits*Tloc||mem.size()>4096)) mem.pop_front(); }
    // RK4 on (x, v, S); S' = v.v
    V2 a1,a2,a3,a4; deriv(x,v,t,S,Tloc,a1);
    V2 x2=x+v*(dt/2), v2=v+a1*(dt/2); double S2=S+dot(v,v)*dt/2; deriv(x2,v2,t+dt/2,S2,Tloc,a2);
    V2 x3=x+v2*(dt/2), v3=v+a2*(dt/2); double S3=S+dot(v2,v2)*dt/2; deriv(x3,v3,t+dt/2,S3,Tloc,a3);
    V2 x4=x+v3*dt, v4=v+a3*dt; double S4=S+dot(v3,v3)*dt; deriv(x4,v4,t+dt,S4,Tloc,a4);
    x=x+(v+v2*2.0+v3*2.0+v4)*(dt/6); S=S+(dot(v,v)+2*dot(v2,v2)+2*dot(v3,v3)+dot(v4,v4))*dt/6;
    v=v+(a1+a2*2.0+a3*2.0+a4)*(dt/6); t+=dt; ++step;
    if(t-lastPrint>T0){ lastPrint=t; const double E=0.5*dot(v,v)-1/norm(x); const double L=x.x*v.y-x.y*v.x;
      std::fprintf(stderr,"t/T0=%.2f n=%.5f L=%.5f\n",t/T0,E<0?1/std::sqrt(-2*E):-1,L); }
  }
  const double E=0.5*dot(v,v)-1/norm(x); const double L=x.x*v.y-x.y*v.x; const double n=E<0?1/std::sqrt(-2*E):-1;
  std::printf("eps=%g kappa=%g n0=%.2f gamma=%g -> final n=%.5f (kappa n=%.5f) L=%.5f  steps=%ld\n",eps,kappa,n0,gam,n,kappa*n,L,step);
}
