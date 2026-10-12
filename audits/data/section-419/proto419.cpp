// Audit 419: 3D prototype of the path-memory wave (k = kappa p/hbar).  Units hbar = mu = k = 1.
// Forces: Coulomb -r/r^3; optional perihelion precession from V = -beta/r^2 (F = -2 beta r/r^4); optional plane
// precession from a uniform field B z (F = v x B, Larmor rate B/2); wave force (direct coupling, phase phic) from
// samples deposited every T/256 with a 3D gaussian envelope sigma = s lambda/4, memory exp(-dt/tau), tau in LOCAL
// orbits, current loop excluded (dt < T/2).
// Output per orbit (stderr): E, n = 1/sqrt(-2E), |L|, L_z, inclination, J_r = 1/sqrt(-2E) - sqrt(L^2 - 2 beta).
// argv: eps n0 e0 incl_deg beta B tauOrbits s orbits [phic=0] [maxSamples=16384]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
struct V3 { double x=0,y=0,z=0; };
static V3 operator+(V3 a,V3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};} static V3 operator-(V3 a,V3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static V3 operator*(V3 a,double s){return {a.x*s,a.y*s,a.z*s};} static double dot(V3 a,V3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static V3 cross(V3 a,V3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};} static double norm(V3 a){return std::sqrt(dot(a,a));}
struct Sample { double t, S; V3 x, p; };
static double eps, beta, Bz, tauO, sScale, phic; static std::size_t maxSamples;
static std::deque<Sample> mem;
static V3 waveForce(V3 x,double t,double S,double Tloc){
  V3 F; if(eps==0.0) return F; const double tmin=0.5*Tloc, tau=tauO*Tloc;
  for(const Sample& m: mem){ const double dt=t-m.t; if(dt<tmin) continue; const double pn=norm(m.p); if(pn<=0) continue;
    const double w=std::exp(-dt/tau); const V3 d=x-m.x; const double sigma=sScale*0.25*2*M_PI/pn; const double d2=dot(d,d);
    if(d2>16*sigma*sigma) continue; const double G=std::exp(-0.5*d2/(sigma*sigma));
    const double phi=m.S+dot(m.p,d)-S;
    const V3 gradG=d*(-G/(sigma*sigma));
    F=F+(gradG*std::cos(phi-phic)-m.p*(G*std::sin(phi-phic)))*(eps*w); }
  return F; }
int main(int argc,char**argv){
  if(argc<10){ std::fprintf(stderr,"usage: eps n0 e0 incl_deg beta B tauOrbits s orbits [phic] [maxSamples]\n"); return 1; }
  eps=std::atof(argv[1]); const double n0=std::atof(argv[2]), e0=std::atof(argv[3]), incl=std::atof(argv[4])*M_PI/180;
  beta=std::atof(argv[5]); Bz=std::atof(argv[6]); tauO=std::atof(argv[7]); sScale=std::atof(argv[8]); const double orbits=std::atof(argv[9]);
  phic=argc>10?std::atof(argv[10]):0.0; maxSamples=argc>11?std::atol(argv[11]):16384;
  // Kepler ellipse (a = n0^2, e0) from apoapsis, in a plane tilted by incl about the x axis
  const double a0=n0*n0, ra=a0*(1+e0), va=std::sqrt((1-e0)/(a0*(1+e0)));
  V3 x{ra,0,0}, v{0,va*std::cos(incl),va*std::sin(incl)}; double S=0, t=0;
  const double T0=2*M_PI*std::pow(n0,3), tEnd=orbits*T0; double nextSample=0, lastPrint=0;
  const auto acc=[&](V3 xx,V3 vv,double tt,double SS,double Tloc){ const double r=norm(xx);
    V3 a=xx*(-1.0/(r*r*r)); if(beta!=0) a=a+xx*(-2*beta/(r*r*r*r)); if(Bz!=0) a=a+cross(vv,V3{0,0,Bz});
    return a+waveForce(xx,tt,SS,Tloc); };
  const auto report=[&](FILE* f,const char* tag){ const double r=norm(x); const double E=0.5*dot(v,v)-1/r-beta/(r*r); const V3 L=cross(x,v);
    const double Ln=norm(L), n=E<0?1/std::sqrt(-2*E):-1, Leff2=Ln*Ln-2*beta;
    std::fprintf(f,"%s t/T0=%.2f n=%.5f L=%.5f Lz=%.5f incl=%.3f Jr=%.5f\n",tag,t/T0,n,Ln,L.z,std::acos(std::max(-1.0,std::min(1.0,L.z/Ln)))*180/M_PI,
      (E<0&&Leff2>0)?n-std::sqrt(Leff2):-1); };
  bool collapsed=false;
  while(t<tEnd){
    const double r=norm(x); { const double E=0.5*dot(v,v)-1/r-beta/(r*r); if(E<0&&1/std::sqrt(-2*E)<0.8){collapsed=true;break;} }
    const double Tloc=2*M_PI*std::pow(r,1.5), dt=Tloc/2000;
    if(t>=nextSample){ mem.push_back({t,S,x,v}); nextSample=t+Tloc/256;
      while(!mem.empty()&&(t-mem.front().t>6*tauO*Tloc||mem.size()>maxSamples)) mem.pop_front(); }
    const V3 a1=acc(x,v,t,S,Tloc);
    const V3 x2=x+v*(dt/2), v2=v+a1*(dt/2); const double S2=S+dot(v,v)*dt/2; const V3 a2=acc(x2,v2,t+dt/2,S2,Tloc);
    const V3 x3=x+v2*(dt/2), v3=v+a2*(dt/2); const double S3=S+dot(v2,v2)*dt/2; const V3 a3=acc(x3,v3,t+dt/2,S3,Tloc);
    const V3 x4=x+v3*dt, v4=v+a3*dt; const double S4=S+dot(v3,v3)*dt; const V3 a4=acc(x4,v4,t+dt,S4,Tloc);
    x=x+(v+v2*2.0+v3*2.0+v4)*(dt/6); S=S+(dot(v,v)+2*dot(v2,v2)+2*dot(v3,v3)+dot(v4,v4))*dt/6; v=v+(a1+a2*2.0+a3*2.0+a4)*(dt/6); t+=dt;
    if(t-lastPrint>T0){ lastPrint=t; report(stderr,"orbit"); } }
  std::printf("eps=%g n0=%.2f e0=%.2f incl=%.1f beta=%g B=%g tau=%g s=%g:%s ",eps,n0,e0,incl*180/M_PI,beta,Bz,tauO,sScale,collapsed?" COLLAPSED":"");
  report(stdout,"final"); }
