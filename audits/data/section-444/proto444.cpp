// Audit 444: three internal action clocks in 3D (radial, polar, azimuthal), Maslov jumps at radial and polar turning points.
// hbar = mu = k = 1.  argv: J0 e0 incl_deg orbits eps kappa_r kappa_th kappa_phi gamma
#include <cmath>
#include <cstdio>
#include <cstdlib>
struct S{double x,y,z,vx,vy,vz,cr,ct,cp;};
double eps,kr,kt,kp,gam; bool maslov=true;
struct Geo{double r,th,ph,rd,thd,phd; double rh[3],th_h[3],ph_h[3]; double Lx,Ly,Lz,phn;};
Geo geo(const S& s){ Geo g; g.r=std::sqrt(s.x*s.x+s.y*s.y+s.z*s.z); g.th=std::acos(s.z/g.r); g.ph=std::atan2(s.y,s.x);
  const double st=std::sin(g.th), ct=std::cos(g.th), sp=std::sin(g.ph), cp=std::cos(g.ph);
  g.rh[0]=st*cp; g.rh[1]=st*sp; g.rh[2]=ct; g.th_h[0]=ct*cp; g.th_h[1]=ct*sp; g.th_h[2]=-st; g.ph_h[0]=-sp; g.ph_h[1]=cp; g.ph_h[2]=0;
  g.rd=s.vx*g.rh[0]+s.vy*g.rh[1]+s.vz*g.rh[2]; g.thd=(s.vx*g.th_h[0]+s.vy*g.th_h[1]+s.vz*g.th_h[2]); g.phd=(s.vx*g.ph_h[0]+s.vy*g.ph_h[1]);
  g.Lx=s.y*s.vz-s.z*s.vy; g.Ly=s.z*s.vx-s.x*s.vz; g.Lz=s.x*s.vy-s.y*s.vx;
  g.phn=std::atan2(g.Lx,-g.Ly);   // node direction n = z x L = (-Ly, Lx, 0)
  return g; }
S der(const S& s){ const Geo g=geo(s);
  const double gr=std::exp(-g.r), dgr=-gr;                          // radial coupling shape and d/dr
  const double ht=std::cos(g.th), dht=-std::sin(g.th);              // polar coupling shape and d/dth
  const double u=g.ph-g.phn, hp=std::exp(std::cos(u)), dhp=-std::sin(u)*hp;   // azimuthal coupling and d/dphi
  const double Fr=-1/(g.r*g.r)+eps*dgr*std::cos(s.cr);
  static const double et=std::getenv("EPS_TH")?std::atof(std::getenv("EPS_TH")):eps, ep=std::getenv("EPS_PHI")?std::atof(std::getenv("EPS_PHI")):eps;
  const double Qt=et*dht*std::cos(s.ct), Qp=ep*dhp*std::cos(s.cp);   // generalized forces -dU/dth, -dU/dphi
  const double Ft=Qt/g.r, Fp=Qp/(g.r*std::sin(g.th));
  double a[3]; for(int i=0;i<3;++i) a[i]=Fr*g.rh[i]+Ft*g.th_h[i]+Fp*g.ph_h[i];
  a[0]-=gam*s.vx; a[1]-=gam*s.vy; a[2]-=gam*s.vz;
  const double cr=g.rd*g.rd-kr*eps*(dgr*g.rd)*std::cos(s.cr);     // radial clock pulled by radial coupling power
  const double ct=g.thd*g.thd-kt*Qt;                                // polar clock pulled by Q_th
  const double cp=g.phd*g.phd-kp*Qp;                                // azimuthal clock pulled by Q_phi
  return {s.vx,s.vy,s.vz,a[0],a[1],a[2],cr,ct,cp}; }
S add(const S& a,const S& b,double h){ return {a.x+h*b.x,a.y+h*b.y,a.z+h*b.z,a.vx+h*b.vx,a.vy+h*b.vy,a.vz+h*b.vz,a.cr+h*b.cr,a.ct+h*b.ct,a.cp+h*b.cp}; }
int main(int argc,char** argv){ setvbuf(stdout,nullptr,_IOLBF,0);
  const double J0=std::atof(argv[1]), e0=std::atof(argv[2]), inc=std::atof(argv[3])*M_PI/180; const int orbits=std::atoi(argv[4]);
  eps=std::atof(argv[5]); kr=std::atof(argv[6]); kt=std::atof(argv[7]); kp=std::atof(argv[8]); gam=std::atof(argv[9]);
  if(std::getenv("NO_MASLOV")) maslov=false;
  const double a=J0*J0, rp=a*(1-e0), vp=std::sqrt((1+e0)/(a*(1-e0)));
  S s{rp,0,0,0,vp*std::cos(inc),vp*std::sin(inc),0,0,0};
  int count=0, nb=0; double sJ=0,sL=0,sLz=0,sR=0,sT=0,sP=0, pr=0,pt=0,pp=0;
  std::printf("# orbit <J> <L> <Lz> <J-L> <L-|Lz|> clk_r clk_th clk_phi\n");
  while(count<orbits){ const Geo g=geo(s); const double E=0.5*(s.vx*s.vx+s.vy*s.vy+s.vz*s.vz)-1/g.r; if(!(E<0)){ std::printf("# unbound\n"); break; }
    const double J=1/std::sqrt(-2*E), dt=2*M_PI*J*J*J/4000;
    const S k1=der(s), k2=der(add(s,k1,dt/2)), k3=der(add(s,k2,dt/2)), k4=der(add(s,k3,dt));
    S n=add(s,k1,dt/6); n=add(n,k2,dt/3); n=add(n,k3,dt/3); n=add(n,k4,dt/6);
    const Geo g1=geo(n);
    if(maslov){ if(g.rd*g1.rd<0) n.cr-=M_PI/2; if(g.thd*g1.thd<0) n.ct-=M_PI/2; }
    const bool peri=(g.rd<0&&g1.rd>=0);
    s=n;
    if(peri){ ++count; const Geo h=geo(s); const double EE=0.5*(s.vx*s.vx+s.vy*s.vy+s.vz*s.vz)-1/h.r, JJ=1/std::sqrt(-2*EE);
      const double L=std::sqrt(h.Lx*h.Lx+h.Ly*h.Ly+h.Lz*h.Lz);
      sJ+=JJ; sL+=L; sLz+=h.Lz; sR+=(s.cr-pr)/(2*M_PI); sT+=(s.ct-pt)/(2*M_PI); sP+=(s.cp-pp)/(2*M_PI); pr=s.cr; pt=s.ct; pp=s.cp; ++nb;
      if(count%(std::getenv("EVERY")?std::atoi(std::getenv("EVERY")):1000)==0){ std::printf("%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",count,sJ/nb,sL/nb,sLz/nb,(sJ-sL)/nb,(sL-std::fabs(sLz))/nb,sR/nb,sT/nb,sP/nb); sJ=sL=sLz=sR=sT=sP=0; nb=0; } } } }
