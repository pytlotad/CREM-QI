// Audit 441: two internal action clocks (radial p_r rdot, angular L phidot), radial coupling g(r) and torque coupling h(f).
// 2D relative coordinate, hbar = mu = k = 1.  argv: J0 e0 orbits eps_r eps_phi kappa_r kappa_phi gamma
#include <cmath>
#include <cstdio>
#include <cstdlib>
struct S{double x,y,vx,vy,tr,tp;};
double er,ep,kr,kp,gam;
static double trueAnomaly(const S& s,double& L){ const double r=std::hypot(s.x,s.y); L=s.x*s.vy-s.y*s.vx;
  const double Ax=s.vy*L-s.x/r, Ay=-s.vx*L-s.y/r; return std::atan2(s.y,s.x)-std::atan2(Ay,Ax); }
S der(const S& s){ const double r=std::hypot(s.x,s.y); double L; const double f=trueAnomaly(s,L);
  const double g=std::exp(-r), gp=-g; static const int hmode=std::getenv("HMODE")?std::atoi(std::getenv("HMODE")):0; const double h=hmode==2?std::cos(2*f):std::exp(std::cos(f)), hp=hmode==2?-2*std::sin(2*f):-std::sin(f)*h;
  const double rdot=(s.x*s.vx+s.y*s.vy)/r;
  const double Fr=-1/(r*r)+er*gp*std::cos(s.tr);                 // Coulomb + radial coupling force
  const double tau=ep*hp*std::cos(s.tp), Ft=tau/r;                // torque coupling -> tangential force
  const double ax=Fr*s.x/r+Ft*(-s.y/r)-gam*s.vx, ay=Fr*s.y/r+Ft*(s.x/r)-gam*s.vy;
  const double tr=rdot*rdot-kr*er*(gp*rdot)*std::cos(s.tr);       // radial clock: p_r rdot, pulled by radial coupling power
  const double tp=L*L/(r*r)-kp*ep*hp*std::cos(s.tp);              // angular clock: L phidot, pulled by the torque
  return {s.vx,s.vy,ax,ay,tr,tp}; }
int main(int argc,char** argv){ setvbuf(stdout,nullptr,_IOLBF,0);
  const double J0=std::atof(argv[1]), e0=std::atof(argv[2]); const int orbits=std::atoi(argv[3]);
  er=std::atof(argv[4]); ep=std::atof(argv[5]); kr=std::atof(argv[6]); kp=std::atof(argv[7]); gam=std::atof(argv[8]);
  const double a=J0*J0, rp=a*(1-e0), vp=std::sqrt((1+e0)/(a*(1-e0)));
  const double th0r=argc>9?std::atof(argv[9]):0, th0p=argc>10?std::atof(argv[10]):0;
  S s{rp,0,0,vp,th0r,th0p}; double angleAcc=0, trPrev=s.tr, tpPrev=s.tp; int count=0; double sJ=0,sL=0,sR=0,sP=0; int nb=0;
  std::printf("# orbit <J> <L> <J_r> <radial clock/2pi> <angular clock/2pi> e\n");
  auto add=[](const S& a,const S& b,double h){ return S{a.x+h*b.x,a.y+h*b.y,a.vx+h*b.vx,a.vy+h*b.vy,a.tr+h*b.tr,a.tp+h*b.tp}; };
  while(count<orbits){ const double r=std::hypot(s.x,s.y), E=0.5*(s.vx*s.vx+s.vy*s.vy)-1/r; if(!(E<0)){ std::printf("# unbound\n"); break; }
    const double J=1/std::sqrt(-2*E), dt=2*M_PI*J*J*J/4000;
    const S k1=der(s), k2=der(add(s,k1,dt/2)), k3=der(add(s,k2,dt/2)), k4=der(add(s,k3,dt));
    S n=add(s,k1,dt/6); n=add(n,k2,dt/3); n=add(n,k3,dt/3); n=add(n,k4,dt/6);
    double d=std::atan2(n.y,n.x)-std::atan2(s.y,s.x); if(d>M_PI) d-=2*M_PI; if(d<-M_PI) d+=2*M_PI; angleAcc+=d; s=n;
    if(angleAcc>=2*M_PI){ angleAcc-=2*M_PI; ++count;
      const double rr=std::hypot(s.x,s.y), EE=0.5*(s.vx*s.vx+s.vy*s.vy)-1/rr, JJ=1/std::sqrt(-2*EE), LL=s.x*s.vy-s.y*s.vx;
      sJ+=JJ; sL+=LL; sR+=(s.tr-trPrev)/(2*M_PI); sP+=(s.tp-tpPrev)/(2*M_PI); ++nb; trPrev=s.tr; tpPrev=s.tp;
      if(count%1000==0){ const double ecc=std::sqrt(std::max(0.0,1-LL*LL/(JJ*JJ)));
        std::printf("%d %.5f %.5f %.5f %.5f %.5f %.4f\n",count,sJ/nb,sL/nb,(sJ-sL)/nb,sR/nb,sP/nb,ecc); sJ=sL=sR=sP=0; nb=0; } } } }
