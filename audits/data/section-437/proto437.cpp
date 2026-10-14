// Audit 437 (from 436): pullable (Adler-type) internal "action clock" theta' = p.v/hbar coupled to the orbit through U = -eps g(r) cos theta.
// 2D relative coordinate, hbar = mu = k = 1.  argv: eps J0 e0 orbits
#include <cmath>
#include <cstdio>
#include <cstdlib>
double gam=1e-4, kappa=0; const double r0=2.0;
struct S{double x,y,vx,vy,th;};
double eps;
S der(const S& s){ const double r=std::hypot(s.x,s.y), g=std::exp(-r/r0), gp=-g/r0;
  const double fr=-1/(r*r)*1.0 + eps*gp*std::cos(s.th);      // radial: Coulomb (-1/r^2) and -dU/dr = eps g'(r) cos theta
  // -dU/dr with U = -eps g cos th: F_r = eps g'(r) cos th
  return {s.vx,s.vy, fr*s.x/r-gam*s.vx, fr*s.y/r-gam*s.vy, s.vx*s.vx+s.vy*s.vy-kappa*eps*g*std::sin(s.th)}; }
int main(int argc,char** argv){
  eps=std::atof(argv[1]); gam=std::atof(argv[5]); kappa=std::atof(argv[6]); const double J0=std::atof(argv[2]), e0=std::atof(argv[3]); const int orbits=std::atoi(argv[4]);
  // start at periapsis of the Kepler orbit with a = J0^2, e0
  const double a=J0*J0, rp=a*(1-e0), vp=std::sqrt((1+e0)/(a*(1-e0)));
  S s{rp,0,0,vp,0};
  double t=0, tPrev=0, thPrev=0; int count=0; double Jsum=0, Rsum=0; int nb=0; double prevAngle=0, angleAcc=0;
  std::printf("# orbit t <J>_1000 L e <theta per orbit/2pi>_1000\n");
  while(count<orbits){
    const double r=std::hypot(s.x,s.y), E=0.5*(s.vx*s.vx+s.vy*s.vy)-1/r;
    if(!(E<0)) { std::printf("# unbound\n"); break; }
    const double J=1/std::sqrt(-2*E), T=2*M_PI*J*J*J, dt=T/2000;
    auto add=[&](const S& a,const S& b,double h){ return S{a.x+h*b.x,a.y+h*b.y,a.vx+h*b.vx,a.vy+h*b.vy,a.th+h*b.th}; };
    const S k1=der(s), k2=der(add(s,k1,dt/2)), k3=der(add(s,k2,dt/2)), k4=der(add(s,k3,dt));
    S n{s.x+dt/6*(k1.x+2*k2.x+2*k3.x+k4.x), s.y+dt/6*(k1.y+2*k2.y+2*k3.y+k4.y), s.vx+dt/6*(k1.vx+2*k2.vx+2*k3.vx+k4.vx),
        s.vy+dt/6*(k1.vy+2*k2.vy+2*k3.vy+k4.vy), s.th+dt/6*(k1.th+2*k2.th+2*k3.th+k4.th)};
    const double ang0=std::atan2(s.y,s.x), ang1=std::atan2(n.y,n.x); double d=ang1-ang0; if(d>M_PI) d-=2*M_PI; if(d<-M_PI) d+=2*M_PI;
    angleAcc+=d; s=n; t+=dt;
    if(angleAcc>=2*M_PI){ angleAcc-=2*M_PI; ++count;
      const double rr=std::hypot(s.x,s.y), EE=0.5*(s.vx*s.vx+s.vy*s.vy)-1/rr, JJ=1/std::sqrt(-2*EE), L=s.x*s.vy-s.y*s.vx;
      const double ecc=std::sqrt(std::max(0.0,1-L*L/(JJ*JJ)));
      const double ratio=(s.th-thPrev)/(2*M_PI);     // theta advance per orbit / 2 pi  (= J/hbar if unperturbed)
      Jsum+=JJ; Rsum+=ratio; ++nb;
      if(count%1000==0){ std::printf("%d %.4f %.6f %.6f %.5f %.6f\n",count,t,Jsum/nb,L,ecc,Rsum/nb); Jsum=Rsum=0; nb=0; }
      thPrev=s.th; tPrev=t; } }
}
