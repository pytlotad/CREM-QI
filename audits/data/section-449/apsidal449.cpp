// Audit 449: apsidal precession vs l-level spacing.  (1) identity on the engine's E_rel; (2) measured periapsis advance per orbit.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
struct St{double x,y,px,py;};
int main(){
  // (1) engine identity (Ps), contact off via env
  const double hP=2*pi*hbar, GHz=1e-9/hP, aB=pairBohrRadius(activePair), alpha=pairCoulombStrength/(hbar*c), a4=std::pow(alpha,4)*firstMass*c*c;
  std::printf("# (1) engine E_rel: E(l+1) - E(l) vs hbar omega_w at L = sqrt((l+1/2)(l+3/2)) hbar\n");
  for(int n=2;n<=4;++n) for(int l=0;l<=n-2;++l){
    const double a=aB*n*n, E0=levelFineStructureEnergy(a,(l+0.5)*hbar), E1=levelFineStructureEnergy(a,(l+1.5)*hbar);
    const double Lm=std::sqrt((l+0.5)*(l+1.5)), dL=1e-6, Ep=levelFineStructureEnergy(a,(Lm+dL)*hbar), Em=levelFineStructureEnergy(a,(Lm-dL)*hbar);
    const double hw=(Ep-Em)/(2*dL), formula=a4/(4.0*n*n*n*Lm*Lm);
    std::printf("n=%d l=%d->%d: dE %.9f GHz, hbar omega_w (numeric) %.9f, (formula) %.9f, rel diff %.2e\n",n,l,l+1,(E1-E0)*GHz,hw*GHz,formula*GHz,(E1-E0)/hw-1);
  }
  // (2) dynamics: Darwin Hamiltonian, units hbar = m_e = k = 1, c = 1/alpha, mu = 1/2
  std::printf("# (2) measured periapsis advance per orbit vs pi alpha^2/L^2\n");
  for(double al: {7.2973525693e-3, 0.05}){ const double cc=1/al;
    auto H_dot=[&](const St& s){ const double r=std::sqrt(s.x*s.x+s.y*s.y), p2=s.px*s.px+s.py*s.py, pr=(s.x*s.px+s.y*s.py)/r;
      // H = p^2 - 1/r - p^4/(4c^2) - (p^2 + pr^2)/(2 c^2 r)   (p^2/(2 mu) = p^2 with mu = 1/2)
      const double dHdp2=1-2*p2/(4*cc*cc)-1/(2*cc*cc*r);
      // derivatives of pr^2 = (x.p)^2/r^2
      const double xp=s.x*s.px+s.y*s.py, r2=r*r;
      const double dpr2_dpx=2*xp*s.x/r2, dpr2_dpy=2*xp*s.y/r2;
      const double dpr2_dx=2*xp*s.px/r2-2*xp*xp*s.x/(r2*r2), dpr2_dy=2*xp*s.py/r2-2*xp*xp*s.y/(r2*r2);
      const double c2r=2*cc*cc*r;
      St d; d.x=dHdp2*2*s.px-dpr2_dpx/c2r; d.y=dHdp2*2*s.py-dpr2_dpy/c2r;
      const double dHdr=1/(r*r)+(p2+pr*pr)/(c2r*r);               // explicit r-derivative of -1/r - (p2+pr2)/(2c^2 r)
      d.px=-(dHdr*s.x/r - dpr2_dx/c2r); d.py=-(dHdr*s.y/r - dpr2_dy/c2r); return d; };
    for(auto [n,L]: {std::pair{2.0,0.8660254},{2.0,1.5},{3.0,2.5}}){
      // start at periapsis of the Kepler orbit with a = 2 n^2, L (relative momentum p = mu v)
      const double a=2*n*n, e=std::sqrt(1-L*L/(n*n)), rp=a*(1-e), pp=L/rp;
      St s{rp,0,0,pp}; const double T=4*M_PI*n*n*n, dt=T/200000; double prevAng=0, firstAng=0; int peri=0; double prd=-1;
      for(long i=0;i<200000L*32 && peri<31;++i){
        auto k1=H_dot(s); auto add=[&](const St& u,const St& v,double h){ return St{u.x+h*v.x,u.y+h*v.y,u.px+h*v.px,u.py+h*v.py}; };
        auto k2=H_dot(add(s,k1,dt/2)), k3=H_dot(add(s,k2,dt/2)), k4=H_dot(add(s,k3,dt));
        St nn{s.x+dt/6*(k1.x+2*k2.x+2*k3.x+k4.x), s.y+dt/6*(k1.y+2*k2.y+2*k3.y+k4.y), s.px+dt/6*(k1.px+2*k2.px+2*k3.px+k4.px), s.py+dt/6*(k1.py+2*k2.py+2*k3.py+k4.py)};
        const double rd0=s.x*s.px+s.y*s.py, rd1=nn.x*nn.px+nn.y*nn.py;
        if(rd0<0&&rd1>=0){ const double fr=-rd0/(rd1-rd0); const double ang=std::atan2(s.y+fr*(nn.y-s.y),s.x+fr*(nn.x-s.x)); if(peri==0) firstAng=ang; prevAng=ang; ++peri; }
        s=nn; }
      double adv=prevAng-firstAng; while(adv<0) adv+=2*M_PI; adv/= (peri-1);
      std::printf("alpha %.6f n %.0f L %.4f: measured %.6e rad/orbit, pi alpha^2/L^2 %.6e, ratio %.5f\n",al,n,L,adv,M_PI*al*al/(L*L),adv/(M_PI*al*al/(L*L)));
    } }
}
