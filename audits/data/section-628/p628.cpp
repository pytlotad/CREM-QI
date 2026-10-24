// Audit 628: rule Q as the interference of the phase S/hbar over successive orbits, on the engine's own spin-free Hamiltonian.
// MODE=levels: E_EBK(n, l) from J_r(E, L = (l + 1/2) hbar) = (n_r + 1/2) hbar, J_r = (1/pi) int p_r dr with the engine's energy function
//   (conservativeParticleEnergy, moments zero) and canonical momenta p_i = gamma_i m_i v_i + dE_D/dv_i (darwinInteractionEnergy, homogeneous of
//   degree 2 in the velocities, so L_D = E_D); the pair in its centre-of-mass frame, p_1 = -p_2.
// MODE=sum: |sum_{N<Nmax} exp(i N (2 pi J_r(E)/hbar - pi))| / Nmax on an energy grid around the levels n = 2 (l = 0, 1).
// MODE=orbit: J_r measured along an engine trajectory (stepReady) over successive radial periods, against the quadrature at the run's E, L.
// env PAIR, MODE, NMAX (levels), NLEV, LHBAR, PERIODS (orbit).
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>
static double envd(const char* k,double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }

static State makeState(double r,const Vec3& v1,const Vec3& v2){
  State s; const double M=firstMass+secondMass;
  s.firstPosition={r*secondMass/M,0,0}; s.secondPosition={-r*firstMass/M,0,0};
  s.firstVelocity=v1; s.secondVelocity=v2; s.firstDipole={0,0,0}; s.secondDipole={0,0,0};
  return s;
}
// dE_D/dv_i by central differences (E_D is bilinear in v1, v2, so the difference is exact up to rounding)
static Vec3 darwinGrad(const State& s,bool first){
  Vec3 g; const double h=1.0e3;   // m/s; E_D is linear in each velocity, so any h is exact up to rounding
  for(int k=0;k<3;++k){ State a=s, b=s; Vec3 e{0,0,0}; (k==0?e.x:(k==1?e.y:e.z))=h;
    if(first){ a.firstVelocity=s.firstVelocity+e; b.firstVelocity=s.firstVelocity-e; }
    else { a.secondVelocity=s.secondVelocity+e; b.secondVelocity=s.secondVelocity-e; }
    const double d=(darwinInteractionEnergy(a)-darwinInteractionEnergy(b))/(2*h);
    (k==0?g.x:(k==1?g.y:g.z))=d; }
  return g;
}
static Vec3 canonical(const State& s,bool first){
  const Vec3 v=first?s.firstVelocity:s.secondVelocity; const double m=first?firstMass:secondMass;
  const double gam=1.0/std::sqrt(1.0-v.squaredNorm()/(c*c));
  return v*(gam*m)+darwinGrad(s,first);
}
// energy of the pair at separation r with relative canonical momentum p (p_1 = p, p_2 = -p), r along x
static double hamiltonian(double r,const Vec3& p){
  Vec3 v1=velocityFromMomentum(p,firstMass), v2=velocityFromMomentum(-p,secondMass);
  for(int it=0;it<60;++it){
    const State s=makeState(r,v1,v2);
    const Vec3 n1=velocityFromMomentum(p-darwinGrad(s,true),firstMass);
    const Vec3 n2=velocityFromMomentum(-p-darwinGrad(s,false),secondMass);
    const double d=(n1-v1).norm()+(n2-v2).norm(); v1=n1; v2=n2;
    if(d<1e-15*(v1.norm()+v2.norm())) break; }
  return conservativeParticleEnergy(makeState(r,v1,v2));
}
static double prOf(double E,double r,double L){
  const double pt=L/r;
  if(hamiltonian(r,{0,pt,0})>=E) return 0.0;
  double lo=0, hi=std::sqrt(2*pairReducedMass*std::max(1e-40,E-hamiltonian(r,{0,pt,0})))*2+1e-40;
  while(hamiltonian(r,{hi,pt,0})<E) hi*=2;
  for(int k=0;k<200;++k){ const double mid=0.5*(lo+hi); if(hamiltonian(r,{mid,pt,0})<E) lo=mid; else hi=mid; if(hi-lo<=1e-16*hi) break; }
  return 0.5*(lo+hi);
}
static std::vector<double> gx, gw;
static void gauss(int n){ gx.resize(n); gw.resize(n);
  for(int i=0;i<n;++i){ double x=std::cos(pi*(i+0.75)/(n+0.5)), dp=0;
    for(int it=0;it<100;++it){ double p0=1,p1=x; for(int k=2;k<=n;++k){ const double p2=((2*k-1)*x*p1-(k-1)*p0)/k; p0=p1; p1=p2; }
      dp=n*(x*p1-p0)/(x*x-1); const double dx=p1/dp; x-=dx; if(std::abs(dx)<1e-16) break; }
    gx[i]=x; gw[i]=2/((1-x*x)*dp*dp); } }
static double turning(double E,double L,double a,double b){   // f(r) = H(r, 0, L/r) - E changes sign in [a, b]
  double fa=hamiltonian(a,{0,L/a,0})-E;
  for(int k=0;k<200;++k){ const double m=0.5*(a+b), fm=hamiltonian(m,{0,L/m,0})-E;
    if((fm<0)==(fa<0)){ a=m; fa=fm; } else b=m; if(b-a<=1e-16*b) break; }
  return 0.5*(a+b);
}
static double radialAction(double E,double L){
  const double k=pairCoulombStrength, mu=pairReducedMass;
  const double a=k/(2*std::abs(E)), rc=L*L/(mu*k);             // Kepler scales for the brackets
  const double rmin=turning(E,L,rc*1e-3,rc), rmax=turning(E,L,rc,4*a+rc);
  const double cen=0.5*(rmin+rmax), half=0.5*(rmax-rmin);
  double s=0;
  for(size_t i=0;i<gx.size();++i){ const double t=0.5*pi*(gx[i]+1); const double r=cen-half*std::cos(t);
    s+=gw[i]*0.5*pi*prOf(E,r,L)*half*std::sin(t); }
  return s/pi;
}
static double levelEBK(int n,int l){
  const double k=pairCoulombStrength, mu=pairReducedMass, E0=-mu*k*k/(2*hbar*hbar*n*n);
  const double L=(l+0.5)*hbar, target=(n-l-1+0.5)*hbar;
  double lo=E0*(1+2e-3), hi=E0*(1-2e-3);                         // E0 < 0: lo is deeper
  for(int it=0;it<200;++it){ const double m=0.5*(lo+hi); if(radialAction(m,L)<target) lo=m; else hi=m; if(hi-lo<=4e-16*std::abs(m)) break; }
  return 0.5*(lo+hi);
}
int main(){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  // uzupelnienie 1: alpha_eff = s alpha (Coulomb and Darwin together), and moments of length zero for the orbit check
  if(const double s=envd("ALPHASCALE",1.0); s!=1.0){ pairCoulombStrength*=s; pairChargeProduct*=s; }
  if(envd("DARWIN",1)==0) pairChargeProduct=0.0;   // uzupelnienie 2: Darwin term off (Coulomb uses pairCoulombStrength)
  if(envd("ZEROMOMENTS",0)!=0){ firstMagneticMoment=0.0; secondMagneticMoment=0.0; }
  gauss((int)envd("NGAUSS",160));
  const std::string mode=std::getenv("MODE")?std::getenv("MODE"):"levels";
  const double k=pairCoulombStrength, mu=pairReducedMass;
  std::printf("CONST m1 %.15e m2 %.15e mu %.15e k %.15e hbar %.15e c %.15e q1q2 %.15e\n",firstMass,secondMass,mu,k,hbar,c,pairChargeProduct);
  // the orbit check also reports the spread of E and L along the trajectory
  if(mode=="check"){   // the probe's Hamiltonian against the analytic two-body form (Darwin off)
    for(double fr: {0.3,1.0,3.0}) for(double fp: {0.2,1.0}){
      const double a0=hbar*hbar/(mu*k), r=fr*a0, P=fp*hbar/a0; const Vec3 p{0.6*P,0.8*P,0};
      const double h=hamiltonian(r,p);
      const double an=std::sqrt(P*P*c*c+firstMass*firstMass*c*c*c*c)+std::sqrt(P*P*c*c+secondMass*secondMass*c*c*c*c)
                      -(firstMass+secondMass)*c*c-k/r;
      std::printf("CHECK r/a0 %.1f p a0/hbar %.1f H %.17e analytic %.17e rel %.3e\n",fr,fp,h,an,(h-an)/std::abs(an)); }
    return 0;
  }
  if(mode=="levels"){
    const int nmax=(int)envd("NMAX",5);
    for(int n=1;n<=nmax;++n) for(int l=0;l<n;++l){
      const double E0=-mu*k*k/(2*hbar*hbar*n*n), E=levelEBK(n,l);
      std::printf("LEVEL n %d l %d E0 %.17e E_EBK %.17e dE %.17e\n",n,l,E0,E,E-E0); std::fflush(stdout); }
  } else if(mode=="sum"){
    const double E0=-mu*k*k/(2*hbar*hbar*4);
    const double El0=levelEBK(2,0), El1=levelEBK(2,1);
    std::printf("EBK n=2: l=0 %.17e l=1 %.17e (E0 %.17e)\n",El0,El1,E0);
    for(double Nmax: {1.0e3,1.0e5,1.0e6}) for(int l=0;l<=1;++l){
      const double L=(l+0.5)*hbar; double best=-1, bestE=0;
      const double span=std::abs(El1-El0)*4, Ec=0.5*(El0+El1);
      for(int i=0;i<=400;++i){ const double E=Ec-span/2+span*i/400;
        const double ph=2*pi*radialAction(E,L)/hbar-pi;           // phase per radial period
        const double x=std::remainder(ph,2*pi);
        const double amp=std::abs(x)<1e-300?1.0:std::abs(std::sin(0.5*Nmax*x)/(Nmax*std::sin(0.5*x)));
        if(amp>best){ best=amp; bestE=E; }
        if(Nmax==1.0e6) std::printf("SUM Nmax %.0e l %d E %.17e amp %.6e\n",Nmax,l,E,amp); }
      std::printf("PEAK Nmax %.0e l %d E_peak %.17e (E_EBK %.17e, grid step %.3e)\n",Nmax,l,bestE,l?El1:El0,span/400); }
  } else {
    // engine trajectory: conservative Coulomb-Darwin between photons, no photons, no dipole forces (env set by the caller)
    gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
    const int n=(int)envd("NLEV",3), l=(int)envd("LHBAR",1);
    gInitialSeparationOverride=n*n*pairBohrRadius(activePair); gInitialAngularMomentumFraction=(l+0.5)/n;
    const double a=n*n*pairBohrRadius(activePair), T=2*pi*std::sqrt(a*a*a*mu/k);
    SimulationOptions o; o.frameCount=2; o.observationTime=envd("PERIODS",4)*T; o.radiatedEnergyBookkeeping=false;
    struct Acc { bool started=false; double rPrev=0, prPrev=0, phi=0; int sign=0, turns=0; double sE=0, sL=0; long n=0;
                 std::vector<double> perPeriod; } acc;
    o.stepReady=[&](const State& s){
      const Vec3 d=s.firstPosition-s.secondPosition; const double r=d.norm(); const Vec3 nh=d*(1.0/r);
      const double M=firstMass+secondMass;
      const Vec3 p=(canonical(s,true)*secondMass-canonical(s,false)*firstMass)*(1.0/M);   // relative canonical momentum
      const double pr=dot(p,nh);
      State z=s; z.firstDipole={0,0,0}; z.secondDipole={0,0,0}; z.dipoleConstraintEnergy=0;
      acc.sE+=conservativeParticleEnergy(z); acc.sL+=cross(d,p).norm(); ++acc.n;
      if(acc.started){ acc.phi+=0.5*(pr+acc.prPrev)*(r-acc.rPrev);
        const int sg=(r>acc.rPrev)?1:-1;
        if(acc.sign!=0&&sg!=acc.sign&&sg==1){   // periapsis passed: one radial period completed since the previous one
          if(acc.turns>0) acc.perPeriod.push_back(acc.phi/(2*pi)); acc.phi=0; ++acc.turns; }
        acc.sign=sg; }
      acc.started=true; acc.rPrev=r; acc.prPrev=pr; };
    simulate(splitMix64(628ULL),1,o);
    const double E=acc.sE/acc.n, L=acc.sL/acc.n, Jq=radialAction(E,L);
    std::printf("ORBIT n %d l %d steps %ld periods %zu <E> %.15e <L>/hbar %.12f J_r quadrature/hbar %.12f\n",n,l,acc.n,acc.perPeriod.size(),E,L/hbar,Jq/hbar);
    for(double j: acc.perPeriod) std::printf("ORBIT period J_r/hbar %.12f  diff %.3e\n",j/hbar,(j-Jq)/hbar);
  }
}
