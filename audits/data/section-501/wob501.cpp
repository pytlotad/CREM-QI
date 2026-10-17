#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <vector>
// Audit 501: spin wobble on a Kepler orbit with the engine's instantaneous Thomas-BMT field -- forced response or resonance?
// argv: pair, n, L/hbar.  Particle 2 (electron for H) is followed; partner dipole fixed.
static Vec3 rot(const Vec3& v,const Vec3& w,double t){ const double a=w.norm()*t; if(a<1e-300) return v; const Vec3 k=w*(1/w.norm());
  return v*std::cos(a)+cross(k,v)*std::sin(a)+k*(dot(k,v)*(1-std::cos(a))); }
static double angle(const Vec3& a,const Vec3& b){ return std::atan2(cross(a,b).norm(),dot(a,b)); }
int main(int argc,char**argv){
  applyPairFromOption(argv[1]); const double nq=std::atof(argv[2]), Lh=std::atof(argv[3]);
  const double aB=pairBohrRadius(activePair), M=firstMass+secondMass, mu=firstMass*secondMass/M, a=nq*nq*aB, k=pairCoulombStrength;
  const double nmot=std::sqrt(k/(mu*a*a*a)), T=2*pi/nmot, ecc=std::sqrt(std::max(0.0,1-std::pow(Lh/nq,2))), ec=std::sqrt(1-ecc*ecc);
  const Vec3 L{0,0,1}, P{1,0,0}, Q{0,1,0};
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  const Vec3 sdir{std::sin(pi/3),0,std::cos(pi/3)};
  const Vec3 m1=sdir*(g1*0.5*hbar), m2=sdir*(g2*0.5*hbar);
  const int NE=8192; std::vector<Vec3> Om(NE); std::vector<double> tE(NE+1); tE[0]=0;
  for(int j=0;j<NE;++j){ const double E=2*pi*j/NE, w=1-ecc*std::cos(E);
    const Vec3 rr=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E)), vv=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w);
    State s{}; s.firstPosition=rr*(secondMass/M); s.secondPosition=rr*(-firstMass/M); s.firstVelocity=vv*(secondMass/M); s.secondVelocity=vv*(-firstMass/M);
    s.firstProperDipole=m1; s.secondProperDipole=m2; synchronizeCovariantDipoles(s);
    const StateHistory h{State{s}}; const LocalElectromagneticFields f=localRelativisticFields(s,h);
    Om[j]=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,secondGFactor)*(-secondCharge/secondMass);
    tE[j+1]=tE[j]+(2*pi/NE)*(1-ecc*std::cos(E+pi/NE))/nmot; }               // time at node j+1 (midpoint rule)
  // orbit average and harmonics in time
  Vec3 avg{}; double Tsum=0; for(int j=0;j<NE;++j){ const double dt=tE[j+1]-tE[j]; avg=avg+Om[j]*dt; Tsum+=dt; } avg=avg*(1/Tsum);
  std::printf("%s n %.0f L %.1f hbar e %.4f: T %.4e s, omega_orb %.4e, |<Omega>| %.4e rad/s, ratio %.3e; r_p %.4e m, (3/4) r_e/r_p %.3e\n",
    argv[1],nq,Lh,ecc,T,nmot,avg.norm(),avg.norm()/nmot,a*(1-ecc),0.75*2.8179403262e-15/(a*(1-ecc)));
  double c1perp=0,c2perp=0; std::printf("  harmonics |c_k|/|<Omega>|:"); for(int kk=1;kk<=16;kk*=2){ std::complex<double> cx(0),cy(0),cz(0);
    for(int j=0;j<NE;++j){ const double t=0.5*(tE[j]+tE[j+1]), dt=tE[j+1]-tE[j]; const std::complex<double> ph=std::exp(std::complex<double>(0,-kk*nmot*t))*dt;
      cx+=Om[j].x*ph; cy+=Om[j].y*ph; cz+=Om[j].z*ph; }
    const Vec3 nh=avg*(1/avg.norm()); const std::complex<double> cpar=cx*nh.x+cy*nh.y+cz*nh.z;
    const double tot=std::sqrt(std::norm(cx)+std::norm(cy)+std::norm(cz)), perp=std::sqrt(std::max(0.0,tot*tot-std::norm(cpar)));
    std::printf(" k%d %.3e (perp %.3e)",kk,tot/Tsum/avg.norm(),perp/Tsum/avg.norm());
    if(kk==1) c1perp=perp/Tsum; if(kk==2) c2perp=perp/Tsum; } std::printf("\n");
  std::printf("  transverse drive: |c_1,perp| %.4e, |c_2,perp| %.4e rad/s; resonant tilt after 50 orbits ~ |c_perp|/2 x 50 T = %.3e / %.3e rad\n",
    c1perp,c2perp,0.5*c1perp*50*T,0.5*c2perp*50*T);
  // linear response amplitude perpendicular to S0
  const Vec3 S0=sdir; Vec3 W{}; double Alin=0;
  for(int j=0;j<NE;++j){ W=W+(Om[j]-avg)*(tE[j+1]-tE[j]); const Vec3 Wp=W-S0*dot(W,S0); Alin=std::max(Alin,Wp.norm()); }
  std::printf("  linear wobble amplitude A_lin = %.4e rad (/(3/4 r_e/r_p) = %.3f)\n",Alin,Alin/(0.75*2.8179403262e-15/(a*(1-ecc))));
  // nonlinear integration over Norb orbits with scale lambda
  auto run=[&](double lam,int Norb,double& first,double& last,double mu=-1){ Vec3 S=S0; first=0; last=0; double t=0;
    const bool shift=mu>0; const Vec3 ref=shift?avg*mu:avg*lam; auto Ov=[&](int j){ return shift?Om[j%NE]-avg+avg*mu:Om[j%NE]*lam; };
    for(int o=0;o<Norb;++o){ double mx=0;
      for(int j=0;j<NE;j+=2){ const double h=tE[j+2]-tE[j]; const Vec3 w0=Ov(j), w1=Ov(j+1), w2=Ov(j+2);
        auto f=[&](const Vec3& w,const Vec3& s){ return cross(w,s); };
        const Vec3 k1=f(w0,S), k2=f(w1,S+k1*(0.5*h)), k3=f(w1,S+k2*(0.5*h)), k4=f(w2,S+k3*h);
        S=S+(k1+k2*2+k3*2+k4)*(h/6); S=S*(1/S.norm()); t+=h;
        mx=std::max(mx,angle(S,rot(S0,ref,t))); }
      if(o==0) first=mx; if(o==Norb-1) last=mx; } };
  double f1,l1; run(1.0,50,f1,l1);
  std::printf("  lambda 1: max deviation orbit 1 %.4e, orbit 50 %.4e rad; /A_lin %.4f; growth 50/1 %.4f\n",f1,l1,f1/Alin,l1/f1);
  const double lk1=nmot/avg.norm();
  for(double lam: {1e2,1e3,1e4,3e4,0.98*lk1,lk1,1.02*lk1,2*lk1}){ double f,l; run(lam,50,f,l);
    std::printf("  lambda %.4e (%.3f lambda_1): orbit 1 %.4e, orbit 50 %.4e rad, growth %.3f\n",lam,lam/lk1,f,l,l/f); }
  // (e) oscillating part unchanged, mean shifted to mu <Omega>
  for(double mu: {1.0,0.5*lk1,0.98*lk1,lk1,1.02*lk1,2*lk1}){ double f,l; run(1.0,50,f,l,mu);
    std::printf("  mean shifted mu %.4e (%.3f mu_1): orbit 1 %.4e, orbit 50 %.4e rad, growth %.3f\n",mu,mu/lk1,f,l,l/f); }
}
