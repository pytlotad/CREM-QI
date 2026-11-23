// Audit 833 (from ../section-570/p570j.cpp): the two Thomas expressions of the engine compared statically on an exact Kepler orbit.
//   C  Thomas part of covariantSpinMomentumForce (retarded field, Tulczyjew momentum, force on the spin carrier only):
//      covariantSpinMomentumForce with minus without gSpinNoThomas (the hidden part cancels exactly);
//   L  thomasBackReactionForces (Lagrangian L_T, kappa 1/2, gradient term and reaction on the partner);
//   N  analytic instantaneous Tulczyjew: -d/dt [q_p (E_Coulomb x S_p)/(m_p c^2)] along the exact orbit, no boost.
// Each for the first particle's spin, the second's (for L: partner dipole zeroed) and L for both.  First-order apse rate: time
// average of the rotation of the Laplace-Runge-Lenz vector about L, Omega = <Lh.(Ah x dA/dt)>/|A|, dA/dt = a x h + v x (r x a),
// a = F1/m1 - F2/m2.  Analytic L: H_T = k q_o q_p (S.L)/(2 m_p^2 c^2 r^3), Omega = -2 <H_T>/|L|.
// argv: pair  cos1  cos2  LFRAC  tag   (cos = spin direction . L; second spin = sign(cos2/cos1) x first, same azimuth)
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  applyPairFromOption(argv[1]);
  const double cos1=std::atof(argv[2]), cos2=std::atof(argv[3]), lfrac=std::atof(argv[4]); const char* tag=argv[5];
  const double nlev=2.0, Lh=nlev*lfrac;
  const double aB=pairBohrRadius(activePair), M=firstMass+secondMass, mu=firstMass*secondMass/M;
  const double a=nlev*nlev*aB, k=pairCoulombStrength;
  const double nmot=std::sqrt(k/(mu*a*a*a)), ecc=std::sqrt(std::max(0.0,1-lfrac*lfrac)), ec=std::sqrt(1-ecc*ecc);
  const double T=2*pi/nmot;
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  const Vec3 P{1,0,0}, Q{0,1,0}, Z{0,0,1};
  const auto relative=[&](double t,Vec3& r,Vec3& v,Vec3& acc){
    const double Mm=nmot*t; double E=Mm;
    for(int it=0;it<60;++it){ const double d=(E-ecc*std::sin(E)-Mm)/(1-ecc*std::cos(E)); E-=d; if(std::abs(d)<1e-15) break; }
    const double w=1-ecc*std::cos(E);
    r=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E));
    v=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w);
    acc=r*(-k/(mu*std::pow(r.norm(),3)));
  };
  const double sin1=std::sqrt(std::max(0.0,1-cos1*cos1)), sgn=(cos1*cos2<0)?-1.0:1.0;
  const int NAZ=8, N=512, H=64;
  // analytic L (sum over spins), rad/s
  const double rm3=1.0/(a*a*a*std::pow(1-ecc*ecc,1.5));
  const double an1=k*(0.5*hbar)*cos1*rm3/(firstMass*firstMass*c*c), an2=k*(0.5*hbar)*cos2*rm3/(secondMass*secondMass*c*c);
  std::printf("SYS %s pair %s LFRAC %.17g e %.12f a %.12e T %.12e an1 %.12e an2 %.12e\n",tag,argv[1],lfrac,ecc,a,T,an1,an2);
  for(int iz=0;iz<NAZ;++iz){
    const double phi=2*pi*iz/NAZ;
    const Vec3 s1{sin1*std::cos(phi),sin1*std::sin(phi),cos1}, s2=s1*sgn;
    const auto build=[&](double t,int which){   // which: 0 both spins, 1 only first, 2 only second
      Vec3 r,v,acc; relative(t,r,v,acc);
      State s{}; s.time=t;
      s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M);
      s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
      s.firstAcceleration=acc*(secondMass/M); s.secondAcceleration=acc*(-firstMass/M);
      s.firstProperDipole=which==2?Vec3{}:s1*(g1*0.5*hbar); s.secondProperDipole=which==1?Vec3{}:s2*(g2*0.5*hbar);
      synchronizeCovariantDipoles(s); return s; };
    const auto hist=[&](double t,int which){ StateHistory h; for(int i=H;i>=1;--i) h.push_back(build(t-i*(T/4096.0),which)); return h; };
    // analytic Tulczyjew momentum of spin i at time t (instantaneous Coulomb field of the partner)
    const auto ptul=[&](double t,bool first){
      Vec3 r,v,acc; relative(t,r,v,acc);
      const Vec3 d=first?r:r*-1.0;   // from partner to spin carrier
      const double qo=first?secondCharge:firstCharge, qp=first?firstCharge:secondCharge, mp=first?firstMass:secondMass;
      const Vec3 E=d*(coulomb*qo/std::pow(d.norm(),3));
      const Vec3 S=(first?s1:s2)*(0.5*hbar);
      return cross(E,S)*(qp/(mp*c*c)); };
    double sw=0, om[8]={0,0,0,0,0,0,0,0};
    for(int j=0;j<N;++j){
      const double Ej=2*pi*(j+0.5)/N, wj=1-ecc*std::cos(Ej);
      const double t=(Ej-ecc*std::sin(Ej))/nmot;
      Vec3 r,v,acc; relative(t,r,v,acc);
      const Vec3 hv=cross(r,v);
      const double Aabs=k*ecc/mu;
      const auto rate=[&](const Vec3& F1,const Vec3& F2){
        const Vec3 ap=F1*(1.0/firstMass)-F2*(1.0/secondMass);
        const Vec3 dA=cross(ap,hv)+cross(v,cross(r,ap));
        return dot(Z,cross(P,dA))/Aabs; };
      const State s=build(t,0); const StateHistory h=hist(t,0);
      gSpinNoThomas=false;
      const Vec3 c1on=covariantSpinMomentumForce(s,h,true), c2on=covariantSpinMomentumForce(s,h,false);
      gSpinNoThomas=true;
      const Vec3 c1off=covariantSpinMomentumForce(s,h,true), c2off=covariantSpinMomentumForce(s,h,false);
      gSpinNoThomas=false;
      const MutualForces lb=thomasBackReactionForces(s,thomasBmtDipoleDerivatives(s,h));
      const State s1o=build(t,1); const StateHistory h1o=hist(t,1);
      const MutualForces l1=thomasBackReactionForces(s1o,thomasBmtDipoleDerivatives(s1o,h1o));
      const State s2o=build(t,2); const StateHistory h2o=hist(t,2);
      const MutualForces l2=thomasBackReactionForces(s2o,thomasBmtDipoleDerivatives(s2o,h2o));
      const double hh=T/1.0e5;
      const Vec3 n1=(ptul(t+hh,true)-ptul(t-hh,true))*(-1.0/(2*hh)), n2=(ptul(t+hh,false)-ptul(t-hh,false))*(-1.0/(2*hh));
      const double x[7]={rate(c1on-c1off,Vec3{}),rate(Vec3{},c2on-c2off),rate(l1.first,l1.second),rate(l2.first,l2.second),
                         rate(lb.first,lb.second),rate(n1,Vec3{}),rate(Vec3{},n2)};
      for(int q=0;q<7;++q) om[q]+=x[q]*wj;
      sw+=wj;
    }
    std::printf("OM %s phi %d C1 %.12e C2 %.12e L1 %.12e L2 %.12e Lb %.12e N1 %.12e N2 %.12e\n",tag,iz,
                om[0]/sw,om[1]/sw,om[2]/sw,om[3]/sw,om[4]/sw,om[5]/sw,om[6]/sw);
    std::fflush(stdout);
  }
}
