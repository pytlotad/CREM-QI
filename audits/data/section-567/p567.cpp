#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 567 (from p565): the torque of the WHOLE dipole sector of the production
// function retardedExternalForces (dipole force on minus off), against -dS/dt.
// Audit 565: the force on the moments against the Frenkel / classical
// Foldy-Wouthuysen Stern-Gerlach force (Wen, Keitel, Bauke, PRA 95, 042102
// (2017), eqs. (21)-(30)) on an exact Kepler orbit with a synthetic history.
//   A engine (retarded sum, production): covariantDipoleGradientForce
//     + hiddenMomentumRateForce
//   B A + thomasBackReactionForces (what CREM_RETARDED_THOMAS adds)
//   C FW: -grad H_i, H_i = -m_i.B_j + ((g_i-1)/g_i) m_i.(v_i x E_j)/c^2
//   D C + hiddenMomentumRateForce
// plus, common to all, the Lorentz force of each moment's field on the
// partner's charge.  Measures: torque against -dS/dt (consistency of force
// and precession), the E1 power interference P_cross/P_E1 (as audit 507),
// the Frenkel effective-mass shift, and |C - A| against the Thomas gradient
// minus the hidden-momentum rate.
// argv: pair  n  L/hbar
int main(int argc,char**argv){
  applyPairFromOption(argv[1]);
  const double nlev=std::atof(argv[2]), Lh=std::atof(argv[3]);
  const double aB=pairBohrRadius(activePair), M=firstMass+secondMass, mu=firstMass*secondMass/M;
  const double a=nlev*nlev*aB, k=pairCoulombStrength, Lm=Lh*hbar;
  const double nmot=std::sqrt(k/(mu*a*a*a)), ecc=std::sqrt(std::max(0.0,1-std::pow(Lm/(nlev*hbar),2))), ec=std::sqrt(1-ecc*ecc);
  const double T=2*pi/nmot;
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  const double gf1=std::abs(firstGFactor), gf2=std::abs(secondGFactor);
  const Vec3 P{1,0,0}, Q{0,1,0}, Z{0,0,1};
  struct Cfg{const char* n; Vec3 s1,s2;};
  const Vec3 S60{std::sin(pi/3),0,std::cos(pi/3)};
  const Cfg cf[]={{"both along L",Z,Z},{"both against L",Z*-1,Z*-1},{"both perpendicular",P,P},
                  {"antiparallel, 1 along L",Z,Z*-1},{"antiparallel, perpendicular",P,P*-1},
                  {"both at 60 deg",S60,S60}};
  const double kL=1.0/(6*pi*8.8541878128e-12*c*c*c);
  // relative orbit at time t (exact Kepler), periapsis along P, t = 0 at periapsis
  const auto relative=[&](double t,Vec3& r,Vec3& v,Vec3& acc){
    const double Mm=nmot*t; double E=Mm;
    for(int it=0;it<60;++it){ const double d=(E-ecc*std::sin(E)-Mm)/(1-ecc*std::cos(E)); E-=d; if(std::abs(d)<1e-15) break; }
    const double w=1-ecc*std::cos(E);
    r=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E));
    v=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w);
    acc=r*(-k/(mu*std::pow(r.norm(),3)));
  };
  std::printf("%s: n %.0f, L %.2f hbar, e %.6f, a %.4e m, T %.4e s, g %.6f / %.6f\n",argv[1],nlev,Lh,ecc,a,T,gf1,gf2);
  for(const Cfg& q: cf){
    const auto build=[&](double t){
      Vec3 r,v,acc; relative(t,r,v,acc);
      State s{}; s.time=t;
      s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M);
      s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
      s.firstAcceleration=acc*(secondMass/M); s.secondAcceleration=acc*(-firstMass/M);
      s.firstProperDipole=q.s1*(g1*0.5*hbar); s.secondProperDipole=q.s2*(g2*0.5*hbar);
      synchronizeCovariantDipoles(s); return s; };
    const int N=1024, H=64;
    Vec3 tau[4], dS; double sw=0, sC=0, sX[4]={0,0,0,0}, maxMass=0, worstCheck=0, scaleCheck=0;
    Vec3 ds1, tauSector;
    for(int j=0;j<N;++j){
      // nodes uniform in the eccentric anomaly, weight dt ~ (1 - e cos E) (audit 565, uzupelnienie 1)
      const double Ej=2*pi*(j+0.5)/N, wj=1-ecc*std::cos(Ej);
      const double t=(Ej-ecc*std::sin(Ej))/nmot;
      const State s=build(t);
      StateHistory h; const double dt=T/4096.0;
      for(int i=H;i>=1;--i) h.push_back(build(t-i*dt));
      // common: charge-charge acceleration for P_E1, Lorentz of moment fields on charges
      const ElectromagneticField lw1=lienardWiechertField(s.firstPosition,t,h,s,false,secondCharge);
      const ElectromagneticField lw2=lienardWiechertField(s.secondPosition,t,h,s,true,firstCharge);
      const Vec3 fC1=lorentzForce(firstCharge,s.firstVelocity,lw1), fC2=lorentzForce(secondCharge,s.secondVelocity,lw2);
      const ElectromagneticField md1=retardedMagneticDipoleField(s.firstPosition,t,h,s,false);
      const ElectromagneticField md2=retardedMagneticDipoleField(s.secondPosition,t,h,s,true);
      const Vec3 fL1=lorentzForce(firstCharge,s.firstVelocity,md1), fL2=lorentzForce(secondCharge,s.secondVelocity,md2);
      // A, B
      const Vec3 grad1=covariantDipoleGradientForce(s,h,true), grad2=covariantDipoleGradientForce(s,h,false);
      const Vec3 hid1=hiddenMomentumRateForce(s,h,true), hid2=hiddenMomentumRateForce(s,h,false);
      const MutualForces th=thomasBackReactionForces(s,thomasBmtDipoleDerivatives(s,h));
      // C: FW gradient
      const auto Hof=[&](bool first,const Vec3& x){
        const ElectromagneticField f=fieldFromOtherParticleAt(x,t,s,h,first,1.0e-5,matchedMomentSoftening());
        const Vec3 m=first?s.firstDipole:s.secondDipole, v=first?s.firstVelocity:s.secondVelocity;
        const double gg=first?gf1:gf2;
        return -dot(m,f.magnetic)+((gg-1)/gg)*dot(m,cross(v,f.electric))/(c*c); };
      const auto fw=[&](bool first){
        const Vec3 x0=first?s.firstPosition:s.secondPosition; const double st=1e-4*(s.firstPosition-s.secondPosition).norm();
        Vec3 g; const Vec3 e[3]={{st,0,0},{0,st,0},{0,0,st}};
        g.x=-(Hof(first,x0+e[0])-Hof(first,x0-e[0]))/(2*st);
        g.y=-(Hof(first,x0+e[1])-Hof(first,x0-e[1]))/(2*st);
        g.z=-(Hof(first,x0+e[2])-Hof(first,x0-e[2]))/(2*st); return g; };
      const Vec3 fw1=fw(true), fw2=fw(false);
      // check: C - A(gradient) = grad[(1/g) m.(v x E)]/c^2  (Thomas part), computed independently
      const auto thomasGrad=[&](bool first){
        const Vec3 x0=first?s.firstPosition:s.secondPosition; const double st=1e-4*(s.firstPosition-s.secondPosition).norm();
        const Vec3 m=first?s.firstDipole:s.secondDipole, v=first?s.firstVelocity:s.secondVelocity; const double gg=first?gf1:gf2;
        const auto U=[&](const Vec3& x){ const ElectromagneticField f=fieldFromOtherParticleAt(x,t,s,h,first,1.0e-5,matchedMomentSoftening());
          return dot(m,cross(v,f.electric))/(gg*c*c); };
        Vec3 g; const Vec3 e[3]={{st,0,0},{0,st,0},{0,0,st}};
        g.x=(U(x0+e[0])-U(x0-e[0]))/(2*st); g.y=(U(x0+e[1])-U(x0-e[1]))/(2*st); g.z=(U(x0+e[2])-U(x0-e[2]))/(2*st); return g; };
      const Vec3 tg1=thomasGrad(true);
      worstCheck=std::max(worstCheck,(fw1-grad1-tg1).norm()); scaleCheck=std::max(scaleCheck,tg1.norm());
      gDipoleForceEnabled=true; const MutualForces fOn=retardedExternalForces(s,h);
      gDipoleForceEnabled=false; const MutualForces fOff=retardedExternalForces(s,h); gDipoleForceEnabled=true;
      tauSector=tauSector+(cross(s.firstPosition,fOn.first-fOff.first)+cross(s.secondPosition,fOn.second-fOff.second))*wj;
      const Vec3 F1[4]={grad1+hid1+fL1, grad1+hid1+th.first+fL1, fw1+fL1, fw1+hid1+fL1};
      const Vec3 F2[4]={grad2+hid2+fL2, grad2+hid2+th.second+fL2, fw2+fL2, fw2+hid2+fL2};
      const Vec3 dC=fC1*(firstCharge/firstMass)+fC2*(secondCharge/secondMass);
      const DipoleDerivatives dm=thomasBmtDipoleDerivatives(s,h);
      dS=dS+(dm.first*(1.0/g1)+dm.second*(1.0/g2))*wj; ds1=ds1+dm.first*(wj/g1);
      for(int mdl=0;mdl<4;++mdl){
        tau[mdl]=tau[mdl]+(cross(s.firstPosition,F1[mdl])+cross(s.secondPosition,F2[mdl]))*wj;
        const Vec3 dD=F1[mdl]*(firstCharge/firstMass)+F2[mdl]*(secondCharge/secondMass);
        sX[mdl]+=2*dot(dC,dD)*wj; }
      sC+=dC.squaredNorm()*wj; sw+=wj;
      // Frenkel effective mass: M - m = -(q/(m c^2)) xi.B' (xi rest-frame spin, B' rest-frame field)
      const Vec3 Bp1=lw1.magnetic+md1.magnetic-cross(s.firstVelocity,lw1.electric+md1.electric)*(1/(c*c));
      maxMass=std::max(maxMass,std::abs(firstCharge/(firstMass*c*c)*dot(q.s1*(0.5*hbar),Bp1))/firstMass);
    }
    const char* name[4]={"A engine","B A+Thomas","C FW","D FW+hidden"};
    const Vec3 dSav=dS*(1.0/sw);
    std::printf("  [%s] <dS/dt> %.4e %.4e %.4e J (|.| T/hbar %.3e per orbit); max|M-m|/m %.3e; check |C-A_grad-Thomas grad| %.2e of %.2e N\n",
      q.n,dSav.x,dSav.y,dSav.z,dSav.norm()*T/hbar,maxMass,worstCheck,scaleCheck);
    {
      const Vec3 ts=tauSector*(1.0/sw); const double ds1n=(ds1*(1.0/sw)).norm();
      const bool prec=dSav.norm()>1e-3*ds1n&&dSav.norm()>0.0;
      std::printf("     %-12s <tau> %+.4e %+.4e %+.4e J  |<tau>+<dS/dt>| %.3e hbar/orbit  rel %.5f %s\n","PRODUCTION",
        ts.x,ts.y,ts.z,(ts+dSav).norm()*T/hbar,prec?(ts+dSav).norm()/dSav.norm():(ds1n>0?(ts+dSav).norm()/ds1n:-1.0),
        prec?"(/|dS|)":"(/|ds1|)");
    }
    for(int mdl=0;mdl<4;++mdl){
      const Vec3 ta=tau[mdl]*(1.0/sw);
      // mismatch of d(L+S)/dt in hbar per orbit, and relative to |<dS/dt>| when the spins precess
      const double ds1n=(ds1*(1.0/sw)).norm();
      const double perOrbit=(ta+dSav).norm()*T/hbar, rel=dSav.norm()>1e-3*ds1n&&dSav.norm()>0.0?(ta+dSav).norm()/dSav.norm()
          :(ds1n>0.0?(ta+dSav).norm()/ds1n:-1.0);
      std::printf("     %-12s <tau> %+.4e %+.4e %+.4e J  |<tau>+<dS/dt>| %.3e hbar/orbit  rel %s%s P_cross/P_E1 %+.4e\n",
        name[mdl],ta.x,ta.y,ta.z,perOrbit,rel<0?"   --  ":(std::to_string(rel).substr(0,7)).c_str(),
        (dSav.norm()>1e-3*ds1n&&dSav.norm()>0.0)?" (/|dS|) ":(ds1n>0.0?" (/|ds1|)":"        "),(kL*sX[mdl]/sw)/(kL*sC/sw));
    }
  }
}
