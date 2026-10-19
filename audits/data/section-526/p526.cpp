#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 526: (C') heavy-spin precession averaged over the fine-structure cone about the engine's J (hydrogen 2P1/2), component along J;
// (D) the secular transport itself over 30 ns: rotation frequency of I^ about J^_e and the spread of I^.J^_e.  Run in three modes (env).
int main(){
  applyPairFromOption("proton,electron");
  const double aB=pairBohrRadius(activePair), a=4*aB, L=1.5*hbar, mu=firstMass*secondMass/(firstMass+secondMass);
  std::printf("mode: %s\n",!std::getenv("CREM_NO_HFS_ADIABATIC")?"adiabatic (526)":(std::getenv("CREM_RULE_DIPOLE_TENSOR")?"rule tensor (525)":"classical"));
  const double cLS=-0.75, Jn=std::sqrt(2.25+0.25+2*1.5*0.5*cLS);
  const double aL=std::acos((1.5+0.5*cLS)/Jn), as=std::acos((0.5+1.5*cLS)/Jn);
  { // (C')
    const int N=72; Vec3 sum;
    for(int k=0;k<N;++k){ const double ph=2*pi*k/N;
      const Vec3 Lh{std::sin(aL)*std::cos(ph),std::sin(aL)*std::sin(ph),std::cos(aL)}, sh{-std::sin(as)*std::cos(ph),-std::sin(as)*std::sin(ph),std::cos(as)};
      const Vec3 me=sh*(-secondMagneticMoment), mp=Vec3{1,0,0}*firstMagneticMoment;
      Vec3 per=cross(Lh,Vec3{0,1,0}); per=per*(1.0/per.norm());
      OrbitAveragedBmtAngularVelocities r=secularRatesWithSStateRule(a,Lh*L,mp,me,mu,0.0,per);
      applyContactAboveHbar(r,a,Lh*L,mp,me,mu); applyBornInverseCube(r,Lh*L); applyRuleDipoleTensor(r,a,Lh*L,mp,me); applyAdiabaticHyperfine(r,a,Lh*L,mp,me);
      sum=sum+r.first*(1.0/N); }
    std::printf("(C') cone about the engine's J (L^ %.2f deg, s^ %.2f deg; |L + s| = %.4f hbar): proton precession along J %.4f MHz (import 22: 59.1687; measured 59.22(14))\n",
      aL*180/pi,as*180/pi,Jn,std::abs(sum.z)/(2*pi)*1e-6);
  }
  { // (D) transport dynamics, 30 ns
    SecularSpinOrbitState s; const Vec3 Lh{std::sin(aL),0,std::cos(aL)}, sh{-std::sin(as),0,std::cos(as)};
    s.orbitalAngularMomentum=Lh*L; s.secondDipole=sh*(-secondMagneticMoment);
    const Vec3 ip=Vec3{0.6,0.3,std::sqrt(1-0.45)}; s.firstDipole=ip*firstMagneticMoment;
    Vec3 per=cross(Lh,Vec3{0,1,0}); s.periapsisDirection=per*(1.0/per.norm());
    const double dt=2.0e-11, T=30e-9; double angle=0, angleF=0, cmin=2, cmax=-2; Vec3 prevPerp, prevPerpF; bool have=false, haveF=false; int Fchanges=0, Fprev=-1;
    for(double t=0;t<T;t+=dt){
      const SecularSpinOrbitAdvance ad=advanceCoupledSecularSpinOrbit(s,a,mu,dt); s=ad.state;
      const Vec3 se=s.secondDipole*(-1.0/s.secondDipole.norm()), Je=s.orbitalAngularMomentum+se*(0.5*hbar), Jh=Je*(1.0/Je.norm());
      const Vec3 Ih=s.firstDipole*(1.0/s.firstDipole.norm()); const double c=dot(Ih,Jh); cmin=std::min(cmin,c); cmax=std::max(cmax,c);
      const Vec3 Lq=s.orbitalAngularMomentum*(1.5/s.orbitalAngularMomentum.norm()*1.0); const Vec3 Jq=Lq+se; const int k=std::clamp((int)std::lround(Jq.norm()),1,2);
      const int F=std::clamp((int)std::lround((Jq*(k/Jq.norm())+Ih).norm()-0.5),(int)std::lround(std::abs(k-1.0)),k);
      if(Fprev>=0&&F!=Fprev) ++Fchanges; Fprev=F;
      { const Vec3 Fv=Je+Ih*(0.5*hbar), Fh=Fv*(1.0/Fv.norm()); Vec3 pf=Ih-Fh*dot(Ih,Fh);
        if(pf.norm()>1e-9){ pf=pf*(1.0/pf.norm()); if(haveF) angleF+=std::atan2(dot(cross(prevPerpF,pf),Fh),dot(prevPerpF,pf)); prevPerpF=pf; haveF=true; } }
      Vec3 perp=Ih-Jh*c; if(perp.norm()>1e-9){ perp=perp*(1.0/perp.norm()); if(have) angle+=std::atan2(dot(cross(prevPerp,perp),Jh),dot(prevPerp,perp)); prevPerp=perp; have=true; }
    }
    std::printf("(D) transport 30 ns: I^ about F^_e = (J_e + I_e)^ %.4f MHz (the beat; import 22 59.1687, measured 59.22(14)); about J^_e %.4f MHz; I^.J^_e in [%+.4f, %+.4f], F-basin changes %d\n",
      std::abs(angleF)/(2*pi*T)*1e-6,std::abs(angle)/(2*pi*T)*1e-6,cmin,cmax,Fchanges);
  }
}
