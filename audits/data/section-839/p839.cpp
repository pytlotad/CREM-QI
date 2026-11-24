// Audit 839 (static probe, no trajectory integration beyond +-2h around each orbit point; engine unchanged, tree facf06c2).
// D = Omega(ZT) - Omega(Z) from the exact conjugacy residual: y_i = x_i + xi_i, xi_i = (S_i x v_i)/(2 m_i c^2),
//   R_i = a_i^L(X) + xi_i'' - a_i^T(x + xi, v + xi'),   D = - <rotation rate of the Laplace-Runge-Lenz vector under R_1 - R_2>.
// a^L: the engine's measurement sum (allExternalForces, 814 conservative mode) with thomasBackReactionForces (L_T form);
// a^T: the same sum with the Thomas term replaced by -2 potentialRate on the spin carrier (copy of the 834 branch).
// Ingredient sets (probe composes the sum itself): bit 1 gamma kinematics, 2 Darwin, 4 dipole-dipole + charge-dipole,
// 8 precession (flow of the proper dipoles and spinRate in Thomas / induced term).  15 = full = allExternalForces.
// argv: pair rp e m1x m1y m1z m2x m2y m2z N hdiv mode   (m: dipole moments in the frame P = periapsis, Q, Z = L; mode R | C)
//   mode R: prints OMR for sets 0,1,3,7,15;  mode C: max relative difference of the set-15 sum (L form, or T form when
//   CREM_THOMAS_TULCZYJEW_INSTANT is set) against allExternalForces at the N orbit points.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
static bool gThomasT=false;
static Vec3 nz(){ return Vec3{}; }
static StateHistory localHist(const State& s){ State l{s}; l.firstAcceleration={}; l.secondAcceleration={}; return StateHistory{l}; }
// Thomas term, both forms (L = engine function; T = copy of the 834 branch of thomasBackReactionForces).
static MutualForces thomasT(const State& s,const DipoleDerivatives& derivatives){
  const Vec3 fms=s.firstPosition-s.secondPosition; const double distance=fms.norm(); const double floor=separationFloor();
  double radialFactor=0.0,gradientFactor=0.0;
  if(floor>0.0){ const double rho2=distance*distance+floor*floor; radialFactor=1.0/std::pow(rho2,1.5); gradientFactor=-3.0/std::pow(rho2,2.5); }
  else { const MagneticRadialProfile p=magneticRadialProfile(distance); radialFactor=p.vectorPotentialFactor; gradientFactor=p.firstDerivative/distance; }
  const auto contribution=[&](bool spinFirst){
    const double g=spinFirst?firstGyromagneticRatioOf():secondGyromagneticRatioOf();
    if(!(std::abs(g)>0.0)) return MutualForces{};
    const double qp=spinFirst?firstCharge:secondCharge, mp=spinFirst?firstMass:secondMass, qo=spinFirst?secondCharge:firstCharge;
    const Vec3 vp=spinFirst?s.firstVelocity:s.secondVelocity, vo=spinFirst?s.secondVelocity:s.firstVelocity;
    const Vec3 spin=(spinFirst?s.firstDipole:s.secondDipole)/g, spinRate=(spinFirst?derivatives.first:derivatives.second)/g;
    const Vec3 d=spinFirst?fms:fms*-1.0; const double coef=coulomb*qo*qp/(2.0*mp*c*c);
    const Vec3 vr=vp-vo;
    const Vec3 pr=(cross(d,spin)*(gradientFactor*dot(d,vr))+cross(vr,spin)*radialFactor+cross(d,spinRate)*radialFactor)*coef;
    const Vec3 on=pr*-2.0;
    return spinFirst?MutualForces{on,Vec3{}}:MutualForces{Vec3{},on}; };
  const MutualForces a=contribution(true), b=contribution(false);
  return {a.first+b.first,a.second+b.second};
}
static int gSet=15;
static DipoleDerivatives derivs(const State& s){ if(!(gSet&8)) return {}; return thomasBmtDipoleDerivatives(s,localHist(s)); }
static MutualForces forces(const State& s,bool tform){
  MutualForces f=coulombForces(s);
  const auto add=[&](const MutualForces& g){ f.first=f.first+g.first; f.second=f.second+g.second; };
  if(gSet&2) add(darwinForces(s));
  const DipoleDerivatives dd=derivs(s);
  if(gSet&4){
    const Vec3 dip=pairDipoleForce(s.firstPosition-s.secondPosition,s.firstDipole,s.secondDipole); add({dip,dip*-1.0});
    const Vec3 fms=s.firstPosition-s.secondPosition, vrel=s.firstVelocity-s.secondVelocity;
    const ChargeDipolePairForces a=chargeDipolePairForces(vrel,firstCharge,s.secondDipole,dd.second,fms);
    const ChargeDipolePairForces b=chargeDipolePairForces(vrel*-1.0,secondCharge,s.firstDipole,dd.first,fms*-1.0);
    add({a.onCharge+b.onDipole,b.onCharge+a.onDipole});
  }
  add(tform?thomasT(s,dd):thomasBackReactionForces(s,dd));
  return f;
}
static MutualForces accel(const State& s,bool tform){
  const MutualForces F=forces(s,tform);
  if(gSet&1) return {relativisticAcceleration(s.firstVelocity,F.first,firstMass),relativisticAcceleration(s.secondVelocity,F.second,secondMass)};
  return {F.first*(1.0/firstMass),F.second*(1.0/secondMass)};
}
struct X { Vec3 x1,x2,v1,v2,p1,p2; };
static State toState(const X& z){ State s{}; s.firstPosition=z.x1; s.secondPosition=z.x2; s.firstVelocity=z.v1; s.secondVelocity=z.v2;
  s.firstProperDipole=z.p1; s.secondProperDipole=z.p2; s.firstDipole=z.p1; s.secondDipole=z.p2; if(!std::getenv("P839_NOBOOST")) synchronizeCovariantDipoles(s); return s; }   // POST HOC (839): P839_NOBOOST=1 lab dipole = proper dipole
static X rate(const X& z){
  const State s=toState(z); const MutualForces a=accel(s,false);
  X d{z.v1,z.v2,a.first,a.second,Vec3{},Vec3{}};
  if(gSet&8){ State q=s; q.firstDipole=z.p1; q.secondDipole=z.p2; const DipoleDerivatives dd=thomasBmtDipoleDerivatives(q,localHist(q)); d.p1=dd.first; d.p2=dd.second; }
  return d; }
static X axpy(const X& z,const X& d,double h){ return {z.x1+d.x1*h,z.x2+d.x2*h,z.v1+d.v1*h,z.v2+d.v2*h,z.p1+d.p1*h,z.p2+d.p2*h}; }
static X rk4(const X& z,double h){ const X k1=rate(z),k2=rate(axpy(z,k1,h/2)),k3=rate(axpy(z,k2,h/2)),k4=rate(axpy(z,k3,h));
  X r=z; const double w[4]={1,2,2,1}; const X* k[4]={&k1,&k2,&k3,&k4}; for(int i=0;i<4;++i) r=axpy(r,*k[i],h*w[i]/6.0); return r; }
static void xis(const X& z,Vec3& a,Vec3& b){
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  a=std::abs(g1)>0?cross(z.p1/g1,z.v1)*(1.0/(2.0*firstMass*c*c)):Vec3{};
  b=std::abs(g2)>0?cross(z.p2/g2,z.v2)*(1.0/(2.0*secondMass*c*c)):Vec3{}; }
int main(int argc,char**argv){
  if(argc<13){ std::fprintf(stderr,"usage\n"); return 1; }
  applyPairFromOption(argv[1]);
  gMeasurementConservativeForce=true;
  const double rp=std::atof(argv[2]), ecc=std::atof(argv[3]);
  const Vec3 m1{std::atof(argv[4]),std::atof(argv[5]),std::atof(argv[6])}, m2{std::atof(argv[7]),std::atof(argv[8]),std::atof(argv[9])};
  const int N=std::atoi(argv[10]); const double hdiv=std::atof(argv[11]); const std::string mode=argv[12];
  const double M=firstMass+secondMass, mu=firstMass*secondMass/M, k=pairCoulombStrength, a=rp/(1-ecc);
  const double nmot=std::sqrt(k/(mu*a*a*a)), T=2*pi/nmot, ec=std::sqrt(1-ecc*ecc), h=T/hdiv;
  const Vec3 P{1,0,0}, Q{0,1,0};
  const auto kepler=[&](double E){ const double w=1-ecc*std::cos(E);
    const Vec3 r=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E)), v=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w);
    return X{r*(secondMass/M),r*(-firstMass/M),v*(secondMass/M),v*(-firstMass/M),m1,m2}; };
  std::printf("SYS pair %s rp %.12e e %.12f a %.12e T %.12e N %d hdiv %.0f\n",argv[1],rp,ecc,a,T,N,hdiv);
  if(mode=="C"){
    gSet=15; double mx=0;
    for(int j=0;j<N;++j){ const double E=2*pi*(j+0.5)/N; const State s=toState(kepler(E));
      const MutualForces mine=forces(s,gThomasT=std::getenv("CREM_THOMAS_TULCZYJEW_INSTANT")!=nullptr), eng=allExternalForces(s);
      const double sc=std::max(mine.first.norm(),mine.second.norm());
      mx=std::max(mx,std::max((mine.first-eng.first).norm(),(mine.second-eng.second).norm())/sc); }
    std::printf("CHECK form %s maxrel %.3e\n",gThomasT?"T":"L",mx); return 0; }
  const int sets[5]={0,1,3,7,15};
  for(int si=0;si<5;++si){
    gSet=sets[si]; double sw=0, om=0, r2=0, f2=0;
    for(int j=0;j<N;++j){
      const double E=2*pi*(j+0.5)/N, wj=1-ecc*std::cos(E);
      const X z0=kepler(E);
      X zp1=rk4(z0,h), zp2=rk4(zp1,h), zm1=rk4(z0,-h), zm2=rk4(zm1,-h);
      Vec3 A[5],B[5]; const X* zz[5]={&zm2,&zm1,&z0,&zp1,&zp2}; for(int q=0;q<5;++q) xis(*zz[q],A[q],B[q]);
      const auto d1=[&](const Vec3* u){ return (u[0]-u[1]*8.0+u[3]*8.0-u[4])*(1.0/(12*h)); };
      const auto d2=[&](const Vec3* u){ return (u[0]*-1.0+u[1]*16.0-u[2]*30.0+u[3]*16.0-u[4])*(1.0/(12*h*h)); };
      const Vec3 xd1=d1(A), xd2=d1(B), xdd1=d2(A), xdd2=d2(B);
      const State s0=toState(z0); const MutualForces aL=accel(s0,false), aT0=accel(s0,true);
      X y=z0; y.x1=z0.x1+A[2]; y.x2=z0.x2+B[2]; y.v1=z0.v1+xd1; y.v2=z0.v2+xd2;
      const MutualForces aT=accel(toState(y),true);
      const Vec3 R1=aL.first+xdd1-aT.first, R2=aL.second+xdd2-aT.second, Rr=R1-R2;
      const Vec3 f1=(aT0.first-aL.first)-(aT0.second-aL.second);
      const Vec3 r=z0.x1-z0.x2, v=z0.v1-z0.v2, hv=cross(r,v);
      const Vec3 dA=cross(Rr,hv)+cross(v,cross(r,Rr));
      om+=dot(Q,dA)/(k*ecc/mu)*wj; r2+=Rr.squaredNorm()*wj; f2+=f1.squaredNorm()*wj; sw+=wj;
    }
    std::printf("OMR set %d OmegaR %.12e rmsR %.6e rmsFirst %.6e\n",sets[si],om/sw,std::sqrt(r2/sw),std::sqrt(f2/sw));
    std::fflush(stdout);
  }
}
