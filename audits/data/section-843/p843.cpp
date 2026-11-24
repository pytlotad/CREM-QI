// Audit 843 (static probe S; = p841 with forms ZD (CREM_SPIN_ORBIT_DIRAC, = 841 ZD) and ZDC (+ CREM_SPIN_ORBIT_CANONICAL)).
// Tree 8a856159 + switch.  (p841 header:) Audit 841 static probe S, = p840 with the engine form taken from allExternalForces.
// (p840 header:)  Hydrogen as a test particle (nucleus mass x 1e6, no proton dipole):
// engine forms X = Z (L_T), ZT (Tulczyjew), Z0, ZT0 (lab dipole = proper dipole) against the classical limit of the Dirac FW
// Hamiltonian  H_D = c sqrt(m^2c^2 + p^2) + q phi + S.Omega_BMT(v(p), E),  S' = Omega x S.
//   R_X = a^D(x, v) + xi'' - a^X(x + xi, v + xi'), xi = 0 (Z, Z0) or (S x v)/(2 m c^2) (ZT, ZT0);  Omega_X - Omega_D = -<LRL[R_X]>.
// argv: lfrac s2x s2y s2z(electron dipole in frame P, Q, Z=L) N hdiv mode(R | C)
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
static bool gNoBoost=false;
static StateHistory localHist(const State& s){ State l{s}; l.firstAcceleration={}; l.secondAcceleration={}; return StateHistory{l}; }
static MutualForces thomasT(const State& s,const DipoleDerivatives& derivatives){   // = p839 (copy of the 834 branch)
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
static MutualForces forces(const State& s,bool tform){   // = allExternalForces (814 measurement mode, no external field)
  MutualForces f=coulombForces(s);
  const auto add=[&](const MutualForces& g){ f.first=f.first+g.first; f.second=f.second+g.second; };
  add(darwinForces(s));
  const DipoleDerivatives dd=thomasBmtDipoleDerivatives(s,localHist(s));
  const Vec3 dip=pairDipoleForce(s.firstPosition-s.secondPosition,s.firstDipole,s.secondDipole); add({dip,dip*-1.0});
  const Vec3 fms=s.firstPosition-s.secondPosition, vrel=s.firstVelocity-s.secondVelocity;
  const ChargeDipolePairForces a=chargeDipolePairForces(vrel,firstCharge,s.secondDipole,dd.second,fms);
  const ChargeDipolePairForces b=chargeDipolePairForces(vrel*-1.0,secondCharge,s.firstDipole,dd.first,fms*-1.0);
  add({a.onCharge+b.onDipole,b.onCharge+a.onDipole});
  add(tform?thomasT(s,dd):thomasBackReactionForces(s,dd));
  return f;
}
static double M_, gyro2;
// engine state from the relative electron coordinate d = x_e - x_p, velocity v, electron proper dipole mu
static State engineState(const Vec3& d,const Vec3& v,const Vec3& mu){
  State s{}; const double M=firstMass+secondMass;
  s.firstPosition=d*(-secondMass/M); s.secondPosition=d*(firstMass/M); s.firstVelocity=v*(-secondMass/M); s.secondVelocity=v*(firstMass/M);
  s.secondProperDipole=mu; s.secondDipole=mu; if(!gNoBoost) synchronizeCovariantDipoles(s); return s; }
static Vec3 engineRelAcc(const Vec3& d,const Vec3& v,const Vec3& mu,bool dirac){
  const State s=engineState(d,v,mu); gSpinOrbitDirac=true; gSpinOrbitCanonical=dirac; const MutualForces F=allExternalForces(s); gSpinOrbitDirac=false; gSpinOrbitCanonical=false;
  return relativisticAcceleration(s.secondVelocity,F.second,secondMass)-relativisticAcceleration(s.firstVelocity,F.first,firstMass); }
// ---- Dirac classical limit (static nucleus at the origin; electron = second particle) ----
static Vec3 Efield(const Vec3& d){ return d*(coulomb*firstCharge/std::pow(d.norm(),3)); }
static Vec3 vOfP(const Vec3& p){ const double m=secondMass; return p*(c/std::sqrt(m*m*c*c+p.squaredNorm())); }
static Vec3 OmegaBMT(const Vec3& p,const Vec3& d){
  return thomasBmtEffectiveField(vOfP(p),ElectromagneticField{Efield(d),Vec3{}},secondGFactor)*(-secondCharge/secondMass); }
static double HS(const Vec3& d,const Vec3& p,const Vec3& S){ return dot(S,OmegaBMT(p,d)); }
static Vec3 gradP(const Vec3& d,const Vec3& p,const Vec3& S){
  const double dp=1e-6*p.norm(); Vec3 g; const Vec3 e[3]={{1,0,0},{0,1,0},{0,0,1}}; double* gc[3]={&g.x,&g.y,&g.z};
  for(int i=0;i<3;++i) *gc[i]=(HS(d,p+e[i]*dp,S)-HS(d,p-e[i]*dp,S))/(2*dp); return g; }
static Vec3 gradX(const Vec3& d,const Vec3& p,const Vec3& S){
  const double dx=1e-6*d.norm(); Vec3 g; const Vec3 e[3]={{1,0,0},{0,1,0},{0,0,1}}; double* gc[3]={&g.x,&g.y,&g.z};
  for(int i=0;i<3;++i) *gc[i]=(HS(d+e[i]*dx,p,S)-HS(d-e[i]*dx,p,S))/(2*dx); return g; }
struct Y { Vec3 d,p,S; };
static Vec3 velD(const Y& y){ return vOfP(y.p)+gradP(y.d,y.p,y.S); }
static Y rateD(const Y& y){
  const Vec3 F=Efield(y.d)*secondCharge;   // -grad(q phi)
  return {velD(y),F-gradX(y.d,y.p,y.S),cross(OmegaBMT(y.p,y.d),y.S)}; }
static Y axpy(const Y& z,const Y& k,double h){ return {z.d+k.d*h,z.p+k.p*h,z.S+k.S*h}; }
static Y rk4(const Y& z,double h){ const Y k1=rateD(z),k2=rateD(axpy(z,k1,h/2)),k3=rateD(axpy(z,k2,h/2)),k4=rateD(axpy(z,k3,h));
  Y r=z; r=axpy(r,k1,h/6); r=axpy(r,k2,h/3); r=axpy(r,k3,h/3); r=axpy(r,k4,h/6); return r; }
static Vec3 pOfV(const Vec3& d,const Vec3& v,const Vec3& S){
  const double m=secondMass; Vec3 p=v*(m/std::sqrt(1-v.squaredNorm()/(c*c)));
  for(int it=0;it<50;++it){ const Vec3 u=v-gradP(d,p,S); const Vec3 pn=u*(m/std::sqrt(1-u.squaredNorm()/(c*c)));
    if((pn-p).norm()<=1e-15*p.norm()){ p=pn; break; } p=pn; }
  return p; }
int main(int argc,char**argv){
  if(argc<8){ std::fprintf(stderr,"usage\n"); return 1; }
  applyPairFromOption("proton,electron");
  firstMass*=1e6; pairReducedMass=firstMass*secondMass/(firstMass+secondMass);   // test particle
  gMeasurementConservativeForce=true;
  const double lfrac=std::atof(argv[1]); const Vec3 mu2{std::atof(argv[2]),std::atof(argv[3]),std::atof(argv[4])};
  const int N=std::atoi(argv[5]); const double hdiv=std::atof(argv[6]); const std::string mode=argv[7];
  gyro2=secondGyromagneticRatioOf();
  const double m=secondMass, mu=pairReducedMass, k=pairCoulombStrength, nlev=2.0;
  const double a=nlev*nlev*hbar*hbar/(mu*k), ecc=std::sqrt(std::max(0.0,1-lfrac*lfrac)), ec=std::sqrt(1-ecc*ecc);
  const double nmot=std::sqrt(k/(mu*a*a*a)), T=2*pi/nmot, h=T/hdiv;
  const Vec3 P{1,0,0}, Q{0,1,0};
  const auto kep=[&](double E,Vec3& r,Vec3& v){ const double w=1-ecc*std::cos(E);
    r=P*(a*(std::cos(E)-ecc))+Q*(a*ec*std::sin(E)); v=(P*(-a*nmot*std::sin(E))+Q*(a*nmot*ec*std::cos(E)))*(1/w); };
  std::printf("SYS lfrac %.17g e %.12f a %.12e T %.12e N %d hdiv %.0f mu2 %.6e %.6e %.6e\n",lfrac,ecc,a,T,N,hdiv,mu2.x,mu2.y,mu2.z);
  if(mode=="C"){
    const bool tf=std::getenv("CREM_THOMAS_TULCZYJEW_INSTANT")!=nullptr; double mx=0;
    for(int j=0;j<N;++j){ Vec3 r,v; kep(2*pi*(j+0.5)/N,r,v); const State s=engineState(r,v,mu2);
      const MutualForces mine=forces(s,tf), eng=allExternalForces(s);
      mx=std::max(mx,std::max((mine.first-eng.first).norm(),(mine.second-eng.second).norm())/std::max(mine.first.norm(),mine.second.norm())); }
    std::printf("CHECK form %s maxrel %.3e\n",tf?"T":"L",mx); return 0; }
  const Vec3 S=mu2/gyro2;
  const char* names[2]={"ZD","ZDC"};
  double om[2]={0,0}, r2[2]={0,0}, f2=0, s0=0, sw=0, oms=0;
  for(int j=0;j<N;++j){
    const double E=2*pi*(j+0.5)/N, wj=1-ecc*std::cos(E); Vec3 r,v; kep(E,r,v);
    const Y z0{r,pOfV(r,v,S),S}; const Y zp1=rk4(z0,h), zp2=rk4(zp1,h), zm1=rk4(z0,-h), zm2=rk4(zm1,-h);
    const Y* zz[5]={&zm2,&zm1,&z0,&zp1,&zp2}; Vec3 V[5],XI[5];
    for(int q=0;q<5;++q){ V[q]=velD(*zz[q]); XI[q]=cross(zz[q]->S,V[q])*(1.0/(2.0*m*c*c)); }
    const auto d1=[&](const Vec3* u){ return (u[0]-u[1]*8.0+u[3]*8.0-u[4])*(1.0/(12*h)); };
    const auto d2=[&](const Vec3* u){ return (u[0]*-1.0+u[1]*16.0-u[2]*30.0+u[3]*16.0-u[4])*(1.0/(12*h*h)); };
    const Vec3 vD=V[2], aD=d1(V), xi=XI[2], xid=d1(XI), xidd=d2(XI);
    // spinless control and first-order scale
    gNoBoost=false; const Vec3 aZ=engineRelAcc(r,vD,mu2,false), aZ0=engineRelAcc(r,vD,Vec3{},false);
    f2+=(aZ-aZ0).squaredNorm()*wj;
    { const Y w{r,pOfV(r,vD,Vec3{}),Vec3{}}; const Y a1=rk4(w,h),a2=rk4(a1,h),b1=rk4(w,-h),b2=rk4(b1,-h);
      const Vec3 VV[5]={velD(b2),velD(b1),velD(w),velD(a1),velD(a2)}; s0+=(d1(VV)-aZ0).squaredNorm()*wj;
      const Vec3 Rs=d1(VV)-aZ0, hs=cross(r,vD), dAs=cross(Rs,hs)+cross(vD,cross(r,Rs)); oms+=dot(Q,dAs)/(k*ecc/mu)*wj; }   // POST HOC: secular rate of the spinless residual
    for(int x=0;x<2;++x){
      gNoBoost=false;
      const Vec3 Rr=aD-engineRelAcc(r,vD,mu2,x==1);
      const Vec3 hv=cross(r,vD), dA=cross(Rr,hv)+cross(vD,cross(r,Rr));
      om[x]+=dot(Q,dA)/(k*ecc/mu)*wj; r2[x]+=Rr.squaredNorm()*wj; }
    sw+=wj;
  }
  std::printf("SCALE rmsFirst %.6e rmsSpinless %.6e\n",std::sqrt(f2/sw),std::sqrt(s0/sw));
  if(std::getenv("P840_POSTHOC")) std::printf("SPINLESS OmegaR %.12e\n",oms/sw);
  for(int x=0;x<2;++x) std::printf("OMX %s OmegaR %.12e rmsR %.6e\n",names[x],om[x]/sw,std::sqrt(r2[x]/sw));
}
