#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 571: do the model's fields satisfy Maxwell's equations away from the
// sources?  The rotating pair of tools/rotating_pair_momentum_balance (moments
// along the axis, exact circular kinematics in the history); for the fields of
// particle 1 -- charge (Lienard-Wiechert) and moment (retardedMagneticDipoleField)
// separately -- the Faraday residual |curl E + dB/dt| / |curl E| and the
// Ampere-Maxwell residual |curl B - dE/dt/c^2| / |curl B| at points on a sphere of
// radius 0.5 sep around particle 1 and on a ring of radius 2 sep about the centre;
// residuals normalised by the full spatial gradient |grad E|, |grad B| (v2).
// argv: beta separation[r*]
int main(int argc,char**argv){
  if(const char* pair=std::getenv("PAIR")) applyPairFromOption(pair);
  const double beta=std::atof(argv[1]), separation=std::atof(argv[2])*comptonBarrierRadius;
  const int nodes=argc>3?std::atoi(argv[3]):512; const double span=argc>4?std::atof(argv[4]):40.0;
  const Vec3 axis{0,0,1}; const double arm=0.5*separation, w=beta*c/arm;
  const auto build=[&](double t){ State s{}; s.time=t; const double ph=w*t;
    const Vec3 o{arm*std::cos(ph),arm*std::sin(ph),0}, v{-arm*w*std::sin(ph),arm*w*std::cos(ph),0};
    s.firstPosition=o; s.secondPosition=o*-1.0; s.firstVelocity=v; s.secondVelocity=v*-1.0;
    s.firstAcceleration=o*(-w*w); s.secondAcceleration=o*(w*w);
    s.firstProperDipole=axis*firstMagneticMoment; s.secondProperDipole=axis*secondMagneticMoment;
    synchronizeCovariantDipoles(s); return s; };
  const double lt=separation/c;
  StateHistory fixedHistory; for(int i=0;i<nodes;++i) fixedHistory.push_back(build(-span*lt*(1.0-double(i)/nodes)));
  const State fixedPresent=build(0.0);
  const bool fixedMode=std::getenv("FIXED_HISTORY")!=nullptr;   // v3 (audit 571)
  const auto fieldAt=[&](const Vec3& x,double t,int part){
    StateHistory h; State s;
    if(fixedMode){ h=fixedHistory; s=fixedPresent; t+=-1.0e-4*separation/c; }
    else { for(int i=0;i<nodes;++i) h.push_back(build(t-span*lt*(1.0-double(i)/nodes))); s=build(t); }
    if(part==0) return lienardWiechertField(x,t,h,s,true,firstCharge);
    if(part==1) return retardedMagneticDipoleField(x,t,h,s,true);
    // parts of retardedMagneticDipoleFieldExact (copied): 2 = dual poles, 3 = moving magnetization
    const bool sourceIsFirst=true;
    if(part==2){
      const ElectromagneticField dual=twoChargeLimitDipoleField(x,t,h,s,sourceIsFirst,1.0e-5,magneticDipoleRadius(),
        [&](double time,Vec3& moment,Vec3& first,Vec3& second){
          const RetardedElectricDipoleKinematics dipole=historicalIntegratedDipoleKinematics(h,s,sourceIsFirst,time,false);
          moment=dipole.moment/(c*c); first=dipole.firstDerivative/(c*c); second=dipole.secondDerivative/(c*c); });
      return ElectromagneticField{dual.magnetic*(-c*c),dual.electric}; }
    if(part==4){   // the electric-dipole field before the dual map (p = m/c^2)
      return twoChargeLimitDipoleField(x,t,h,s,sourceIsFirst,1.0e-5,magneticDipoleRadius(),
        [&](double time,Vec3& moment,Vec3& first,Vec3& second){
          const RetardedElectricDipoleKinematics dipole=historicalIntegratedDipoleKinematics(h,s,sourceIsFirst,time,false);
          moment=dipole.moment/(c*c); first=dipole.firstDerivative/(c*c); second=dipole.secondDerivative/(c*c); }); }
    if(part==5){   // the same construction with a CONSTANT moment (no moment kinematics from the history)
      const Vec3 mc=s.firstDipole*(1.0/(c*c));
      return twoChargeLimitDipoleField(x,t,h,s,sourceIsFirst,1.0e-5,magneticDipoleRadius(),
        [&](double,Vec3& moment,Vec3& first,Vec3& second){ moment=mc; first={}; second={}; }); }
    double mt=t; ChargeKinematics src=historicalCharge(h,s,sourceIsFirst,mt);
    for(int it=0;it<16;++it){ const Vec3 dd=x-src.position; const double dist=dd.norm();
      const double res=mt+dist/c-t; const double der=std::max(1.0e-8,1.0-dot(dd/std::max(dist,1e-300),src.velocity/c));
      const double ref=mt-res/der; const bool conv=std::abs(ref-mt)<=1e-30+1e-14*std::abs(mt); mt=ref; src=historicalCharge(h,s,sourceIsFirst,mt); if(conv) break; }
    const State rs=historicalState(h,s,mt); Vec3 pm=rs.firstProperDipole; if(pm.squaredNorm()==0.0) pm=s.firstProperDipole;
    const double lag=t-mt; const Vec3 xp=src.position+src.velocity*lag+src.acceleration*(0.5*lag*lag);
    const double sg=1.0/std::sqrt(std::max(1.0-src.velocity.squaredNorm()/(c*c),1e-300));
    const Vec3 u4=src.velocity*sg+(src.acceleration*sg+src.velocity*(sg*sg*sg*dot(src.velocity,src.acceleration)/(c*c)))*lag;
    const Vec3 vp=u4/std::sqrt(1.0+u4.squaredNorm()/(c*c));
    return movingMagnetizationField(x-xp,vp,pm,magneticDipoleRadius()); };
  const double hx=1e-4*separation, ht=hx/c;
  const auto residuals=[&](const Vec3& x,int part,double& far,double& amp){
    ElectromagneticField fx[3][2];
    const Vec3 e[3]={{hx,0,0},{0,hx,0},{0,0,hx}};
    for(int k=0;k<3;++k){ fx[k][0]=fieldAt(x+e[k],0.0,part); fx[k][1]=fieldAt(x-e[k],0.0,part); }
    const auto d=[&](int k,bool electric){ return ((electric?fx[k][0].electric:fx[k][0].magnetic)-(electric?fx[k][1].electric:fx[k][1].magnetic))*(1.0/(2*hx)); };
    const Vec3 dEx=d(0,true),dEy=d(1,true),dEz=d(2,true),dBx=d(0,false),dBy=d(1,false),dBz=d(2,false);
    const Vec3 curlE{dEy.z-dEz.y,dEz.x-dEx.z,dEx.y-dEy.x}, curlB{dBy.z-dBz.y,dBz.x-dBx.z,dBx.y-dBy.x};
    const ElectromagneticField fp=fieldAt(x,ht,part), fm=fieldAt(x,-ht,part);
    const Vec3 dB=(fp.magnetic-fm.magnetic)*(1.0/(2*ht)), dE=(fp.electric-fm.electric)*(1.0/(2*ht));
    // normalised by the full spatial gradient of the field (curl alone vanishes for static fields)
    const double gE=std::sqrt(dEx.squaredNorm()+dEy.squaredNorm()+dEz.squaredNorm());
    const double gB=std::sqrt(dBx.squaredNorm()+dBy.squaredNorm()+dBz.squaredNorm());
    far=(curlE+dB).norm()/std::max(gE,1e-300);
    amp=(curlB-dE*(1.0/(c*c))).norm()/std::max(gB,1e-300); };
  const State s0=build(0.0);
  const int partMax=std::getenv("PARTS6")?6:4;
  for(int part=0;part<partMax;++part){
    double fs=0,as=0,fmx=0,amx=0; int n=0;
    for(int i=0;i<64;++i){
      Vec3 x;
      if(i<32){ const double ct=1-2*(i+0.5)/32, st=std::sqrt(1-ct*ct), ph=2.399963*i;
        x=s0.firstPosition+Vec3{st*std::cos(ph),st*std::sin(ph),ct}*(0.5*separation); }
      else { const double ph=2*pi*(i-32+0.5)/32; x=Vec3{std::cos(ph),std::sin(ph),0.1}*(2.0*separation); }
      double f,a; residuals(x,part,f,a); fs+=f*f; as+=a*a; fmx=std::max(fmx,f); amx=std::max(amx,a); ++n; }
    std::printf("beta %.2f sep %.0f r*  %-6s  Faraday rms %.3e max %.3e   Ampere-Maxwell rms %.3e max %.3e\n",beta,separation/comptonBarrierRadius,
      part==0?"charge":part==1?"moment":part==2?"poles":part==3?"magnet":part==4?"p-dip":"p-const",std::sqrt(fs/n),fmx,std::sqrt(as/n),amx);
  }
}
