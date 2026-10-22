#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
// Audit 570 (from p568): the dipole sector WITHOUT the engine's hidden-momentum
// terms and without the Thomas back-reaction (run with CREM_NO_HIDDEN_MOMENTUM_FORCE=1
// CREM_NO_RETARDED_THOMAS=1), plus a hidden-momentum force computed here:
//   lab form (control): -d/dt [m_lab x (E + v x B)]/c^2, as the engine since 566;
//   covariant: p' = m_0 x E'/c^2 in the particle's rest frame (exact E'), boosted
//   as a four-vector, force = -dp/dtau projected orthogonal to u, spatial / gamma.
// Audit 568: covariance of the engine's retarded forces.  A pair on a
// circular orbit in its own frame S' (moments along the orbit axis z); the lab
// S sees the pair move with V along y.  For the event of particle 1 at a given
// orbit phase the four-force K = gamma_u (F.u/c, F) from S' is boosted to S and
// compared with the engine's own lab force at the same event.  Whole force,
// dipole sector (gDipoleForceEnabled on - off) and charges alone (control).
// argv: pair  separation[r*]  V/c  [nodes span(d/c)]
int main(int argc,char**argv){
  applyPairFromOption(argv[1]);
  const double d=std::atof(argv[2])*comptonBarrierRadius, V=std::atof(argv[3])*c;
  const int nodes=argc>4?std::atoi(argv[4]):512; const double span=argc>5?std::atof(argv[5]):40.0;
  const double M=firstMass+secondMass, mu=firstMass*secondMass/M, k=pairCoulombStrength;
  const double w=std::sqrt(k/(mu*d*d*d)), R1=d*secondMass/M, R2=d*firstMass/M;
  const double G=1.0/std::sqrt(1.0-V*V/(c*c));
  const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
  const Vec3 Z{0,0,1};
  // rest-frame worldlines
  const auto restX=[&](int i,double t)->Vec3{ const double R=i==1?R1:-R2; return Vec3{R*std::cos(w*t),R*std::sin(w*t),0}; };
  const auto restU=[&](int i,double t)->Vec3{ const double R=i==1?R1:-R2; return Vec3{-R*w*std::sin(w*t),R*w*std::cos(w*t),0}; };
  const auto restA=[&](int i,double t)->Vec3{ return restX(i,t)*(-w*w); };
  const auto restState=[&](double t){
    State s{}; s.time=t;
    s.firstPosition=restX(1,t); s.secondPosition=restX(2,t);
    s.firstVelocity=restU(1,t); s.secondVelocity=restU(2,t);
    s.firstAcceleration=restA(1,t); s.secondAcceleration=restA(2,t);
    s.firstProperDipole=Z*(g1*0.5*hbar); s.secondProperDipole=Z*(g2*0.5*hbar);
    synchronizeCovariantDipoles(s); return s; };
  // lab worldline of particle i at lab time t
  const auto restTimeOf=[&](int i,double t){
    double tp=t/G-V*restX(i,t/G).y/(c*c);
    for(int it=0;it<50;++it){ const double g=G*(tp+V*restX(i,tp).y/(c*c))-t; const double dg=G*(1+V*restU(i,tp).y/(c*c));
      tp-=g/dg; if(std::abs(g)<1e-30) break; }
    return tp; };
  const auto labX=[&](int i,double t)->Vec3{ const double tp=restTimeOf(i,t); const Vec3 x=restX(i,tp); return Vec3{x.x,G*(x.y+V*tp),x.z}; };
  const auto labU=[&](int i,double t)->Vec3{ const double tp=restTimeOf(i,t); const Vec3 u=restU(i,tp); const double den=1+V*u.y/(c*c);
    return Vec3{u.x/(G*den),(u.y+V)/den,u.z/(G*den)}; };
  const auto labState=[&](double t){
    State s{}; s.time=t;
    s.firstPosition=labX(1,t); s.secondPosition=labX(2,t);
    s.firstVelocity=labU(1,t); s.secondVelocity=labU(2,t);
    const double h=1e-4/w;
    s.firstAcceleration=(labU(1,t+h)-labU(1,t-h))*(1/(2*h)); s.secondAcceleration=(labU(2,t+h)-labU(2,t-h))*(1/(2*h));
    s.firstProperDipole=Z*(g1*0.5*hbar); s.secondProperDipole=Z*(g2*0.5*hbar);
    synchronizeCovariantDipoles(s); return s; };
  const double lt=d/c;
  // hidden momentum of particle i at time t of a frame (four-vector: p0 (as p0*c? no: momentum units), spatial)
  struct H4{ double p0; Vec3 p; };
  const auto hiddenAt=[&](auto build,double t,int i,bool covariant)->H4{
    StateHistory h; for(int k=0;k<nodes;++k) h.push_back(build(t-span*lt*(1.0-double(k)/nodes)));
    const State s=build(t);
    const bool first=i==1;
    const Vec3 x=first?s.firstPosition:s.secondPosition, v=first?s.firstVelocity:s.secondVelocity;
    const ElectromagneticField f=fieldFromOtherParticleAt(x,t,s,h,first,1.0e-5,matchedMomentSoftening());
    if(!covariant){ const Vec3 m=first?s.firstDipole:s.secondDipole;
      return H4{0.0,cross(m,f.electric+cross(v,f.magnetic))*(1.0/(c*c))}; }
    const Vec3 m0=first?s.firstProperDipole:s.secondProperDipole;
    const double gu=1.0/std::sqrt(1-v.squaredNorm()/(c*c));
    const Vec3 Ep=(f.electric+cross(v,f.magnetic))*gu-v*(gu*gu/(gu+1)*dot(v,f.electric)/(c*c));
    const Vec3 pp=cross(m0,Ep)*(1.0/(c*c));
    const double vn=v.norm(); const Vec3 vh=vn>0?v*(1.0/vn):Vec3{};
    return H4{gu*dot(v,pp)/c, pp+vh*((gu-1)*dot(pp,vh))}; };
  // the force from it on particle i at time t (spatial), lab form unprojected (as the engine), covariant projected
  const auto hiddenForce=[&](auto build,double t,int i,bool covariant)->Vec3{
    const double hh=1e-4/w;
    const H4 a=hiddenAt(build,t+hh,i,covariant), b=hiddenAt(build,t-hh,i,covariant);
    const Vec3 dps=(a.p-b.p)*(1.0/(2*hh)); const double dp0=(a.p0-b.p0)/(2*hh);
    if(!covariant) return dps*-1.0;
    const State s=build(t); const Vec3 v=i==1?s.firstVelocity:s.secondVelocity;
    const double gu=1.0/std::sqrt(1-v.squaredNorm()/(c*c));
    // K = -gamma dp/dt; u = gamma (c, v); project orthogonal to u
    const double K0=-gu*dp0; const Vec3 Ks=dps*(-gu);
    const double Ku=K0*gu*c-dot(Ks,v*gu);
    const Vec3 Kperp=Ks-v*(gu*Ku/(c*c));
    return Kperp*(1.0/gu); };
  const bool covariantMode=std::getenv("COVARIANT")!=nullptr;
  // stage 2: Thomas momentum p'_T = kappa q (E' x S_0)/(m c^2) in the particle rest frame, as a four-vector
  const double kappa=std::getenv("THOMAS_KAPPA")?std::atof(std::getenv("THOMAS_KAPPA")):0.0;
  const auto thomasAt=[&](auto build,double t,int i)->H4{
    StateHistory h; for(int k=0;k<nodes;++k) h.push_back(build(t-span*lt*(1.0-double(k)/nodes)));
    const State s=build(t); const bool first=i==1;
    const Vec3 x=first?s.firstPosition:s.secondPosition, v=first?s.firstVelocity:s.secondVelocity;
    const ElectromagneticField f=fieldFromOtherParticleAt(x,t,s,h,first,1.0e-5,0.0);
    const double gu=1.0/std::sqrt(1-v.squaredNorm()/(c*c));
    const Vec3 Ep=(f.electric+cross(v,f.magnetic))*gu-v*(gu*gu/(gu+1)*dot(v,f.electric)/(c*c));
    const Vec3 S0=(first?s.firstProperDipole:s.secondProperDipole)*(1.0/(first?g1:g2));
    const double q=first?firstCharge:secondCharge, m=first?firstMass:secondMass;
    const Vec3 pp=cross(Ep,S0)*(kappa*q/(m*c*c));
    const double vn=v.norm(); const Vec3 vh=vn>0?v*(1.0/vn):Vec3{};
    return H4{gu*dot(v,pp)/c, pp+vh*((gu-1)*dot(pp,vh))}; };
  const auto thomasForce=[&](auto build,double t,int i)->Vec3{
    if(kappa==0.0) return {};
    const double hh=1e-4/w;
    const H4 a=thomasAt(build,t+hh,i), b=thomasAt(build,t-hh,i);
    const Vec3 dps=(a.p-b.p)*(1.0/(2*hh)); const double dp0=(a.p0-b.p0)/(2*hh);
    const State s=build(t); const Vec3 v=i==1?s.firstVelocity:s.secondVelocity;
    const double gu=1.0/std::sqrt(1-v.squaredNorm()/(c*c));
    const double K0=-gu*dp0; const Vec3 Ks=dps*(-gu);
    const double Ku=K0*gu*c-dot(Ks,v*gu);
    return (Ks-v*(gu*Ku/(c*c)))*(1.0/gu); };
  int which=1;   // the particle whose event is compared (audit 568, uzupelnienie 1)
  const auto forcesOf=[&](auto build,double t0,Vec3& full,Vec3& charge){
    StateHistory h; for(int i=0;i<nodes;++i) h.push_back(build(t0-span*lt*(1.0-double(i)/nodes)));
    const State s=build(t0);
    gDipoleForceEnabled=true; const MutualForces fo=retardedExternalForces(s,h); full=which==1?fo.first:fo.second;
    full=full+hiddenForce(build,t0,which,covariantMode)+thomasForce(build,t0,which);
    gDipoleForceEnabled=false; const MutualForces fc=retardedExternalForces(s,h); charge=which==1?fc.first:fc.second;
    gDipoleForceEnabled=true; };
  const auto transform=[&](const Vec3& Fp,const Vec3& up){
    const double gu=1.0/std::sqrt(1-up.squaredNorm()/(c*c));
    const double K0=gu*dot(Fp,up)/c; const Vec3 K=Fp*gu;
    const double K0l=G*(K0+(V/c)*K.y); (void)K0l;
    const Vec3 Kl{K.x,G*(K.y+(V/c)*K0),K.z};
    const double gul=G*gu*(1+V*up.y/(c*c));
    return Kl*(1.0/gul); };
  std::printf("[%s hidden momentum, Thomas kappa %.2f] ",covariantMode?"covariant":"lab-form",kappa);
  std::printf("%s: d %.1f r*, V/c %.2f, beta_orbit(1) %.4e, nodes %d span %.0f d/c\n",argv[1],d/comptonBarrierRadius,V/c,R1*w/c,nodes,span);
  for(which=1;which<=2;++which){
  std::printf(" particle %d\n",which);
  double sF=0,sD=0,sC=0, nF=0,nD=0,nC=0;
  for(int p=0;p<8;++p){
    const double tp0=(2*pi*p/8)/w;
    Vec3 fFr,fCr; forcesOf(restState,tp0,fFr,fCr);
    const double t0=G*(tp0+V*restX(which,tp0).y/(c*c));
    Vec3 fFl,fCl; forcesOf(labState,t0,fFl,fCl);
    const Vec3 up=restU(which,tp0);
    const Vec3 tF=transform(fFr,up), tC=transform(fCr,up), tD=transform(fFr-fCr,up);
    const Vec3 dl=fFl-fCl;
    sF+=(fFl-tF).squaredNorm(); nF+=fFl.squaredNorm();
    sC+=(fCl-tC).squaredNorm(); nC+=fCl.squaredNorm();
    sD+=(dl-tD).squaredNorm(); nD+=dl.squaredNorm();
    std::printf("  phase %d  |F_dip lab| %.4e  lab %+.4e %+.4e %+.4e  transf %+.4e %+.4e %+.4e\n",p,dl.norm(),dl.x,dl.y,dl.z,tD.x,tD.y,tD.z);
  }
  std::printf("  RMS non-covariance: charges %.3e  whole %.3e  dipole sector %.4e (of |F_dip|)\n",
    std::sqrt(sC/nC),std::sqrt(sF/nF),std::sqrt(sD/nD));
  }
}
