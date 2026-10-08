#include "modules/crem_collapse.hpp"
#include <atomic>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <random>
#include <thread>
// Audit 376: does a 50 uT (Earth-scale) field quantize free spins at the engine's final state (n = 1, L = hbar/2)?
// Rates from the engine (orbitAveragedBmtAngularVelocities + the s-state isotropy rule's own decomposition):
//   omega_i = f_i(|L|) Lhat + (tr A_i/3) m_j + G_i B,
// then dm_i/dt = omega_i x m_i, dL/dt = -sum_i (internal omega_i x m_i)/gamma_i (J changes only through the field).
// Mode "check": reduced equations vs advanceCoupledSecularSpinOrbit over 10 ps (B = 0, 1 T).
// Mode "run": 200 free draws per channel, isotropic field direction, B = 0 / 50 uT / 1 T, to 1 us.

static double mred, aBohr, g1, g2;
struct Coeff { double fL[3]; double Ls[3]; double tr1, tr2; double G1[3][3], G2[3][3]; };
static Coeff C;

static Vec3 axisOf(int k){ return k==0?Vec3{1,0,0}:(k==1?Vec3{0,1,0}:Vec3{0,0,1}); }

static void extract(){
  const double tiny=1e-12;
  const Vec3 m1=Vec3{1,0,0}*firstMagneticMoment, m2=Vec3{1,0,0}*secondMagneticMoment;
  // spin-orbit part f_i(|L|) along Lhat at three |L|
  const double Ls[3]={0.49,0.50,0.51};
  for(int k=0;k<3;++k){
    const Vec3 L{0,0,Ls[k]*hbar};
    const auto o=orbitAveragedBmtAngularVelocities(aBohr,L,m1*tiny,m2*tiny,mred,0.0,Vec3{1,0,0});
    C.Ls[k]=Ls[k]; C.fL[k]=o.first.z;
    if(k==1) std::printf("orbital part at L = hbar/2: first (%.4e %.4e %.4e) second (%.4e %.4e %.4e) rad/s\n",
        o.first.x,o.first.y,o.first.z,o.second.x,o.second.y,o.second.z);
  }
  // f for the second particle stored by symmetry check below
  const Vec3 L{0,0,0.5*hbar};
  const auto orb=orbitAveragedBmtAngularVelocities(aBohr,L,m1*tiny,m2*tiny,mred,0.0,Vec3{1,0,0});
  double t1=0,t2=0;
  for(int k=0;k<3;++k){
    const Vec3 ax=axisOf(k);
    const auto p=orbitAveragedBmtAngularVelocities(aBohr,L,ax*firstMagneticMoment,ax*secondMagneticMoment,mred,0.0,Vec3{1,0,0});
    t1+=dot(ax,p.first-orb.first)/secondMagneticMoment;
    t2+=dot(ax,p.second-orb.second)/firstMagneticMoment;
  }
  C.tr1=t1/3; C.tr2=t2/3;
  // external: G_i from B = 1 T along each axis (the engine's own firstExternal/secondExternal)
  const Vec3 saved=gExternalMagneticField;
  for(int k=0;k<3;++k){
    gExternalMagneticField=axisOf(k)*1.0;
    const auto e=orbitAveragedBmtAngularVelocities(aBohr,L,m1*tiny,m2*tiny,mred,0.0,Vec3{1,0,0});
    for(int r=0;r<3;++r){ C.G1[r][k]=(r==0?e.firstExternal.x:r==1?e.firstExternal.y:e.firstExternal.z);
                          C.G2[r][k]=(r==0?e.secondExternal.x:r==1?e.secondExternal.y:e.secondExternal.z); }
  }
  gExternalMagneticField=saved;
  std::printf("isotropic partner coefficient: tr A1/3 %.6e, tr A2/3 %.6e (rad/s per A m^2)\n",C.tr1,C.tr2);
  std::printf("  -> partner precession rates |tr/3| m: first %.4e, second %.4e rad/s\n",std::abs(C.tr1)*secondMagneticMoment,std::abs(C.tr2)*firstMagneticMoment);
  std::printf("G1 diag %.6e %.6e %.6e offdiag max %.2e; G2 diag %.6e %.6e %.6e (rad/s/T)\n",C.G1[0][0],C.G1[1][1],C.G1[2][2],
      std::max({std::abs(C.G1[0][1]),std::abs(C.G1[0][2]),std::abs(C.G1[1][2]),std::abs(C.G1[1][0])}),C.G2[0][0],C.G2[1][1],C.G2[2][2]);
  std::printf("gamma1 %.6e gamma2 %.6e; eB/m %.6e rad/s/T\n",g1,g2,elementaryCharge/firstMass);
}

static double fOf(double Lh,const double* f){ // linear in |L|/hbar through the three points
  const double s=(f[2]-f[0])/(C.Ls[2]-C.Ls[0]); return f[1]+s*(Lh-C.Ls[1]); }
static double f2tab[3];
static double gRotation=1.0;

static Vec3 mat(const double G[3][3],const Vec3& B){
  return Vec3{G[0][0]*B.x+G[0][1]*B.y+G[0][2]*B.z,G[1][0]*B.x+G[1][1]*B.y+G[1][2]*B.z,G[2][0]*B.x+G[2][1]*B.y+G[2][2]*B.z}; }

struct Y { Vec3 L,m1,m2; };
static Y deriv(const Y& y,const Vec3& B){
  const double Ln=y.L.norm(), Lh=Ln/hbar; const Vec3 Lhat=Ln>0?y.L*(1.0/Ln):Vec3{0,0,1};
  const Vec3 w1i=Lhat*fOf(Lh,C.fL)+y.m2*C.tr1, w2i=Lhat*fOf(Lh,f2tab)+y.m1*C.tr2;
  const Vec3 w1=w1i+mat(C.G1,B), w2=w2i+mat(C.G2,B);
  Y d; d.m1=cross(w1,y.m1); d.m2=cross(w2,y.m2);
  d.L=-(cross(w1i,y.m1)*(1.0/g1)+cross(w2i,y.m2)*(1.0/g2));
  return d;
}
static Y axpy(const Y& y,const Y& d,double h){ return Y{y.L+d.L*h,y.m1+d.m1*h,y.m2+d.m2*h}; }
static void rk4(Y& y,const Vec3& B,double h){
  const Y k1=deriv(y,B),k2=deriv(axpy(y,k1,h/2),B),k3=deriv(axpy(y,k2,h/2),B),k4=deriv(axpy(y,k3,h),B);
  y.L+= (k1.L+k2.L*2+k3.L*2+k4.L)*(h/6); y.m1+=(k1.m1+k2.m1*2+k3.m1*2+k4.m1)*(h/6); y.m2+=(k1.m2+k2.m2*2+k3.m2*2+k4.m2)*(h/6);
  y.m1=y.m1*(firstMagneticMoment/y.m1.norm()); y.m2=y.m2*(secondMagneticMoment/y.m2.norm());
}
static Vec3 rot(const Vec3& v,const Vec3& axis,double ang){ // Rodrigues, axis unit
  const double c=std::cos(ang),sn=std::sin(ang); return v*c+cross(axis,v)*sn+axis*(dot(axis,v)*(1-c)); }
// rest of the flow: isotropic partner (Heisenberg: no torque on L) + external field
static Y derivRest(const Y& y,const Vec3& B){ Y d; d.L=Vec3{};
  d.m1=cross(y.m2*C.tr1+mat(C.G1,B),y.m1); d.m2=cross(y.m1*C.tr2+mat(C.G2,B),y.m2); return d; }
static void rk4Rest(Y& y,const Vec3& B,double h){
  const Y k1=derivRest(y,B),k2=derivRest(axpy(y,k1,h/2),B),k3=derivRest(axpy(y,k2,h/2),B),k4=derivRest(axpy(y,k3,h),B);
  y.m1+=(k1.m1+k2.m1*2+k3.m1*2+k4.m1)*(h/6); y.m2+=(k1.m2+k2.m2*2+k3.m2*2+k4.m2)*(h/6);
  y.m1=y.m1*(firstMagneticMoment/y.m1.norm()); y.m2=y.m2*(secondMagneticMoment/y.m2.norm()); }
// exact L.S flow (kappa L.S, kappa = f/|L|): S_i(t) = R_J(kappa|J|t) R_S0(-kappa|S|t) S_i(0), L(t) = R_J L(0)
static double kappaOf(const Y& y){ return fOf(y.L.norm()/hbar,C.fL)/y.L.norm(); }
static void exactLS(Y& y,double h){
  const double k=kappaOf(y);
  Vec3 S1=y.m1*(1.0/g1), S2=y.m2*(1.0/g2); const Vec3 S=S1+S2, J=y.L+S;
  const double Sn=S.norm(), Jn=J.norm();
  if(Sn>0){ const Vec3 a=S*(1.0/Sn); S1=rot(S1,a,-k*Sn*h); S2=rot(S2,a,-k*Sn*h); }
  if(Jn>0){ const Vec3 a=J*(1.0/Jn); S1=rot(S1,a,k*Jn*h); S2=rot(S2,a,k*Jn*h); y.L=rot(y.L,a,k*Jn*h); }
  y.m1=S1*g1; y.m2=S2*g2; }
static double splitStepSize(const Y& y,const Vec3& B,double rotation){
  const double k=std::abs(kappaOf(y)); const Vec3 J=y.L+y.m1*(1.0/g1)+y.m2*(1.0/g2);
  const double rest=std::max((y.m2*C.tr1+mat(C.G1,B)).norm(),(y.m1*C.tr2+mat(C.G2,B)).norm());
  return std::min(rotation/(k*J.norm()),0.02/std::max(rest,1.0)); }
static void splitStep(Y& y,const Vec3& B,double h){ rk4Rest(y,B,h/2); exactLS(y,h); rk4Rest(y,B,h/2); }
static double maxRate(const Y& y,const Vec3& B){
  const Vec3 Lhat=y.L*(1.0/y.L.norm());
  const double a=(Lhat*fOf(y.L.norm()/hbar,C.fL)+y.m2*C.tr1+mat(C.G1,B)).norm();
  const double b=(Lhat*fOf(y.L.norm()/hbar,f2tab)+y.m1*C.tr2+mat(C.G2,B)).norm();
  return std::max(a,b);
}
static Vec3 iso(std::mt19937_64& r){ std::uniform_real_distribution<double> u(0,1);
  const double c=2*u(r)-1,p=2*pi*u(r),s=std::sqrt(1-c*c); return Vec3{s*std::cos(p),s*std::sin(p),c}; }

int main(int argc,char** argv){
  mred=firstMass*secondMass/(firstMass+secondMass); aBohr=pairBohrRadius(activePair);
  g1=firstGyromagneticRatioOf(); g2=secondGyromagneticRatioOf();
  extract();
  { const double tiny=1e-12; const Vec3 m1=Vec3{1,0,0}*firstMagneticMoment*tiny, m2=Vec3{1,0,0}*secondMagneticMoment*tiny;
    for(int k=0;k<3;++k){ const auto o=orbitAveragedBmtAngularVelocities(aBohr,Vec3{0,0,C.Ls[k]*hbar},m1,m2,mred,0.0,Vec3{1,0,0}); f2tab[k]=o.second.z; } }
  const double eps=orePowellSuppression(), sv=4*pi*classicalElectronRadius*classicalElectronRadius*c;
  const char* mode=argc>1?argv[1]:"check";
  if(!std::strcmp(mode,"check")){
    for(double Bmag: {0.0,1.0}){
      for(int trial=0;trial<3;++trial){
        std::mt19937_64 r(100+trial);
        const Vec3 m1=iso(r)*firstMagneticMoment, m2=iso(r)*secondMagneticMoment, Bd=iso(r);
        gExternalMagneticField=Bd*Bmag;
        const SecularSpinOrbitState st{Vec3{0,0,0.5*hbar},m1,m2,0.0,Vec3{1,0,0}};
        const double T=1e-11;
        const double sub=argc>2?std::atof(argv[2]):0.05;
        const auto e=advanceCoupledSecularSpinOrbit(st,aBohr,mred,T,sub,1<<24);
        Y y{Vec3{0,0,0.5*hbar},m1,m2}; const Vec3 B=Bd*Bmag;
        const int n=(int)std::ceil(T*maxRate(y,B)/0.002); const double h=T/n;
        for(int i=0;i<n;++i) rk4(y,B,h);
        const auto cosOf=[](const Vec3&a,const Vec3&b){return dot(a,b)/(a.norm()*b.norm());};
        const auto ang=[&](const Vec3&a,const Vec3&b){return std::acos(std::clamp(cosOf(a,b),-1.0,1.0));};
        std::printf("check B=%.0f T trial %d: engine completed %d substeps %d; cos0 %.6f engine %.6f reduced %.6f; angle m1 %.2e m2 %.2e rad; |L| engine %.6f reduced %.6f hbar\n",
          Bmag,trial,e.completed,e.substeps,cosOf(m1,m2),cosOf(e.state.firstDipole,e.state.secondDipole),cosOf(y.m1,y.m2),
          ang(e.state.firstDipole,y.m1),ang(e.state.secondDipole,y.m2),e.state.orbitalAngularMomentum.norm()/hbar,y.L.norm()/hbar);
      }
    }
    gExternalMagneticField=Vec3{};
    return 0;
  }
  if(!std::strcmp(mode,"splitcheck")){
    const double T=2e-9;
    for(double Bmag: {50e-6,1.0}) for(double R: {0.3,1.0}) for(int trial=0;trial<3;++trial){
      std::mt19937_64 r(200+trial); const Vec3 m1=iso(r)*firstMagneticMoment, m2=iso(r)*secondMagneticMoment, B=iso(r)*Bmag;
      Y a{Vec3{0,0,0.5*hbar},m1,m2}, b=a; double t=0;
      while(t<T){ const double h=std::min(0.005/maxRate(a,B),T-t); rk4(a,B,h); t+=h; }
      t=0; while(t<T){ const double h=std::min(splitStepSize(b,B,R),T-t); splitStep(b,B,h); t+=h; }
      std::printf("splitcheck B %.0e rot %.1f trial %d: w0 %.9f plainRK4 %.9f split %.9f  |L| %.6f %.6f\n",Bmag,R,trial,
        contactTwoPhotonWeight(m1,m2),contactTwoPhotonWeight(a.m1,a.m2),contactTwoPhotonWeight(b.m1,b.m2),a.L.norm()/hbar,b.L.norm()/hbar);
    }
    std::printf("f1 %.9e f2 %.9e (L = hbar/2); tr1/g1 %.9e tr2/g2 %.9e\n",C.fL[1],f2tab[1],C.tr1/g1,C.tr2/g2);
    return 0;
  }
  if(!std::strcmp(mode,"conv1T")){
    for(double T: {1e-10,5e-10,2e-9}) for(int trial=0;trial<3;++trial){
      std::mt19937_64 r(200+trial); const Vec3 m1=iso(r)*firstMagneticMoment, m2=iso(r)*secondMagneticMoment, B=iso(r)*1.0;
      std::printf("T %.0e trial %d:",T,trial);
      for(double q: {0.005,0.0025,0.00125}){ Y a{Vec3{0,0,0.5*hbar},m1,m2}; double t=0;
        while(t<T){ const double h=std::min(q/maxRate(a,B),T-t); rk4(a,B,h); t+=h; } std::printf(" rk4(%g) %.9f",q,contactTwoPhotonWeight(a.m1,a.m2)); }
      for(double R: {0.1,0.03}){ Y b{Vec3{0,0,0.5*hbar},m1,m2}; double t=0;
        while(t<T){ const double h=std::min(splitStepSize(b,B,R),T-t); splitStep(b,B,h); t+=h; } std::printf(" split(%g) %.9f",R,contactTwoPhotonWeight(b.m1,b.m2)); }
      std::printf("\n"); }
    return 0;
  }
  if(!std::strcmp(mode,"conv50")){
    const double T=1e-7;
    for(int trial=0;trial<4;++trial){
      std::mt19937_64 r(300+trial); const Vec3 m1=iso(r)*firstMagneticMoment, m2=iso(r)*secondMagneticMoment, B=iso(r)*50e-6;
      std::printf("T 1e-7 trial %d: w0 %.9f",trial,contactTwoPhotonWeight(m1,m2));
      for(double q: {0.02,0.01}){ Y a{Vec3{0,0,0.5*hbar},m1,m2}; double t=0;
        while(t<T){ const double h=std::min(q/maxRate(a,B),T-t); rk4(a,B,h); t+=h; } std::printf(" rk4(%g) %.9f",q,contactTwoPhotonWeight(a.m1,a.m2)); }
      for(double R: {1.0,0.3}){ Y b{Vec3{0,0,0.5*hbar},m1,m2}; double t=0;
        while(t<T){ const double h=std::min(splitStepSize(b,B,R),T-t); splitStep(b,B,h); t+=h; } std::printf(" split(%g) %.9f",R,contactTwoPhotonWeight(b.m1,b.m2)); }
      std::printf("\n"); std::fflush(stdout); }
    return 0;
  }
  if(!std::strcmp(mode,"proj")){
    // Audit 376f: projections of the net moment (m1+m2) and of J = L + S1 + S2 on the field direction, start vs end,
    // free spins, 50 uT, 200 draws per channel, to T; histograms of the end/start values in 10 bins.
    const int N=argc>2?std::atoi(argv[2]):200; const double T=argc>3?std::atof(argv[3]):1e-6; const double Bmag=50e-6;
    for(int ch=0;ch<2;++ch){
      int h0[10]={0},h1[10]={0},hJ0[10]={0},hJ1[10]={0}; double maxdm=0,maxdJ=0,netmean=0;
      for(int i=0;i<N;++i){
        std::mt19937_64 r(42+1000*ch+i);
        const Vec3 m1=iso(r)*firstMagneticMoment; Vec3 m2;
        for(int a=0;a<10000;++a){ m2=iso(r)*secondMagneticMoment; const double cc=dot(m1,m2)/(firstMagneticMoment*secondMagneticMoment); if((cc>=0.5)==(ch==0)) break; }
        std::mt19937_64 rb(777+i); const Vec3 Bh=iso(rb), B=Bh*Bmag;
        Y y{Vec3{0,0,0.5*hbar},m1,m2};
        const auto proj=[&](const Y& z){ return dot(z.m1+z.m2,Bh)/(2*firstMagneticMoment); };      // in units of 2 mu
        const auto Jp=[&](const Y& z){ return dot(z.L+z.m1*(1.0/g1)+z.m2*(1.0/g2),Bh)/hbar; };
        const double p0=proj(y), j0=Jp(y); netmean+=(y.m1+y.m2).norm()/(2*firstMagneticMoment);
        double t=0; while(t<T){ const double h=std::min(splitStepSize(y,B,1.0),T-t); splitStep(y,B,h); t+=h; }
        const double p1=proj(y), j1=Jp(y);
        maxdm=std::max(maxdm,std::abs(p1-p0)); maxdJ=std::max(maxdJ,std::abs(j1-j0));
        const auto bin=[](double x,double lo,double hi){ return std::clamp((int)((x-lo)/(hi-lo)*10),0,9); };
        ++h0[bin(p0,-1,1)]; ++h1[bin(p1,-1,1)]; ++hJ0[bin(j0,-2,2)]; ++hJ1[bin(j1,-2,2)];
      }
      std::printf("%s (free spins, 50 uT, T = %.0e s, N = %d): <|m1+m2|>/2mu %.4f; max|d proj m| %.3e (2mu), max|d J_B| %.3e hbar\n",ch==0?"p-Ps":"o-Ps",T,N,netmean/N,maxdm,maxdJ);
      std::printf("  (m1+m2).B/2mu bins [-1,1]: start"); for(int b=0;b<10;++b) std::printf(" %d",h0[b]); std::printf("  end"); for(int b=0;b<10;++b) std::printf(" %d",h1[b]); std::printf("\n");
      std::printf("  J.B/hbar bins [-2,2]:      start"); for(int b=0;b<10;++b) std::printf(" %d",hJ0[b]); std::printf("  end"); for(int b=0;b<10;++b) std::printf(" %d",hJ1[b]); std::printf("\n");
    }
    return 0;
  }
  // run
  if(argc>4) gRotation=std::atof(argv[4]);
  const int N=argc>2?std::atoi(argv[2]):200; const double T=argc>3?std::atof(argv[3]):1e-6;
  const double fields[3]={0.0,50e-6,1.0};
  const double marks[5]={0.0,1e-9,1e-8,1e-7,1e-6};
  struct Out { double w[5]; double maxdw, tau; };
  std::vector<Out> out(2*3*N);
  std::atomic<int> next{0}; std::mutex mu;
  auto work=[&](){ while(true){ const int job=next.fetch_add(1); if(job>=2*3*N) break;
      const int ch=job/(3*N), fi=(job/N)%3, i=job%N;
      std::mt19937_64 r(42+1000*ch+i);
      const Vec3 m1=iso(r)*firstMagneticMoment; Vec3 m2;
      for(int a=0;a<10000;++a){ m2=iso(r)*secondMagneticMoment; const double cc=dot(m1,m2)/(firstMagneticMoment*secondMagneticMoment); if((cc>=0.5)==(ch==0)) break; }
      std::mt19937_64 rb(777+i); const Vec3 B=iso(rb)*fields[fi];
      Y y{Vec3{0,0,0.5*hbar},m1,m2};
      Out o{}; const double w0=contactTwoPhotonWeight(m1,m2); o.w[0]=w0;
      double t=0,S=1,tauAcc=0; int mk=1;
      const auto gam=[&](const Y& yy){ const double w=contactTwoPhotonWeight(yy.m1,yy.m2);
        const double n=1.0/(2*pi*aBohr*aBohr*aBohr*(yy.L.norm()/hbar)); return sv*n*(w+(1-w)*eps); };
      // 1 T control: chaotic (conv1T.txt), plain RK4 (step 0.005 rad) to 10 ns, statistics only.
      const bool control=fi==2; const double Tend=control?std::min(T,1e-8):T;
      while(t<Tend){
        const double h=control?std::min(0.005/maxRate(y,B),Tend-t):std::min(splitStepSize(y,B,gRotation),Tend-t);
        const double g0=gam(y); if(control) rk4(y,B,h); else splitStep(y,B,h); t+=h; const double g=gam(y);
        const double rate=0.5*(g0+g), dec=std::exp(-rate*h); tauAcc+=S*(1-dec)/rate; S*=dec;
        const double w=contactTwoPhotonWeight(y.m1,y.m2); o.maxdw=std::max(o.maxdw,std::abs(w-w0));
        while(mk<5&&t>=marks[mk]*(1-1e-12)){ o.w[mk]=w; ++mk; }
      }
      for(;mk<5;++mk) o.w[mk]=contactTwoPhotonWeight(y.m1,y.m2);
      o.tau=tauAcc+S/gam(y);
      out[job]=o;
      { std::lock_guard<std::mutex> l(mu); std::printf("TRAJ ch %d B %.3e i %d w %.9f %.9f %.9f %.9f %.9f maxdw %.6e tau %.9e\n",ch,fields[fi],i,o.w[0],o.w[1],o.w[2],o.w[3],o.w[4],o.maxdw,o.tau); std::fflush(stdout); }
      if(i%50==49){ std::lock_guard<std::mutex> l(mu); std::fprintf(stderr,"ch %d B %.0e i %d\n",ch,fields[fi],i+1); }
    } };
  std::vector<std::thread> th; for(unsigned k=0;k<std::max(1u,std::thread::hardware_concurrency());++k) th.emplace_back(work);
  for(auto& x:th) x.join();
  std::printf("\nrun: N = %d per channel, T = %.1e s, field direction isotropic (same per draw index for every B)\n",N,T);
  for(int ch=0;ch<2;++ch) for(int fi=0;fi<3;++fi){
    std::vector<double> md; double tau=0,ends0=0,endsT=0,wm[5]={0,0,0,0,0};
    for(int i=0;i<N;++i){ const Out& o=out[ch*3*N+fi*N+i]; md.push_back(o.maxdw); tau+=o.tau;
      for(int k=0;k<5;++k) wm[k]+=o.w[k];
      ends0+=(o.w[0]<0.01||o.w[0]>0.99); endsT+=(o.w[4]<0.01||o.w[4]>0.99); }
    std::sort(md.begin(),md.end());
    std::printf("%s B = %-7g T: <tau> %.6g ns; <w> at 0/1ns/10ns/100ns/1us %.6f %.6f %.6f %.6f %.6f; max|dw| median %.3e p90 %.3e max %.3e; ends (w<0.01|w>0.99) start %.0f end %.0f /%d\n",
      ch==0?"p-Ps":"o-Ps",fields[fi],tau/N*1e9,wm[0]/N,wm[1]/N,wm[2]/N,wm[3]/N,wm[4]/N,md[N/2],md[(9*N)/10],md[N-1],ends0,endsT,N);
  }
  // paired lifetime ratios
  for(int ch=0;ch<2;++ch){ double r50=0,r1=0; for(int i=0;i<N;++i){ const double t0=out[ch*3*N+i].tau;
      r50+=out[ch*3*N+N+i].tau/t0; r1+=out[ch*3*N+2*N+i].tau/t0; }
    std::printf("%s paired <tau(B)/tau(0)>: 50 uT %.8f, 1 T %.6f\n",ch==0?"p-Ps":"o-Ps",r50/N,r1/N); }
}
