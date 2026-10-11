#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <string>
#include <vector>
// Audit 407 (from h0406): argv: dump [pair or -] [scan].
// Audit 406: M1 power and moment derivatives at a dumped state (CREM_DUMP_FAIL + CREM_DUMP_BELOW_R).
// Compares the history-differenced m', m'', m''' (historicalDipoleKinematics, what production uses for the M1 power)
// with the analytic ones from the Thomas-BMT precession, Omega = -(q/m) B_eff:
//   full:      m' = W x m, m'' = W' x m + W x (W x m), m''' = W'' x m + 2 W' x (W x m) + W x (W' x m) + W x (W x (W x m))
//   no-Wdot:   m'' = W x (W x m), m''' = W x (W x (W x m))              (audit 404's torque)
// W' and W'' from W evaluated on the state extrapolated by +-delta, +-2 delta along its own motion.
// argv: dump [pair or "-"].
static State start; static StateHistory hist; static double dtStep=0.0;
static void load(const char* path){
  std::FILE* f=std::fopen(path,"rb"); if(!f){std::perror(path); std::exit(1);}
  std::uint64_t nodes=0;
  if(std::fread(&dtStep,sizeof dtStep,1,f)!=1||std::fread(&start,sizeof(State),1,f)!=1||std::fread(&nodes,sizeof nodes,1,f)!=1) std::exit(1);
  for(std::uint64_t i=0;i<nodes;++i){State s; if(std::fread(&s,sizeof s,1,f)!=1) std::exit(1); hist.push_back(s);}
  std::fclose(f);
}
struct Omegas { Vec3 first, second; };
static Omegas omegaAt(const State& s){
  const LocalElectromagneticFields f=localRelativisticFields(s,hist);
  return {thomasBmtEffectiveField(s.firstVelocity,f.atFirst,firstGFactor)*(-firstCharge/firstMass),
          thomasBmtEffectiveField(s.secondVelocity,f.atSecond,secondGFactor)*(-secondCharge/secondMass)};
}
static State shifted(double d,const Omegas& w){
  State s=start; s.time+=d;
  s.firstPosition+=s.firstVelocity*d+s.firstAcceleration*(0.5*d*d); s.secondPosition+=s.secondVelocity*d+s.secondAcceleration*(0.5*d*d);
  s.firstVelocity+=s.firstAcceleration*d; s.secondVelocity+=s.secondAcceleration*d;
  s.firstDipole+=cross(w.first,s.firstDipole)*d; s.secondDipole+=cross(w.second,s.secondDipole)*d;
  return s;
}
int main(int argc,char**argv){
  if(argc<2){std::fprintf(stderr,"usage: h0406 dump [pair]\n");return 1;}
  if(argc>2&&std::string(argv[2])!="-") applyPairFromOption(argv[2]);
  load(argv[1]);
  const double r=separation(start), v=(start.firstVelocity-start.secondVelocity).norm();
  const Omegas w0=omegaAt(start);
  const double wmax=std::max(w0.first.norm(),w0.second.norm());
  const double delta=1.0e-3*std::min(1.0/std::max(wmax,1e-300),r/std::max(v,1e-300));
  const Omegas wp=omegaAt(shifted(delta,w0)), wm=omegaAt(shifted(-delta,w0));
  const Omegas wp2=omegaAt(shifted(2*delta,w0)), wm2=omegaAt(shifted(-2*delta,w0));
  std::printf("state: t=%.9e r=%.6e m v_rel/c=%.4e omega_orb=v/r=%.4e |W1|=%.4e |W2|=%.4e history=%zu delta=%.3e\n",start.time,r,v/c,v/r,
              w0.first.norm(),w0.second.norm(),hist.size(),delta);
  if(argc>3&&std::string(argv[3])=="scan"){
    // Omega of the FIRST particle along the state's own extrapolated motion, split into the Lienard-Wiechert
    // and the magnetic-dipole parts of the partner's field; second differences at the base step and the
    // history nodes the retarded time / the dipole-derivative stencil (t_ret +- h, +- 2h) cross.
    const double hDer=historyDerivativeStep(hist,2.0);
    const auto parts=[&](const State& s){
      const ElectromagneticField lw=lienardWiechertField(s.firstPosition,s.time,hist,s,false,secondCharge);
      const ElectromagneticField dp=retardedMagneticDipoleField(s.firstPosition,s.time,hist,s,false);
      const double k=-firstCharge/firstMass;
      return std::pair<Vec3,Vec3>{thomasBmtEffectiveField(s.firstVelocity,lw,firstGFactor)*k,
                                  thomasBmtEffectiveField(s.firstVelocity,dp,firstGFactor)*k}; };
    const auto tret=[&](const State& s){ double tr=s.time-separation(s)/c;
      for(int it=0;it<60;++it){ const ChargeKinematics q=historicalCharge(hist,s,false,tr); tr=s.time-(s.firstPosition-q.position).norm()/c; }
      return tr; };
    const auto nodeIndex=[&](double t){ std::size_t k=0; while(k<hist.size()&&hist[k].time<t) ++k; return k; };
    const int N=128; std::vector<std::pair<Vec3,Vec3>> om(2*N+1); std::vector<double> tr(2*N+1);
    for(int i=-N;i<=N;++i){ const State s=shifted(i*delta,w0); om[i+N]=parts(s); tr[i+N]=tret(s); }
    std::printf("scan: base delta=%.3e s, history derivative step h=%.3e s, node spacing ~%.3e s\n",delta,hDer,
      hist.size()>2?(hist.back().time-hist.front().time)/(hist.size()-1):0.0);
    for(int i=-N+1;i<N;++i){
      const Vec3 d2lw=(om[i+N+1].first-om[i+N].first*2.0+om[i+N-1].first)*(1.0/(delta*delta));
      const Vec3 d2dp=(om[i+N+1].second-om[i+N].second*2.0+om[i+N-1].second)*(1.0/(delta*delta));
      std::string cross;
      for(double off: {0.0,-hDer,hDer,-2*hDer,2*hDer}){ const std::size_t a=nodeIndex(tr[i+N-1]+off), b=nodeIndex(tr[i+N+1]+off);
        if(a!=b){ char buf[64]; std::snprintf(buf,sizeof buf," node(t_ret%+.0fh)",off/hDer); cross+=buf; } }
      std::printf("d/delta=%+4d |W_LW''|=%.4e |W_dip''|=%.4e%s\n",i,d2lw.norm(),d2dp.norm(),cross.c_str()); }
    // least-squares quadratic fits of W (total) over +-K delta: W'' = 2 a2
    for(int K: {2,4,8,16,32,64,128}){ double s0=0,s2=0,s4=0; Vec3 y0,y2;
      for(int i=-K;i<=K;++i){ const double x=i*delta; const Vec3 w=om[i+N].first+om[i+N].second;
        s0+=1; s2+=x*x; s4+=x*x*x*x; y0+=w; y2+=w*(x*x); }
      const double det=s0*s4-s2*s2; const Vec3 a2=(y2*s0-y0*s2)*(1.0/det);
      Vec3 y1; double s2b=0; for(int i=-K;i<=K;++i){ const double x=i*delta; y1+=(om[i+N].first+om[i+N].second)*x; s2b+=x*x; }
      std::printf("LSQ fit +-%3d delta: |W'|=%.6e |W''|=%.6e\n",K,(y1*(1.0/s2b)).norm(),(a2*2.0).norm()); }
    { // production analyticDipoleKinematics vs harness LSQ (+-32 base steps), first particle
      const int K=32; double s0=0,s2=0,s4=0; Vec3 y0,y1,y2;
      for(int i=-K;i<=K;++i){ const double x=i*delta; const Vec3 w=om[i+N].first+om[i+N].second;
        s0+=1; s2+=x*x; s4+=x*x*x*x; y0+=w; y1+=w*x; y2+=w*(x*x); }
      const Vec3 Wd=y1*(1.0/s2), Wdd=(y2*s0-y0*s2)*(2.0/(s0*s4-s2*s2)), W=w0.first, m=start.firstDipole;
      const Vec3 m3=cross(Wdd,m)+cross(Wd,cross(W,m))*2.0+cross(W,cross(Wd,m))+cross(W,cross(W,cross(W,m)));
      const Vec3 m2=cross(Wd,m)+cross(W,cross(W,m));
      const AnalyticDipoleKinematicsPair k=analyticDipoleKinematics(start,hist);
      std::printf("first particle: m'' production %.6e harness-LSQ %.6e; m''' production %.6e harness-LSQ %.6e\n",
        k.first.secondDerivative.norm(),m2.norm(),k.first.thirdDerivative.norm(),m3.norm()); }
    return 0;
  }
  // delta convergence of W' and W'' (first particle)
  for(double f: {0.25,1.0,4.0,16.0,64.0,256.0}){ const double d=delta*f;
    const Omegas p=omegaAt(shifted(d,w0)), q=omegaAt(shifted(-d,w0));
    std::printf("  delta x %-5g |W1'|=%.6e |W1''|=%.6e\n",f,((p.first-q.first)*(1.0/(2*d))).norm(),((p.first-w0.first*2.0+q.first)*(1.0/(d*d))).norm()); }
  constexpr double coef=mu0/(6.0*pi*c*c*c);
  Vec3 sumHist2, sumFull2, sumNo2, sumHist3, sumFull3, sumNo3;
  for(int k=0;k<2;++k){ const bool first=k==0;
    const Vec3 m=first?start.firstDipole:start.secondDipole;
    const Vec3 W=first?w0.first:w0.second;
    const Vec3 Wp=first?wp.first:wp.second, Wm=first?wm.first:wm.second, Wp2=first?wp2.first:wp2.second, Wm2=first?wm2.first:wm2.second;
    const Vec3 Wd=(Wp-Wm)*(1.0/(2*delta));
    const Vec3 Wdd=(Wp-W*2.0+Wm)*(1.0/(delta*delta));
    const Vec3 Wdd4=(Wp2-W*2.0+Wm2)*(1.0/(4*delta*delta));
    const RetardedDipoleKinematics h=historicalDipoleKinematics(hist,start,first,start.time);
    const Vec3 a1=cross(W,m);
    const Vec3 a2=cross(Wd,m)+cross(W,cross(W,m)), n2=cross(W,cross(W,m));
    const Vec3 a3=cross(Wdd,m)+cross(Wd,cross(W,m))*2.0+cross(W,cross(Wd,m))+cross(W,cross(W,cross(W,m))), n3=cross(W,cross(W,cross(W,m)));
    std::printf("%s: |m|=%.4e |W'|/|W|=%.4e 1/s |W''|=%.4e (4delta: %.4e)\n   m'  hist %.4e analytic %.4e\n   m'' hist %.4e full %.4e noWdot %.4e\n   m''' hist %.4e full %.4e noWdot %.4e\n",
      first?"first":"second",m.norm(),Wd.norm()/W.norm(),Wdd.norm(),Wdd4.norm(),h.firstDerivative.norm(),a1.norm(),
      h.secondDerivative.norm(),a2.norm(),n2.norm(),h.thirdDerivative.norm(),a3.norm(),n3.norm());
    sumHist2+=h.secondDerivative; sumFull2+=a2; sumNo2+=n2;
    sumHist3+=h.thirdDerivative; sumFull3+=a3; sumNo3+=n3; }
  { // the production path (particleMultipoleRadiation), whatever CREM_M1_ANALYTIC says
    const MutualForces forces=retardedExternalForces(start,hist);
    const ParticleMultipoleRadiation rad=particleMultipoleRadiation(start,forces,hist,false,
        ChargeRadiationReactionModel::stochasticElectricDipole,true);
    std::printf("production path: M1 power %.6e W, torque rate gamma|T|/|mu| %.4e %.4e 1/s\n",rad.magneticDipoleFlux.energy,
      std::abs(firstGyromagneticRatioOf())*rad.firstDipoleTorque.norm()/start.firstDipole.norm(),
      std::abs(secondGyromagneticRatioOf())*rad.secondDipoleTorque.norm()/start.secondDipole.norm());
    std::printf("harness full-m3 torque rate (first): %.4e 1/s\n",
      std::abs(firstGyromagneticRatioOf())*coef*cross(start.firstDipole,sumFull3).norm()/start.firstDipole.norm()); }
  std::printf("M1 power coef|sum m''|^2: hist %.6e  full %.6e  noWdot %.6e W\n",coef*sumHist2.squaredNorm(),coef*sumFull2.squaredNorm(),coef*sumNo2.squaredNorm());
  std::printf("|sum m'''|: hist %.6e  full %.6e  noWdot %.6e\n",sumHist3.norm(),sumFull3.norm(),sumNo3.norm());
}
