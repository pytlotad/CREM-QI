#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <string>
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
