// Audit 478: one close e+ e- pass on the engine's mechanical path (runMechanicalTrajectory: retardation, Thomas-BMT,
// dipole coupling, radiation reaction) -- change of the spins and of their mutual angle vs the closest approach.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
int main(int argc,char**argv){
  applyPairFromOption("electron,positron");
  const double bfm=std::atof(argv[1]); const int N=std::atoi(argv[2]); const double EkeV=std::atof(argv[3]);
  const double fm=1e-15, D0=5000*fm, mu=firstMass*secondMass/(firstMass+secondMass);
  const double v=c*std::sqrt(1-1/std::pow(1+EkeV*1e3*1.602176634e-19/(mu*c*c),2));   // relative speed at infinity (rough)
  const double s1=firstCharge<0?-1:1, s2=secondCharge<0?-1:1;
  for(int i=0;i<N;++i){
    std::mt19937_64 rng(478000+i); std::normal_distribution<double> g(0,1);
    auto unit=[&]{ Vec3 u{g(rng),g(rng),g(rng)}; return u*(1.0/u.norm()); };
    const Vec3 u1=unit(), u2=unit();
    State s{}; const Vec3 rel{-D0,bfm*fm,0}, vrel{v,0,0};
    s.firstPosition=rel*(secondMass/(firstMass+secondMass)); s.secondPosition=rel*(-firstMass/(firstMass+secondMass));
    s.firstVelocity=vrel*(secondMass/(firstMass+secondMass)); s.secondVelocity=vrel*(-firstMass/(firstMass+secondMass));
    s.firstDipole=u1*(s1*firstMagneticMoment); s.secondDipole=u2*(s2*secondMagneticMoment);
    SimulationOptions opt; opt.collectFrames=false; opt.radiatedEnergyBookkeeping=false;
    const auto run=runMechanicalTrajectory(s,2.0*D0/v*1.2,0.5*fm,opt);
    const State& f=run.finalState;
    const Vec3 f1=f.firstDipole*(s1/f.firstDipole.norm()), f2=f.secondDipole*(s2/f.secondDipole.norm());
    const double c0=dot(u1,u2), cf=dot(f1,f2), rot1=std::acos(std::clamp(dot(u1,f1),-1.0,1.0)), rot2=std::acos(std::clamp(dot(u2,f2),-1.0,1.0));
    std::printf("b %.1f fm  rmin %.4f fm  beta_max %.4f  outcome %d stop %d  cos12 %+.6f -> %+.6f (d %+.3e)  rot1 %.4e rot2 %.4e rad  sep_end %.1f fm\n",
      bfm,run.minimumSeparation/fm,run.maximumBeta,int(run.outcome),int(run.stopReason),c0,cf,cf-c0,rot1,rot2,separation(f)/fm);
    std::fflush(stdout);
  }
}
