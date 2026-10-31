// Audit 661 = audit 660 probe (copy of audit 623's p623.cpp + a "BHAT" line with the seeded field direction): C, P, CP image runs.  env MODE = traj | est, PAIR, PHEN (1 para, 2 ortho), NLEV, LFRAC, QSPIN (1 = quantized
// spins), BEXT [T] along a seeded direction times BSIGN, OBS [s] (traj window), FRAMES; CREM_SYMMETRY_TEST (P, M, PM) is
// read by the engine.  argv: run indices.  Output in %a (hex floats) so image runs can be compared bit for bit.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
static double envd(const char* k,double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static void pv(const char* tag,const Vec3& v){ std::printf(" %s %a %a %a",tag,v.x,v.y,v.z); }
int main(int argc,char**argv){
  if(const char* p=std::getenv("PAIR")) applyPairFromOption(p);
  const std::string mode=std::getenv("MODE")?std::getenv("MODE"):"traj";
  const int phen=(int)envd("PHEN",1);
  gSpinQuantization=envd("QSPIN",1)!=0.0;
  gInitialAngularMomentumFraction=envd("LFRAC",0.75);
  for(int a=1;a<argc;++a){ const int i=std::atoi(argv[a]);
    const double nl=envd("NLEV",2.0); gInitialSeparationOverride=nl*nl*pairBohrRadius(activePair);
    { std::mt19937_64 r(623000ULL+i); std::uniform_real_distribution<double> u(0,1); const double cz=2*u(r)-1, ph=2*pi*u(r), sz=std::sqrt(1-cz*cz);
      gExternalMagneticField=Vec3{sz*std::cos(ph),sz*std::sin(ph),cz}*(envd("BEXT",0.0)*envd("BSIGN",1.0)); }
    const std::uint64_t seed=splitMix64(623000ULL+i);
    if(mode=="traj"){
      SimulationOptions o; o.frameCount=(int)envd("FRAMES",200); o.observationTime=envd("OBS",2.0e-14);
      const SimulationResult r=simulate(seed,phen,o);
      std::printf("RUN %d outcome %d frames %zu minsep %a elapsed %a radiated %a\n",i,(int)r.outcome,r.frames.size(),
                  r.minimumSeparation,r.elapsedTime,r.finalRadiatedEnergy);
      for(size_t k=0;k<r.frames.size();++k){ const Frame& f=r.frames[k];
        std::printf("F %d %zu t %a",i,k,f.time); pv("r1",f.first); pv("r2",f.second); pv("m1",f.firstDipole); pv("m2",f.secondDipole);
        pv("P",f.noetherMomentum); pv("J",f.noetherAngularMomentum); pv("Prad",f.radiatedMomentum); pv("Jrad",f.radiatedAngularMomentum);
        std::printf(" Emech %a Erad %a rad %a E1 %a E2 %a\n",f.mechanicalEnergy,f.radiatedEnergy,f.radius,
                    f.firstMechanicalEnergy,f.secondMechanicalEnergy); }
    } else {
      gDeterministicEmission=true;
      gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole; gMeasureCollapseTransit=false;
      std::fflush(stdout); std::fprintf(stderr,"START idx %d\n",i); { const double Bn=gExternalMagneticField.norm(); if(Bn>0.0){ const Vec3 b=gExternalMagneticField*(1.0/Bn); std::fprintf(stderr,"BHAT %.15e %.15e %.15e\n",b.x,b.y,b.z);} } std::fflush(stderr);
      const CremCollapseEstimate e=estimateCremCollapse(seed,phen,envd("BUDGET",900.0));
      std::fprintf(stderr,"END idx %d stop %d photons %lld tau %a born %a chosen %d vcm %a %a %a\n",i,(int)e.stopCause,
                   (long long)e.emittedPhotonCount,e.lifetimeSeconds,e.bornSingletProbability,int(e.bornSingletChosen),
                   e.centreOfMassVelocityAtStop.x,e.centreOfMassVelocityAtStop.y,e.centreOfMassVelocityAtStop.z);
      std::fflush(stderr);
    }
  }
}
