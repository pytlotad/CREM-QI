// Audit 245: the retarded-history startup transient, measured the
// way the ZPF one was.
//
// causalInitialHistory fabricates the past before t = 0, so for the
// first few light-crossing times the retarded forces are built from
// an extrapolation.  The objection is that no analogue of the ZPF's
// 750-orbit settling measurement exists for it, and that any window
// shorter than tens of orbits from the start is therefore suspect.
//
// The quantity tracked is the one 215 showed to be the collapse
// time's integrand: orbitalRadiatedEnergy accumulated over each
// successive orbit, from a fresh start, reaction off so nothing
// physical changes it between orbits.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const double sep0=a;
    std::printf("# 0.25 a_pair: T = %.4e s, r/c = %.4e s = T/%.0f\n",
                period,sep0/c,period/(sep0/c));
    std::printf("# retentionTime = 4r/c = T/%.0f\n",
                period/(4.0*sep0/c));
    const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}*firstMagneticMoment;
    const State start=osculatingPeriapsisState(
        el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-8; acc.maximumDepth=22;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State s=start;
    const int orbits=40,perOrbit=128;
    double previous=-conservativeParticleEnergy(start);
    std::vector<double> perOrbitLoss;
    for(int o=0;o<orbits;++o){
        bool ok=true;
        for(int i=0;i<perOrbit&&ok;++i)
            ok=engine.advance(s,period/perOrbit);
        if(!ok){ std::printf("# przerwane na obiegu %d\n",o+1); break; }
        const double now=-conservativeParticleEnergy(s);
        perOrbitLoss.push_back(now-previous);
        previous=now;
    }
    if(perOrbitLoss.size()<10){ std::printf("za malo obiegow\n"); return 1; }
    // The asymptote: mean of the last ten orbits.
    double tail=0.0;
    for(std::size_t i=perOrbitLoss.size()-8;i<perOrbitLoss.size();++i)
        tail+=perOrbitLoss[i];
    tail/=8.0;
    std::printf("\n# dryf energii zachowawczej na obieg (reakcja OFF,\n");
    std::printf("# sily retardowane ON), wzgledem sredniej z 8 ostatnich\n");
    std::printf("%7s %18s %14s\n","obieg","dE [J]","odchylka");
    for(std::size_t i=0;i<perOrbitLoss.size();++i){
        const bool show=i<12||i==14||i==19||i==29||i==39
                       ||i+1==perOrbitLoss.size();
        if(show)
            std::printf("%7zu %18.10e %14.3e\n",
                        i+1,perOrbitLoss[i],perOrbitLoss[i]/tail-1.0);
    }
    std::printf("\n# asymptota (10 ostatnich) = %.10e J\n",tail);
}
