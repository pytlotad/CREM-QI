#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <vector>
#include <algorithm>
struct Prep { double e2, J0, mm, vr; };
static Prep prepare(std::uint64_t seed) {
    State first{}; bool got=false;
    SimulationOptions o; o.collectFrames=false; o.radiatedEnergyBookkeeping=false;
    o.observationTime=1e-20;
    o.stepReady=[&](const State& s){ if(!got){first=s; got=true;} };
    o.stopRequested=[&]{ return got; };
    simulate(seed,1,o);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const Vec3 r=first.firstPosition-first.secondPosition, v=first.firstVelocity-first.secondVelocity;
    const double E=0.5*mu*dot(v,v)-pairCoulombStrength/r.norm();
    const double L=mu*cross(r,v).norm();
    const double e2=1.0+2.0*E*L*L/(mu*pairCoulombStrength*pairCoulombStrength);
    const double Lc=std::sqrt(mu*pairCoulombStrength*r.norm());
    return {e2, L/Lc, (first.firstDipole+first.secondDipole).norm()/first.firstDipole.norm(), dot(r,v)/(r.norm()*v.norm())};
}
int main() {
    gRadiationReactionModel=ChargeRadiationReactionModel::disabled;
    const int N=2000; std::vector<double> e2s; int paired=0, inward=0; double maxe2circ=0;
    for(int i=0;i<N;++i) {
        gMicrocanonicalStart=false; const Prep c=prepare(1000+i);
        gMicrocanonicalStart=true;  const Prep m=prepare(1000+i);
        maxe2circ=std::max(maxe2circ,std::abs(c.e2));
        if(c.mm==m.mm) ++paired;
        if(m.vr<0) ++inward;
        e2s.push_back(m.e2);
    }
    std::sort(e2s.begin(),e2s.end()); double D=0;
    for(int i=0;i<N;++i) D=std::max({D,std::abs(e2s[i]-(double)i/N),std::abs(e2s[i]-(double)(i+1)/N)});
    double mean=0; for(double x:e2s) mean+=x; mean/=N;
    std::printf("N=%d: circle max|e^2| %.2e; micro <e^2> %.4f (uniform 0.5, sd %.4f); KS D=%.4f (5%% crit %.4f)\n",N,maxe2circ,mean,std::sqrt(1.0/12/N),D,1.36/std::sqrt(N));
    std::printf("moments identical circle vs micro: %d/%d; inward radial: %d/%d\n",paired,N,inward,N);
    std::printf("quartiles e^2: %.3f %.3f %.3f; min %.2e max %.6f\n",e2s[N/4],e2s[N/2],e2s[3*N/4],e2s[0],e2s[N-1]);
}
