#include "modules/crem_trajectory.hpp"
#include <cstdio>
static Frame first(std::uint64_t seed) {
    Frame f{}; bool got=false;
    SimulationOptions o; o.collectFrames=false; o.radiatedEnergyBookkeeping=false; o.observationTime=1e-20;
    o.frameReady=[&](const Frame& fr){ if(!got){f=fr; got=true;} };
    simulate(seed,1,o); return f;
}
int main() {
    gRadiationReactionModel=ChargeRadiationReactionModel::disabled;
    int same=0; double worst=0;
    for(int i=0;i<500;++i) {
        gMicrocanonicalStart=false; const Frame c=first(1000+i);
        gMicrocanonicalStart=true;  const Frame m=first(1000+i);
        const double d=(c.firstDipole-m.firstDipole).norm()+(c.secondDipole-m.secondDipole).norm();
        if(d==0.0 && c.time==m.time) ++same; worst=std::max(worst,d/c.firstDipole.norm());
    }
    std::printf("t=0 moments bit-identical circle vs micro: %d/500, worst rel diff %.2e\n",same,worst);
}
