// Audit 314: jedno przejscie przez perycentrum, minimum separacji z silnika.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
    const std::uint64_t seed=std::strtoull(argv[1],nullptr,10);
    const int ph=std::atoi(argv[2]);
    gRadiationReactionModel=ChargeRadiationReactionModel::disabled;
    SimulationOptions o; o.collectFrames=true; o.frameCount=2;
    o.radiatedEnergyBookkeeping=false;
    o.observationTime=1.05*3.0397e-16;
    o.terminalSeparation=0.05*comptonBarrierRadius;
    const SimulationResult r=simulate(seed,ph,o);
    const auto& f=r.frames.front();
    const Vec3 z{0,0,1};
    const double c12=dot(f.firstDipole,f.secondDipole)/(f.firstDipole.norm()*f.secondDipole.norm());
    std::printf("%llu %d %.6f %d %+.5f %+.5f %+.5f\n",(unsigned long long)seed,ph,
        r.minimumSeparation/comptonBarrierRadius,static_cast<int>(r.outcome),
        c12,dot(f.firstDipole,z)/f.firstDipole.norm(),dot(f.secondDipole,z)/f.secondDipole.norm());
}
