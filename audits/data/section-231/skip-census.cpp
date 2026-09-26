#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(int argc,char** argv){
    gSpinQuantization=true;
    const std::uint64_t seed=argc>1?std::strtoull(argv[1],nullptr,10):1;
    const CremCollapseEstimate e=estimateCremCollapse(seed,1,240.0,
        ChargeRadiationReactionModel::stochasticElectricDipole);
    std::printf("seed %llu  t=%.6e s  photons %zu  rev %.1f\n",
        (unsigned long long)seed,e.lifetimeSeconds,
        e.labFramePhotons.size(),e.revolutions);
}
