// Audit 232: the refusal counter 231f said the code did not expose.
// 231c inferred that the three lifetime levels are zero, one and two
// REFUSED emission attempts, each costing one full deterministic
// hazard unit of 186.72 ps.  Inferred from the level structure and
// the code path, never counted.  Counted here.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    std::printf("%6s %14s %8s %9s %9s %9s %7s\n",
        "seed","lifetime [ps]","photons","ceiling","kinematic","recoil",
        "level");
    for(std::uint64_t seed=1;seed<=12;++seed){
        const CremCollapseEstimate e=estimateCremCollapse(seed,1,240.0,
            ChargeRadiationReactionModel::stochasticElectricDipole);
        std::printf("%6llu %14.3f %8zu %9llu %9llu %9llu %7.0f\n",
            (unsigned long long)seed,e.lifetimeSeconds*1e12,
            e.labFramePhotons.size(),e.refusedByCeiling,
            e.refusedByKinematics,e.refusedByRecoil,
            e.lifetimeSeconds*1e12/186.72);
    }
    std::printf("# 231c predicts level = 1 + total refusals.\n");
}
