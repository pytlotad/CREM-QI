// Audit 228: does making the spin magnitude the default change the
// production path at all, and does it censor?  213b measured zero
// stochastic photons under the default reaction model, and the switch
// only acts where a photon is emitted, so the first question is
// whether anything moves.  Both reaction models are run.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <cstdlib>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    std::printf("# CREM_NO_SPIN_MAGNITUDE=%s\n",
        std::getenv("CREM_NO_SPIN_MAGNITUDE")?"set (old default)":"unset (new)");
    std::printf("%6s %10s %20s %8s %10s\n",
                "seed","model","lifetime [s]","photons","outcome");
    for(std::uint64_t seed=1;seed<=4;++seed)
        for(int m=0;m<2;++m){
            const ChargeRadiationReactionModel model=m==0
                ?ChargeRadiationReactionModel::individualLandauLifshitz
                :ChargeRadiationReactionModel::stochasticElectricDipole;
            const CremCollapseEstimate e=
                estimateCremCollapse(seed,1,300.0,model);
            std::printf("%6llu %10s %20.12e %8zu %10d\n",
                (unsigned long long)seed,m==0?"LL":"stochastic",
                e.lifetimeSeconds,e.labFramePhotons.size(),
                static_cast<int>(e.calibrationOutcome));
        }
    std::printf("# outcome 0 = ReachedCutoff, 1 = ObservationLimit,"
                " 2 = NumericalFailure\n");
}
