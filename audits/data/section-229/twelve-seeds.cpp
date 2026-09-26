// Audit 228: what does the orbit DO after the photon that takes a full
// hbar and leaves L = 0?  Classically that is a radial orbit through
// the centre -- the only configuration with contact, whose absence
// 13c calls structural because radiation reaction circularizes.  The
// configuration is the spin magnitude (default since 228) plus
// CREM_AXIAL_SPIN, which the panel records as removing exactly
// 1.000000 hbar and ending the cascade after ONE photon at L = 0.
// "Ends" is not a measurement of what the orbit does.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <cstdlib>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    std::printf("# CREM_AXIAL_SPIN=%s  CREM_NO_SPIN_MAGNITUDE=%s\n",
        std::getenv("CREM_AXIAL_SPIN")?"set":"unset",
        std::getenv("CREM_NO_SPIN_MAGNITUDE")?"set":"unset");
    std::printf("%6s %20s %8s %8s %16s %12s\n",
                "seed","lifetime [s]","photons","outcome",
                "<P> [W]","revolutions");
    for(std::uint64_t seed=5;seed<=12;++seed){
        const CremCollapseEstimate e=estimateCremCollapse(seed,1,240.0,
            ChargeRadiationReactionModel::stochasticElectricDipole);
        std::printf("%6llu %20.12e %8zu %8d %16.6e %12.4f\n",
            (unsigned long long)seed,e.lifetimeSeconds,
            e.labFramePhotons.size(),
            static_cast<int>(e.calibrationOutcome),
            e.meanRadiatedPowerWatts,e.revolutions);
    }
    std::printf("# outcome 0 = ReachedCutoff, 1 = ObservationLimit,"
                " 2 = NumericalFailure\n");
}
