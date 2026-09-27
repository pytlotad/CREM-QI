// Audit 233: the refusal probability, on an ensemble rather than the
// twelve trajectories the counter was first read on.  Each refused
// attempt redraws the photon's direction and helicity, so attempts
// should be independent Bernoulli trials with a common p and the
// refusal count per trajectory geometric: P(k) = p^k (1-p).  That is
// testable, and the MLE is p = kbar/(1+kbar).
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    const std::uint64_t lo=argc>1?std::strtoull(argv[1],nullptr,10):1;
    const std::uint64_t hi=argc>2?std::strtoull(argv[2],nullptr,10):25;
    const double budget=argc>3?std::atof(argv[3]):240.0;
    for(std::uint64_t seed=lo;seed<=hi;++seed){
        const CremCollapseEstimate e=estimateCremCollapse(seed,1,budget,
            ChargeRadiationReactionModel::stochasticElectricDipole);
        std::printf("ROW %llu %.6f %zu %llu %llu %llu %d\n",
            (unsigned long long)seed,e.lifetimeSeconds*1e12,
            e.labFramePhotons.size(),e.refusedByCeiling,
            e.refusedByKinematics,e.refusedByRecoil,
            static_cast<int>(e.calibrationOutcome));
    }
}
