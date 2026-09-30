// Audit 261 end-to-end: the wiring, exercised through the real
// estimator.  Two questions, and the first one matters more.
//
//   OFF must be INERT.  With CREM_COMPUTED_EMISSION_PATTERN unset the
//   Cardano draw and the azimuth draw both still happen, patternDirection
//   stays zero, usePatternDirection is false, and both expressions fall
//   back to exactly what they were -- so the run must match the
//   pre-wiring code bit for bit.  Checked against a worktree at HEAD,
//   not against an argument.
//
//   ON must DIFFER, and differ in the emission direction.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const std::uint64_t seed=argc>1?std::strtoull(argv[1],nullptr,10):7;
    const double budget=argc>2?std::atof(argv[2]):25.0;
    const int phenomenon=argc>3?std::atoi(argv[3]):1;
    const CremCollapseEstimate e=estimateCremCollapse(seed,phenomenon,budget,
        ChargeRadiationReactionModel::stochasticElectricDipole);
    double angleSum=0.0; int finite=0;
    for(const LabFramePhoton& p:e.labFramePhotons)
        if(std::isfinite(p.angleFromRecoilAxisRadians)){
            angleSum+=p.angleFromRecoilAxisRadians; ++finite;
        }
    std::printf("ziarno=%llu zjaw=%d zycie=%.12e fotony=%zu katow=%d sumaKatow=%.12e "
        "odmowy=%llu/%llu/%llu wzor=%llu/%llu M1pominiete=%llu\n",
        (unsigned long long)seed,phenomenon,
        e.lifetimeSeconds,e.labFramePhotons.size(),finite,angleSum,
        e.refusedByCeiling,e.refusedByKinematics,e.refusedByRecoil,
        e.patternDirectionDraws,e.prescribedDirectionDraws,
        e.magneticPatternSkips);
    int shown=0; (void)shown;
    if(false)
    for(const LabFramePhoton& p:e.labFramePhotons){
        if(shown++>=3) break;
        std::printf("  foton %d: E=%.10e kat=%.10e beta=%.6e\n",
            shown,p.energyJoules,p.angleFromRecoilAxisRadians,p.sourceBeta);
    }
}
