// Stop cause of single trajectories of the audit 327 ensemble B, same seeds
// as runCremCollapseExperiment (splitMix64(master + index)).
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
static const char* name(CollapseStopCause c){
    switch(c){case CollapseStopCause::ComptonBarrier:return "barrier";
    case CollapseStopCause::RetardationLimit:return "retardation";
    case CollapseStopCause::GroundStateFloor:return "floor";
    case CollapseStopCause::EmissionChannelClosed:return "closed";
    default:return "none";}
}
int main(int argc,char** argv){
    gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
    gMeasureCollapseTransit=false; gInitialPrincipalLevel=1;
    gMicrocanonicalStart=std::getenv("CREM_MICROCANONICAL_START")!=nullptr;
    for(int a=1;a<argc;++a){
        const int index=std::atoi(argv[a]);
        for(int ph=1;ph<=(std::getenv("ONLY_PARA")?1:2);++ph){
            const CremCollapseEstimate e=estimateCremCollapse(splitMix64(42ULL+index),ph,600.0);
            std::printf("idx %2d %s: t = %.6f ps, stop %s, photons %llu, Lcap events %llu largest %.3e hbar unresolved %.3e hbar\n",index,ph==1?"para ":"ortho",
                e.calibrationSecondsLab*1e12,name(e.stopCause),(unsigned long long)e.emittedPhotonCount,e.orbitalCapEvents,e.orbitalCapLargestExcessHbar,e.orbitalCapUnresolvedHbar);
        }
    }
}
