// Audit 301: czy rozdzielenie perycentrum para/orto przy malym J0 jest
// wlasnoscia TRAJEKTORII, a nie estymatora perycentrum?  Jedna pelna
// orbita przez simulate() (produkcyjne przygotowanie, przelacznik
// CREM_INITIAL_ANGULAR_MOMENTUM), reakcja WYLACZONA, minimum separacji z
// silnika (adaptacyjne, zapisywane na kazdym kroku -- lekcja 291).
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
    const int ph=std::atoi(argv[1]);
    gRadiationReactionModel=ChargeRadiationReactionModel::disabled;
    SimulationOptions o;
    o.collectFrames=false;
    o.radiatedEnergyBookkeeping=false;
    o.observationTime=1.05*3.0397e-16;      // nieco ponad okres przy n = 1
    o.terminalSeparation=0.05*comptonBarrierRadius;  // nie przerywac przy r*
    const SimulationResult r=simulate(40,ph,o);
    std::printf("%d %.6f %d %.6e\n",ph,r.minimumSeparation/comptonBarrierRadius,
        static_cast<int>(r.outcome),r.elapsedTime);
}
