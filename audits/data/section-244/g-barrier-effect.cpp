#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    std::printf("# g(e-) = %.14f, r* = %.10e m, |mu| = %.10e J/T\n",
        electron.gFactor,comptonBarrierRadius,firstMagneticMoment);
    for(std::uint64_t s=1;s<=3;++s){
        const CremCollapseEstimate e=estimateCremCollapse(s,1,240.0);
        std::printf("seed %llu  t = %.12e s\n",(unsigned long long)s,
                    e.lifetimeSeconds);
    }
}
