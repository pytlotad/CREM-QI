// Audit 247: which results survive the regulator, and which are the
// regulator.  matchedMomentSoftening() returns magneticDipoleRadius()
// -- about 0.967 r* by default -- and is applied as the charge field's
// Plummer floor wherever that field acts on a moment.  At a periapsis
// of 3.9 r* the vector Plummer form r/(r^2+a^2)^{3/2} is 8.6% weaker
// than the bare Coulomb, so every deep-passage number depends on a
// choice.  CREM_MAGNETIC_RADIUS_SCALE sets it in units of r*.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <cstdlib>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    const char* sc=std::getenv("CREM_MAGNETIC_RADIUS_SCALE");
    std::printf("# CREM_MAGNETIC_RADIUS_SCALE=%s  ->  a = %.6e m"
                " = %.4f r*\n",sc?sc:"(domyslny)",magneticDipoleRadius(),
                magneticDipoleRadius()/comptonBarrierRadius);
    std::printf("# separationFloor = %.4f r*\n",
                separationFloor()/comptonBarrierRadius);
    std::printf("%6s %20s %20s %14s\n",
                "seed","para [s]","ortho [s]","ortho/para-1");
    for(std::uint64_t s=1;s<=3;++s){
        const CremCollapseEstimate p=estimateCremCollapse(s,1,300.0);
        const CremCollapseEstimate o=estimateCremCollapse(s,2,300.0);
        std::printf("%6llu %20.12e %20.12e %14.6e\n",
            (unsigned long long)s,p.lifetimeSeconds,o.lifetimeSeconds,
            o.lifetimeSeconds/p.lifetimeSeconds-1.0);
    }
}
