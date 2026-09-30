// Audit 264, part one: every fallback return, exercised.
//
// 263e recorded two unexercised production branches.  This covers the
// machinery side: the six places where the draw gives up and hands back
// a zero vector, which every caller reads as "use the prescribed draw".
// An unexercised give-up path is the kind that returns NaN or spins
// rather than returning zero, and nothing would notice until it did.
//
// Contract under test: the result is EITHER exactly zero OR a unit
// vector.  Never NaN, never a non-unit vector, never a hang.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <random>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    std::mt19937_64 rng(4242);
    std::uniform_real_distribution<double> u(0.0,1.0);
    const auto draw=[&]{ return u(rng); };
    const Vec3 axis{0.0,0.0,1.0};
    const Vec3 peri{1.0,0.0,0.0};
    const auto good=[]{
        AngularPatternMoments m;
        m.a00=4.0*pi; m.a2[2]=(4.0*pi/5.0)*0.5; return m;
    };
    struct Case{const char* name;AngularPatternMoments m;Vec3 ax;Vec3 pe;
                bool expectZero;};
    AngularPatternMoments zeroMonopole; zeroMonopole.a00=0.0;
    AngularPatternMoments negMonopole;  negMonopole.a00=-1.0;
    AngularPatternMoments nanMonopole;
    nanMonopole.a00=std::numeric_limits<double>::quiet_NaN();
    AngularPatternMoments nanMoment=good();
    nanMoment.a2[0]=std::numeric_limits<double>::quiet_NaN();
    // Shape driven far negative over most of the sphere: the monopole is
    // tiny against the quadrupole, so most of the sphere clamps to zero
    // and the acceptance probability collapses.  This is the 64-rejection
    // give-up, and the only one that depends on the RNG.
    AngularPatternMoments starved;
    starved.a00=4.0*pi*1.0e-9;
    starved.a2[2]=(4.0*pi/5.0)*1.0;
    const Case cases[]={
        {"poprawny wzor",           good(),       axis, peri, false},
        {"monopol zerowy",          zeroMonopole, axis, peri, true },
        {"monopol ujemny",          negMonopole,  axis, peri, true },
        {"monopol NaN",             nanMonopole,  axis, peri, true },
        {"moment NaN",              nanMoment,    axis, peri, true },
        {"os zerowa",               good(),       Vec3{}, peri, true },
        {"perycentrum || os",       good(),       axis, axis, true },
        {"perycentrum zerowe",      good(),       axis, Vec3{}, true },
        {"wzor zaglodzony",         starved,      axis, peri, false},
    };
    std::printf("%22s %9s %9s %11s %9s\n",
        "przypadek","zerowych","niezer.","norma max bl","wynik");
    for(const Case& cs:cases){
        int zeros=0,nonzeros=0; double worstNorm=0.0; bool sane=true;
        for(int i=0;i<20000;++i){
            const Vec3 d=drawLabDirectionFromOrbitalPattern(
                cs.m,cs.ax,cs.pe,draw);
            if(!std::isfinite(d.x)||!std::isfinite(d.y)
               ||!std::isfinite(d.z)) { sane=false; break; }
            if(d.squaredNorm()==0.0) ++zeros;
            else {
                ++nonzeros;
                worstNorm=std::max(worstNorm,std::abs(d.norm()-1.0));
            }
        }
        const bool ok=sane
            &&(cs.expectZero?(nonzeros==0):(zeros+nonzeros==20000))
            &&worstNorm<1.0e-12;
        std::printf("%22s %9d %9d %11.2e %9s\n",cs.name,zeros,nonzeros,
            worstNorm,!sane?"NaN!":(ok?"ok":"UWAGA"));
    }
}
