// Audit 259 null test for drawDirectionFromPattern: feed it the moments
// of an exactly known pattern and check the draws come back with that
// pattern.  For f ~ (1+z^2) the P2 coefficient is c2/c0 = 1/2, and for
// a normalized pdf <P2(z)> = (c2/c0)/5 = 0.1, so 5<P2> recovers 0.5.
// Second case: a pattern with a genuine m=2 azimuthal harmonic, which
// is the whole reason this exists -- checked through <cos 2phi>.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <random>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> u(0.0,1.0);
    const auto draw=[&]{ return u(rng); };
    const double f4pi=4.0*pi;
    std::printf("%34s %12s %12s %10s\n",
        "wzor","5<P2>","<cos2phi>","odrzucone");
    struct Case{const char* n;double c2;double cm2;double oczP2;double oczC2;};
    for(const Case& cs:{
            Case{"(1+cos^2), oczekiwane 0.5 i 0",0.5,0.0,0.5,0.0},
            Case{"izotropowy, oczekiwane 0 i 0",0.0,0.0,0.0,0.0},
            Case{"z m=2: c2=0.5, cm2=0.3",0.5,0.3,0.5,0.0}}){
        AngularPatternMoments m;
        m.a00=f4pi*1.0;
        m.a2[2]=(f4pi/5.0)*cs.c2;
        m.a2[4]=(f4pi/15.0)*cs.cm2;
        const int n=400000;
        double sp2=0.0,sc2=0.0; int got=0,tries=0;
        for(int i=0;i<n;++i){
            ++tries;
            const Vec3 d=drawDirectionFromPattern(m,draw);
            if(d.squaredNorm()==0.0) continue;
            ++got;
            sp2+=0.5*(3.0*d.z*d.z-1.0);
            const double phi=std::atan2(d.y,d.x);
            sc2+=std::cos(2.0*phi);
        }
        std::printf("%34s %12.5f %12.5f %10.4f\n",cs.n,
            5.0*sp2/got,sc2/got,1.0-double(got)/tries);
    }
    // The m=2 case checked against its own analytic mean.
    // For f = 1 + c2 P2(z) + cm2 (x^2-y^2)/2, <cos2phi> over the pdf is
    // cm2 * <(x^2-y^2)cos2phi/2> / <f> ; (x^2-y^2)/2 = (1-z^2)cos2phi/2,
    // so the numerator is cm2/2 * <(1-z^2) cos^2 2phi> = cm2/2 * 1/2 *
    // <1-z^2> over the l=0 part, i.e. cm2/4 * (2/3) = cm2/6.
    std::printf("# dla cm2=0.3 oczekiwane <cos2phi> = cm2/6 = %.5f\n",0.3/6.0);
}
