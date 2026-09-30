// Audit 265 control: is the inversion sampler's polar marginal biased?
// The null test at 4e+05 samples returned 5<P2> = 0.4919 where the
// rejection sampler had 0.5008.  Two values below 0.5 is either chance
// or a 1.6 percent bias, and the difference matters: 5<P2> IS c2/c0,
// the quantity production assumes to be exactly 1/2.  Settled by raising
// the sample count tenfold and running both samplers on independent
// streams.  Expected standard error at 4e+06 is 5*0.45/2000 = 1.1e-03.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <random>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    AngularPatternMoments m;
    m.a00=4.0*pi*(4.0/3.0);
    m.a2[2]=(4.0*pi/5.0)*(2.0/3.0);   // dokladnie (1+cos^2)
    std::printf("# wzor (1+cos^2): c2/c0 = 0.5 dokladnie, 5<P2> tego szuka\n");
    std::printf("%16s %10s %12s %12s\n","sampler","probki","5<P2>","<cos2phi>");
    for(int variant=0;variant<2;++variant){
        for(int rep=0;rep<3;++rep){
            std::mt19937_64 rng(777+rep*31);
            std::uniform_real_distribution<double> u(0.0,1.0);
            const auto dr=[&]{ return u(rng); };
            const int n=4000000;
            double sp=0.0,sc=0.0; long long got=0;
            for(int i=0;i<n;++i){
                const Vec3 d=variant==0
                    ?drawDirectionFromPatternByInversion(m,dr)
                    :drawDirectionFromPatternByRejection(m,dr);
                if(d.squaredNorm()==0.0) continue;
                ++got;
                sp+=0.5*(3.0*d.z*d.z-1.0);
                sc+=std::cos(2.0*std::atan2(d.y,d.x));
            }
            std::printf("%16s %10lld %12.6f %12.6f\n",
                variant==0?"odwrotnosciowy":"odrzuceniowy",got,
                5.0*sp/got,sc/got);
        }
    }
}
