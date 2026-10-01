// Audit 291d: nadwyzka czasu +3.39e-05 nie zalezy ani od podlogi
// zmiekczania (291c), ani od glebokosci (bit-identycznie 12..18).
// Zostaje zakodowane na sztywno relativeTolerance=1.0e-5, a 3.39e-05
// to ~3.4x tej tolerancji.
// PREDYKCJA ZAREJESTROWANA PRZED URUCHOMIENIEM: nadwyzka spadnie z
// tolerancja.  Jesli stanie na podlodze, to nie blad calkowania.
//
// Sonda uzywa KROKU ADAPTACYJNEGO (0.003 r/v), bo 291c pokazalo, ze
// staly krok przeskakuje dno: przy 0.084 r* jedno 1.34e-22 s pokrywa
// 0.17 r*.  Przy okazji jest to niezalezny pomiar minimum.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double rStar=comptonBarrierRadius;
    const double a=0.5*A;
    const double exact=pi*std::sqrt(a*a*a/kk);
    const double rc=nuclearCutoff;
    const double Ec=std::acos(1.0-rc/a);
    const double tail=std::sqrt(a*a*a/kk)*(Ec-std::sin(Ec));
    const double w1=secondMass/(firstMass+secondMass);
    const double w2=-firstMass/(firstMass+secondMass);
    std::printf("# polokres=%.9e s, odciecie=%.5f r*, ogon=%.4e s\n",
        exact,rc/rStar,tail);
    std::printf("%12s %14s %12s %11s %10s\n","tolerancja","czas+ogon [s]",
        "nadwyzka","min/r*","krokow");
    for(double tol:{1.0e-5,1.0e-7,1.0e-9}){
        State start;
        start.firstPosition=Vec3{A,0,0}*w1;
        start.secondPosition=Vec3{A,0,0}*w2;
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol;
        acc.maximumDepth=30;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start;
        double elapsed=0.0,minimum=1e300,previous=A,previousTime=0.0;
        long steps=0; bool ok=true;
        while(ok&&elapsed<1.2*exact){
            const Vec3 rel=st.firstPosition-st.secondPosition;
            const Vec3 vrel=st.firstVelocity-st.secondVelocity;
            const double r=rel.norm();
            const double v=std::max(vrel.norm(),1.0);
            // krok adaptacyjny: 0.003 czasu przelotu lokalnego
            double step=0.003*r/v;
            if(step>exact/2000.0) step=exact/2000.0;
            previous=r; previousTime=elapsed;
            ok=engine.advance(st,step);
            if(!ok) break;
            elapsed+=step; ++steps;
            const double now=
                (st.firstPosition-st.secondPosition).norm();
            minimum=std::min(minimum,now);
            if(now<=rc){
                // liniowa interpolacja chwili przejscia przez odciecie
                const double frac=(previous-rc)/(previous-now);
                elapsed=previousTime+frac*(elapsed-previousTime);
                break;
            }
        }
        std::printf("%12.0e %14.9e %+12.4e %11.5f %10ld\n",
            tol,elapsed+tail,(elapsed+tail)/exact-1.0,
            minimum/rStar,steps);
    }
}
