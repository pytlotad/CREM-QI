// Audit 291e: nadwyzka czasu spadku +3.61e-05 nie jest ani zmiekczaniem
// (291c: niezalezna od podlogi do 16x), ani bledem calkowania (291d:
// bit-identyczna od 1e-5 do 1e-9, i glebokosc 12..18 bez zmiany).
// KANDYDAT: poprawka relatywistyczna.  Orbita siega 0.83 c, a wzor
// pi*sqrt(a^3/k) jest nierelatywistyczny, i znak sie zgadza (wieksza
// bezwladnosc = dluzszy spadek).
// PREDYKCJA ZAREJESTROWANA PRZED URUCHOMIENIEM: efekt siedzi WYLACZNIE
// w glebokim obszarze, ktorego czas bezwzgledny jest ustalony przez
// odciecie i k, a polokres rosnie jak A^(3/2).  Nadwyzka ma wiec malec
// jak A^(-3/2): dwukrotne apocentrum = spadek 2.83x, czterokrotne = 8x.
// Jesli nadwyzka nie zalezy od apocentrum, efekt jest rozlozony po
// calej orbicie i kandydat pada.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A0=pairBohrRadius({electron,positron});
    const double rStar=comptonBarrierRadius;
    const double rc=nuclearCutoff;
    const double w1=secondMass/(firstMass+secondMass);
    const double w2=-firstMass/(firstMass+secondMass);
    std::printf("%10s %16s %12s %10s %12s\n","apo/a_pair","polokres [s]",
        "nadwyzka","vmax/c","nadw*apo^1.5");
    for(double f:{1.0,2.0,4.0,8.0}){
        const double A=A0*f;
        const double a=0.5*A;
        const double exact=pi*std::sqrt(a*a*a/kk);
        const double Ec=std::acos(1.0-rc/a);
        const double tail=std::sqrt(a*a*a/kk)*(Ec-std::sin(Ec));
        State start;
        start.firstPosition=Vec3{A,0,0}*w1;
        start.secondPosition=Vec3{A,0,0}*w2;
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-7;
        acc.maximumDepth=30;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start;
        double elapsed=0.0,previous=A,previousTime=0.0,vmax=0.0;
        bool ok=true;
        while(ok&&elapsed<1.2*exact){
            const Vec3 rel=st.firstPosition-st.secondPosition;
            const Vec3 vrel=st.firstVelocity-st.secondVelocity;
            const double r=rel.norm();
            const double v=std::max(vrel.norm(),1.0);
            vmax=std::max(vmax,vrel.norm());
            double step=0.003*r/v;
            if(step>exact/2000.0) step=exact/2000.0;
            previous=r; previousTime=elapsed;
            ok=engine.advance(st,step);
            if(!ok) break;
            elapsed+=step;
            const double now=
                (st.firstPosition-st.secondPosition).norm();
            if(now<=rc){
                const double frac=(previous-rc)/(previous-now);
                elapsed=previousTime+frac*(elapsed-previousTime);
                break;
            }
        }
        const double excess=(elapsed+tail)/exact-1.0;
        std::printf("%10.0f %16.9e %+12.4e %10.5f %12.4e\n",
            f,exact,excess,vmax/positronium::kinematics::kSpeedOfLight,
            excess*std::pow(f,1.5));
    }
    std::printf("# kolumna nadw*apo^1.5 ma byc STALA, jesli predykcja"
        " jest trafiona\n");
}
