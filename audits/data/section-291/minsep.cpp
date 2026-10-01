// Audit 291: co wyznacza minimum 0.084 r* na orbicie radialnej?
//
// Przy L=0 orbita radialna powinna przechodzic przez zero.  285
// przypisalo minimum zmiekczaniu (blednie), skan przypisal je krokowi
// (tez blednie -- jest zbiezne i niezalezne od podlogi).
//
// Kandydat wiodacy: orbita NIE jest L=0.  Startuje z L=0 dokladnie,
// ale sily dipolowe sa NIECENTRALNE, wiec moment pedu narasta podczas
// spadku.  Dla malego l perycentrum to r_p = l^2/(2k), wiec
// r_p = 0.084 r* wymaga l = sqrt(2 k r_p) = 0.0175 hbar/mu.
//
// Ablacja rozstrzygajaca: przy momentach DOKLADNIE zerowych nie ma
// zadnej sily niecentralnej (Coulomb centralny, Darwin dla ruchu
// radialnego nie daje momentu, bo v jest rownolegle do r), wiec L
// musi zostac zerem i orbita MUSI przejsc przez zero.
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
    const double halfPeriod=pi*std::sqrt(std::pow(0.5*A,3)/kk);
    std::printf("# k=%.6e m^3/s^2, hbar/mu=%.6e m^2/s\n",
        kk,hbar/mu_);
    std::printf("%10s %11s %13s %13s %13s\n","momenty","min/r*",
        "l max [m^2/s]","l max [hbar]","r_p(l) /r*");
    for(double scale:{0.0,1.0e-6,1.0e-3,1.0}){
        const double w1=secondMass/(firstMass+secondMass);
        const double w2=-firstMass/(firstMass+secondMass);
        State start;
        start.firstPosition=Vec3{A,0,0}*w1;
        start.secondPosition=Vec3{A,0,0}*w2;
        const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*scale);
        start.firstDipole=m1;
        start.secondDipole=m1*(secondMagneticMoment/firstMagneticMoment);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-9; acc.maximumDepth=30;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        double minimumSeparation=1e300,maximumAngular=0.0;
        const double step=halfPeriod/400000.0;
        double elapsed=0.0;
        while(ok&&elapsed<1.2*halfPeriod){
            ok=engine.advance(st,step);
            if(!ok) break;
            elapsed+=step;
            const Vec3 rel=st.firstPosition-st.secondPosition;
            const Vec3 vrel=st.firstVelocity-st.secondVelocity;
            minimumSeparation=std::min(minimumSeparation,rel.norm());
            maximumAngular=std::max(maximumAngular,cross(rel,vrel).norm());
        }
        const double impliedPeriapsis=
            maximumAngular*maximumAngular/(2.0*kk);
        std::printf("%10.0e %11.5f %13.4e %13.4e %13.5f\n",scale,
            minimumSeparation/rStar,maximumAngular,
            maximumAngular/(hbar/mu_),impliedPeriapsis/rStar);
        if(!ok) std::printf("           (przerwane)\n");
    }
    std::printf("# r_p(l) ma odtworzyc min/r*, jesli minimum jest"
        " nabytym momentem pedu\n");
}
