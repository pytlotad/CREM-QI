// Audit 241: is the unsourced 54% of 122d the RADIATED momentum?
// 122d's balance was F1+F2 = d(P_field + P_hidden)/dt with no
// radiation term in it at all.  240f measured P_field = -P_hidden to
// one part in 6300, so their sum carries nothing and the balance
// should reduce to F1+F2 = -(radiated momentum flux).  Both sides
// are available: retardedExternalForces for the left,
// electromagneticFieldFluxRates(...).momentum for the right.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <cstdlib>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    std::printf("# CREM_NO_HIDDEN_MOMENTUM_FORCE=%s\n",
        std::getenv("CREM_NO_HIDDEN_MOMENTUM_FORCE")?"set":"unset");
    std::printf("%9s %14s %14s %12s %12s\n",
                "r/a_pair","|F1+F2| [N]","F_Coulomb [N]","niezbil./FC",
                "|dPrad|/FC");
    for(double f:{0.25,0.1,0.05,0.02,0.01,0.005}){
        const double a=f*A;
        const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}*firstMagneticMoment;
        const State start=osculatingPeriapsisState(
            el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start; bool ok=true;
        const int n=400;
        for(int i=0;i<n&&ok;++i) ok=engine.advance(s,0.5*period/n);
        if(!ok){ std::printf("%9.3f  failed\n",f); continue; }
        const MutualForces F=retardedExternalForces(s,engine.history());
        const Vec3 net=F.first+F.second;
        FarFieldSampling sampling; sampling.directionCount=194;
        const FieldFluxRates flux=
            electromagneticFieldFluxRates(s,engine.history(),sampling);
        const double cosang=net.norm()>0&&flux.momentum.norm()>0
            ?dot(net,flux.momentum)/(net.norm()*flux.momentum.norm()):0.0;
        const double sep=(s.firstPosition-s.secondPosition).norm();
        const double FC=pairCoulombStrength/(sep*sep);
        (void)cosang;
        std::printf("%9.3f %14.6e %14.6e %12.3e %12.3e\n",
            f,net.norm(),FC,net.norm()/FC,flux.momentum.norm()/FC);
    }
    std::printf("# jesli bilans domyka sie, stosunek = 1 i cos = -1\n");
}
