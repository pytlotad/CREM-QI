// Audit 256 control: is the net momentum flux physical or quadrature
// residue?  The radiation field's flux through a control sphere is
// radius-independent; the velocity field's contribution is not.  If
// |p| moves with the control radius, or moves when the velocity
// field is dropped, it is not a converged physical quantity.  One
// trajectory, one phase, so only the sampling varies.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A, ecc=0.5;
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const Vec3 mz=Vec3{0.0,0.0,1.0};
    const State start=osculatingPeriapsisState(el,k,
        mz*firstMagneticMoment,mz*secondMagneticMoment,
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-10; acc.maximumDepth=26;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start; bool ok=true;
    for(int i=0;i<50&&ok;++i) ok=engine.advance(st,0.25*period/50);
    if(!ok){ std::printf("nieudane\n"); return 1; }
    std::printf("# jedna trajektoria, tol 1e-10, nd=302, cwierc okresu\n");
    std::printf("%12s %8s %14s %14s %12s\n",
        "R/a0","tylko rad","energia [W]","|ped| [N]","|p|c/E");
    Vec3 ref{};
    for(double rr:{1.0e4,1.0e5,1.0e6,1.0e7}){
        for(bool radOnly:{false,true}){
            FarFieldSampling smp;
            smp.directionCount=302;
            smp.controlRadius=rr*bohrRadius;
            smp.radiationFieldOnly=radOnly;
            const FieldFluxRates r=
                electromagneticFieldFluxRates(st,engine.history(),smp);
            const double pn=r.momentum.norm();
            if(ref.squaredNorm()==0.0&&pn>0.0) ref=r.momentum*(1.0/pn);
            std::printf("%12.0e %8s %14.6e %14.6e %12.4e\n",
                rr,radOnly?"tak":"nie",r.energy,pn,
                pn*c/std::max(r.energy,1e-300));
        }
    }
}
