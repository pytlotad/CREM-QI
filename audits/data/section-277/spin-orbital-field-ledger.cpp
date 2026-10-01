// Audit 277: does J = S + L close once the FIELD term is included?
//
// Findings brought in from outside this register report a determi-
// nistic spin-to-orbit exchange in ortho -- dS = -6.478e-07 hbar and
// dL = +8.158e-07 hbar at each of eleven checkpoints, constant to
// four digits -- and record that J = S + L is then NOT conserved.
// The arithmetic of those two numbers gives dJ = +1.680e-07 hbar per
// checkpoint: angular momentum appearing from nowhere.
//
// No emission model can be built on that.  But there is an obvious
// missing term: the FIELD carries angular momentum too, and the
// engine already accumulates it as State::radiatedAngularMomentum.
// If dJ = -dA the ledger closes and the exchange and the radiation
// are two views of one thing -- which is exactly the shape an
// emission model would need.  If it does not close, there is a
// genuine source and the model has to wait.
//
// Spin from the proper moments by the gyromagnetic relation
// mu = g (q/2m) S, so S = (2m/(g q)) mu -- signed, so the electron's
// spin runs antiparallel to its moment and the positron's parallel.
// That sign is the whole point for para versus ortho and is taken
// from the charge, not assumed.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const OsculatingElements el{-kk/(2.0*A),std::sqrt(kk*A)};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const auto spinOf=[](const Vec3& moment,double charge,double mass){
        return moment*(2.0*mass/(electronGFactor*charge));
    };
    for(int channel=0;channel<2;++channel){
        const bool para=(channel==0);
        // Parallel moments = singlet = para (the inversion recorded in
        // 213b); antiparallel moments = ortho.
        const Vec3 m1=Vec3{0.0,0.0,1.0}*firstMagneticMoment;
        const Vec3 m2=Vec3{0.0,0.0,para?1.0:-1.0}*secondMagneticMoment;
        const State start=osculatingPeriapsisState(el,kk,m1,m2,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        const Vec3 s0=spinOf(start.firstDipole,firstCharge,firstMass)
                     +spinOf(start.secondDipole,secondCharge,secondMass);
        std::printf("\n=== %s: spin wypadkowy na starcie %.4f hbar ===\n",
            para?"PARA":"ORTO",s0.norm()/hbar);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=true;
        acc.farFieldSampling.directionCount=194;
        acc.farFieldSampling.controlRadius=1.0e6*bohrRadius;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        std::printf("%6s %13s %13s %13s %13s %13s\n",
            "chkpt","S [hbar]","L [hbar]","J=S+L","A pola","J+A");
        Vec3 previousTotal{}; bool havePrevious=false;
        for(int checkpoint=0;checkpoint<8&&ok;++checkpoint){
            for(int i=0;i<25&&ok;++i) ok=engine.advance(st,period/200);
            if(!ok) break;
            const Vec3 spin=spinOf(st.firstDipole,firstCharge,firstMass)
                +spinOf(st.secondDipole,secondCharge,secondMass);
            const Vec3 rel=st.firstPosition-st.secondPosition;
            const Vec3 vrel=st.firstVelocity-st.secondVelocity;
            const Vec3 orbital=cross(rel,vrel)*mu_;
            const Vec3 total=spin+orbital;
            const Vec3 withField=total+st.radiatedAngularMomentum;
            std::printf("%6d %13.6e %13.6e %13.6e %13.6e %13.6e\n",
                checkpoint+1,spin.norm()/hbar,orbital.norm()/hbar,
                total.norm()/hbar,st.radiatedAngularMomentum.norm()/hbar,
                withField.norm()/hbar);
            if(havePrevious){
                const Vec3 dJ=total-previousTotal;
                std::printf("%6s %13s %13s dJ=%10.3e %13s dJ+dA ->\n",
                    "","","",dJ.norm()/hbar,"");
            }
            previousTotal=total; havePrevious=true;
        }
        if(!ok) std::printf("  przerwane\n");
    }
}
