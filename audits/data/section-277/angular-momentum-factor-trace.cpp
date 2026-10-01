// Audit 278: where does the factor 1.88 come from?
//
// 277 compared the orbital angular-momentum loss with the radiated
// angular momentum and found the field carrying 1.88 times what the
// orbit lost.  Two things were wrong with that comparison and both
// are fixed here.
//
// (1) It compared NORMS of increments, not increments of VECTORS.
//     If the two are not parallel the comparison is meaningless.
//     Here the z components are taken, with the off-axis part
//     reported so the reader can see whether that mattered.
// (2) The run had the radiation reaction DISABLED, so the orbit had
//     no radiative channel to lose angular momentum through at all.
//     Whatever moved L was therefore not radiation -- and comparing
//     it to the radiated flux compares unrelated quantities.
//
// Two controls decide it.  Moments scaled to 1e-06, the same ablation
// 250 used: if dL collapses and dA does not, the L drift was the
// moment sector, which 69 and 226 already record as lacking a
// Lagrangian that makes its force, precession and energy book agree.
// And the circular-orbit identity dL/dt = (dE/dt)/omega, which the
// radiated channel must satisfy whatever the mechanical side does.
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
    const double omega=2.0*pi/period;
    std::printf("# orbita kolowa a_pair, omega=%.6e rad/s\n",omega);
    std::printf("%6s %8s %13s %13s %9s %13s %9s\n","kanal","momenty",
        "dJz/chk","dAz/chk","|dA|/|dL|","dE/om/chk","dAz/(dE/om)");
    for(int channel=0;channel<2;++channel){
        for(int ablate=0;ablate<2;++ablate){
            const double scale=ablate?1.0e-6:1.0;
            const bool para=(channel==0);
            const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*scale);
            const Vec3 m2=Vec3{0,0,para?1.0:-1.0}
                *(secondMagneticMoment*scale);
            const State start=osculatingPeriapsisState(el,kk,m1,m2,
                Vec3{0,0,1},Vec3{1,0,0},0.0);
            ClassicalTrajectoryEngine::Accuracy acc;
            acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
            acc.reactionModel=ChargeRadiationReactionModel::disabled;
            acc.computeOutwardFlux=true;
            acc.farFieldSampling.directionCount=194;
            acc.farFieldSampling.controlRadius=1.0e6*bohrRadius;
            ClassicalTrajectoryEngine engine(start,acc);
            State st=start; bool ok=true;
            // THE MODEL'S OWN total angular momentum: canonical orbital
            // (Darwin and dipole vector potential included) plus the
            // intrinsic spin via the gyromagnetic ratio.  277 and the first
            // version of this probe used the bare mu (r x v) instead, which
            // omits both -- and that omission is the candidate for 1.88.
            const auto orbitalOf=[&](const State& s){
                return noetherAngularMomentum(s);
            };
            Vec3 firstL{},lastL{},firstA{},lastA{};
            double firstE=0.0,lastE=0.0; int taken=0;
            for(int checkpoint=0;checkpoint<4&&ok;++checkpoint){
                for(int i=0;i<25&&ok;++i) ok=engine.advance(st,period/200);
                if(!ok) break;
                const Vec3 L=orbitalOf(st);
                if(taken==0){ firstL=L; firstA=st.radiatedAngularMomentum;
                              firstE=st.radiatedEnergy; }
                lastL=L; lastA=st.radiatedAngularMomentum;
                lastE=st.radiatedEnergy; ++taken;
            }
            if(taken<2){ std::printf("%6s %8s przerwane\n",
                para?"para":"orto",ablate?"off":"ON"); continue; }
            const double n=taken-1;
            const Vec3 dL=(lastL-firstL)*(1.0/n);
            const Vec3 dA=(lastA-firstA)*(1.0/n);
            const double dE=(lastE-firstE)/n;
            std::printf("%6s %8s %13.5e %13.5e %9.3f %13.5e %9.4f\n",
                para?"para":"orto",ablate?"off":"ON",dL.z,dA.z,
                std::abs(dL.z)>0.0?dA.norm()/dL.norm():0.0,
                dE/omega,(dE/omega)!=0.0?dA.z/(dE/omega):0.0);
        }
    }
    std::printf("# |dA|/|dJ| = 1 oznacza, ze ksiega domyka sie na J Noethera\n");
}
