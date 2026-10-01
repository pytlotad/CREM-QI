// Audit 278: where do the remaining 3.1 percent come from?
//
// 277e offered transit as the candidate: the control sphere sits 580
// orbital periods away, so a steady state might not be reached.  That
// candidate is REFUTED by reading the code, before any run.
// farZoneChargeField takes the emission time as wavefrontTime +
// n.(x_src - centre)/c with wavefrontTime = state.time -
// sourceExtent/c: the R/c delay is absent by construction, which is
// the standard far-zone treatment.  Were it otherwise, a sphere 580
// periods away and a run half a period long would give identically
// zero flux instead of a linearly growing one.
//
// So the remaining candidates are the sampling itself and the time
// quadrature, and they are scanned here: control radius, direction
// count, whether the velocity field is included, and the number of
// checkpoints the rate is averaged over.
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
    std::printf("# przewidywanie: dMz == luka dJz+dAz, iloraz 1.0\n");
    std::printf("%9s %13s %13s %12s %12s %9s\n",
        "obiegi","dJz/chk","dAz/chk","luka","dMz/chk","dM/luka");
    // SCHOTT TEST.  The order-reduced Landau-Lifshitz force omits the
    // Schott term, a total derivative that exchanges angular momentum
    // reversibly with the bound field.  Over a whole number of orbital
    // periods it averages to zero; over half a period -- which is what
    // every run above was -- it does not.  Prediction written before the
    // run: the ratio goes to 1 at integer periods.
    struct Case{double orbits;int points;};
    for(const Case& cs:{
            Case{1.0,4},Case{2.0,4}}){
        const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
        const State start=osculatingPeriapsisState(el,kk,m1,
            m1*(secondMagneticMoment/firstMagneticMoment),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
        acc.reactionModel=
            ChargeRadiationReactionModel::individualLandauLifshitz;
        acc.computeOutwardFlux=true;
        acc.farFieldSampling.directionCount=194;
        acc.farFieldSampling.controlRadius=1.0e6*bohrRadius;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        Vec3 firstJ{},lastJ{},firstA{},lastA{},firstM{},lastM{};
        double firstE=0.0,lastE=0.0; int taken=0;
        const int stepsPerCheckpoint=
            int(200.0*cs.orbits)/cs.points;
        for(int checkpoint=0;checkpoint<cs.points&&ok;++checkpoint){
            for(int i=0;i<stepsPerCheckpoint&&ok;++i)
                ok=engine.advance(st,period/200);
            if(!ok) break;
            const Vec3 J=noetherAngularMomentum(st);
            if(taken==0){ firstJ=J; firstA=st.radiatedAngularMomentum;
                          firstE=st.radiatedEnergy;
                          firstM=st.reactionAngularMomentumMismatch; }
            lastJ=J; lastA=st.radiatedAngularMomentum;
            lastE=st.radiatedEnergy;
            lastM=st.reactionAngularMomentumMismatch; ++taken;
        }
        if(taken<2){ std::printf("%9.2f %6d %8s %7d przerwane\n",
            cs.orbits,194,"nie",cs.points); continue; }
        const double n=taken-1;
        const Vec3 dJ=(lastJ-firstJ)*(1.0/n);
        const Vec3 dA=(lastA-firstA)*(1.0/n);
        const double dE=(lastE-firstE)/n;
        const Vec3 dM=(lastM-firstM)*(1.0/n);
        const double gap=dJ.z+dA.z;
        std::printf("%9.2f %13.5e %13.5e %12.5e %12.5e %9.4f\n",
            cs.orbits,dJ.z,dA.z,gap,dM.z,
            std::abs(gap)>0.0?dM.z/gap:0.0);
    }
}
