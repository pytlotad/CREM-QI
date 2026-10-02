// Audit 297b: czy srodek masy w para dryfuje ze STALA predkoscia?
// PREDYKCJA (z wykladnika 0,57 -> 0,84 -> 1): v_cm na koncu kolejnych
// orbit stale do kilku procent, a |<F1+F2>_okres| << |F1|.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const OsculatingElements el{-kk/(2.0*a),std::sqrt(kk*a)};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
    const Vec3 m2=Vec3{0,0,1}*secondMagneticMoment;      // para
    const State start=osculatingPeriapsisState(el,kk,m1,m2,
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false; acc.useRetardedExternalForces=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start; bool ok=true;
    const double M=firstMass+secondMass;
    {
        const Vec3 v0=(start.firstVelocity*firstMass+start.secondVelocity*secondMass)/M;
        const Vec3 r0=(start.firstPosition*firstMass+start.secondPosition*secondMass)/M;
        std::printf("# t=0: r_cm=(%.4e,%.4e,%.4e)  v_cm=(%.4e,%.4e,%.4e)  T=%.4e s\n",
            r0.x,r0.y,r0.z,v0.x,v0.y,v0.z,period);
    }
    std::printf("%3s %13s %13s %13s %13s\n","orb","|v_cm| [m/s]","|<2F1>|/<|F1|>","|r_cm| [m]","d|r_cm|/dt");
    double prevR=0.0;
    for(int orbit=1;orbit<=8&&ok;++orbit){
        Vec3 sumF{}; double sumAbs=0.0;
        for(int i=0;i<200&&ok;++i){
            ok=engine.advance(st,period/200); if(!ok) break;
            State loc{st}; loc.firstAcceleration={}; loc.secondAcceleration={};
            const StateHistory h{loc};
            const MutualForces f=chargeDipoleForces(st,h);
            sumF=sumF+f.first+f.second; sumAbs+=f.first.norm();
        }
        const Vec3 vcm=(st.firstVelocity*firstMass+st.secondVelocity*secondMass)/M;
        const Vec3 rcm=(st.firstPosition*firstMass+st.secondPosition*secondMass)/M;
        std::printf("%3d %13.5e %13.5e %13.5e %13.5e   r_cm=(%.3e,%.3e,%.3e)\n",orbit,vcm.norm(),
            sumF.norm()/sumAbs,rcm.norm(),(rcm.norm()-prevR)/period,rcm.x,rcm.y,rcm.z);
        prevR=rcm.norm();
    }
}
