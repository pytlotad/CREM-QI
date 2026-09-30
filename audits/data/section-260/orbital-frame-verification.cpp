// Audit 260 verification: the orbital frame, and the accumulation test
// 258d asked for -- run through the PRODUCTION path this time, with
// computeOutwardFlux on, reading State::radiatedPattern after a whole
// orbit rather than reconstructing anything.
//
// Four checks:
//   os perycentrum  -- at a periapsis start, dot(separation, e1) must be
//                      1: the eccentricity vector must point at the
//                      periapsis and not somewhere else.
//   5 a2[2]/a00     -- the P2 coefficient about the angular momentum,
//                      which must stay 0.500 whatever the basis does to
//                      the azimuth.  Frame-change invariant.
//   m2 chwilowe     -- azimuthal m=2 of the instantaneous pattern.
//   m2 po obiegu    -- the same after a full orbit of trapezoidal
//                      accumulation.  258a measured the lab-frame
//                      version collapsing by about 45x; the orbital
//                      frame must do the same over ONE orbit, because
//                      over one orbit the two frames barely differ.
//                      A frame bug would show up as no collapse.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    std::printf("%22s %11s %11s %11s %11s %7s\n","konfiguracja",
        "os peri","5a2[2]/a00","m2 chwil.","m2 po obieg","spadek");
    struct Case{const char* n;double e;double ms;};
    for(const Case& cs:{Case{"kolowa, momenty off",0.0,1.0e-6},
                        Case{"e=0.50, momenty off",0.50,1.0e-6},
                        Case{"e=0.50, momenty ON",0.50,1.0}}){
        const OsculatingElements el{-kk/(2.0*a),
            std::sqrt(kk*a*(1.0-cs.e*cs.e))};
        const double period=osculatingPeriod(el.specificEnergy,kk);
        const Vec3 mz=Vec3{0.0,0.0,1.0};
        const State start=osculatingPeriapsisState(el,kk,
            mz*(firstMagneticMoment*cs.ms),mz*(secondMagneticMoment*cs.ms),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        // Frame at the periapsis start: e1 must be the separation.
        const OrbitalFrame f0=pairOrbitalFrame(start);
        const Vec3 sep0=start.firstPosition-start.secondPosition;
        const double axisCheck=f0.valid
            ?dot(sep0*(1.0/sep0.norm()),f0.first):0.0;
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-9; acc.maximumDepth=24;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=true;
        acc.farFieldSampling.directionCount=194;
        acc.farFieldSampling.controlRadius=1.0e7*bohrRadius;
        acc.farFieldSampling.radiationFieldOnly=true;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        const int steps=48;
        // Instantaneous pattern one step in, for the "chwilowe" column.
        double instM2=0.0,instP2=0.0;
        for(int i=0;i<steps&&ok;++i){
            ok=engine.advance(st,period/steps);
            if(i==0){
                const FieldFluxRates r=electromagneticFieldFluxRates(
                    st,engine.history(),acc.farFieldSampling);
                instM2=15.0*std::sqrt(r.pattern.a2[0]*r.pattern.a2[0]
                    +r.pattern.a2[4]*r.pattern.a2[4])/r.pattern.a00;
                instP2=5.0*r.pattern.a2[2]/r.pattern.a00;
            }
        }
        if(!ok){ std::printf("%22s nieudane\n",cs.n); continue; }
        const AngularPatternMoments& m=st.radiatedPattern;
        const double accM2=m.a00>0.0
            ?15.0*std::sqrt(m.a2[0]*m.a2[0]+m.a2[4]*m.a2[4])/m.a00:0.0;
        std::printf("%22s %11.6f %11.5f %11.4e %11.4e %7.1f\n",cs.n,
            axisCheck,instP2,instM2,accM2,
            accM2>0.0?instM2/accM2:0.0);
    }
}
