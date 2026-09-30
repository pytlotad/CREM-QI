// Audit 258: does the m=2 azimuthal structure SURVIVE the emission
// window?  The premise check 257 skipped.
//
// 257c measured the INSTANTANEOUS pattern and found m=2 at 0.2475 of
// the power, phase locked to the separation direction.  But a CREM
// photon is not emitted instantaneously: energy accumulates until a
// threshold, and 229 measured that window at 1 to 3 orbits.  Over a
// full orbit the separation direction sweeps 2pi, so a structure
// locked to it sweeps 4pi in the m=2 harmonic and averages toward
// zero in the LAB frame.  If it averages to zero, production's
// uniform azimuth is CORRECT for whole-orbit accumulation and the
// stage-1 proposal collapses.
//
// Two columns, which is the whole point:
//   LAB      -- harmonic accumulated in a fixed basis.  This is what
//               the emitted photon's direction distribution actually
//               is, and what production must match.
//   OBROTOWY -- harmonic accumulated relative to the instantaneous
//               separation direction.  This is what 257c measured and
//               it cannot average away by construction.
// Reported against window length, as a fraction of accumulated c0.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const double ecc=std::getenv("ECC")?std::atof(std::getenv("ECC")):0.5;
    const int steps=48, nd=194;
    const OsculatingElements el{-kk/(2.0*a),std::sqrt(kk*a*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const Vec3 mz=Vec3{0.0,0.0,1.0};
    const State start=osculatingPeriapsisState(el,kk,
        mz*(firstMagneticMoment*1.0e-6),mz*(secondMagneticMoment*1.0e-6),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-9; acc.maximumDepth=24;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start;
    const auto quad=sphereQuadratureView(nd);
    const double R=1.0e7*bohrRadius;
    // Fixed lab basis, taken from the initial orbit normal.
    const Vec3 rel0=start.firstPosition-start.secondPosition;
    const Vec3 v0=start.firstVelocity-start.secondVelocity;
    const Vec3 L0=cross(rel0,v0);
    const Vec3 axis=L0*(1.0/L0.norm());
    const Vec3 lab1=rel0*(1.0/rel0.norm());
    const Vec3 lab2=cross(axis,lab1);
    double c0=0,labC=0,labS=0,rotC=0,rotS=0;
    std::printf("# a=0.25 a_pair, e=%.2f, momenty off, nd=%d, %d krokow "
        "na obieg\n",ecc,nd,steps);
    std::printf("%10s %13s %13s %13s\n",
        "okno/obieg","m=2 LAB","m=2 OBROTOWY","c0 [J]");
    for(int i=0;i<steps;++i){
        if(!engine.advance(st,period/steps)){
            std::printf("nieudane w kroku %d\n",i); return 1; }
        const StateHistory& h=engine.history();
        const Vec3 rel=st.firstPosition-st.secondPosition;
        const Vec3 e1=rel*(1.0/rel.norm());
        const Vec3 e2=cross(axis,e1);
        const Vec3 centre=(st.firstPosition+st.secondPosition)*0.5;
        const double ext=std::max((st.firstPosition-centre).norm(),
                                  (st.secondPosition-centre).norm());
        const double wave=st.time-ext/c;
        const double dt=period/steps;
        for(const SphereQuadraturePoint& p:quad){
            const Vec3 n=p.direction, obs=centre+n*R;
            ElectromagneticField f=farZoneChargeField(
                obs,n,wave,centre,h,st,true,firstCharge,true);
            const ElectromagneticField f2=farZoneChargeField(
                obs,n,wave,centre,h,st,false,secondCharge,true);
            f.electric=f.electric+f2.electric;
            f.magnetic=f.magnetic+f2.magnetic;
            const double w=p.solidAngleWeight
                *dot(cross(f.electric,f.magnetic)*(1.0/mu0),n)*dt;
            c0+=w;
            const double pl=std::atan2(dot(n,lab2),dot(n,lab1));
            const double pr=std::atan2(dot(n,e2),dot(n,e1));
            labC+=w*std::cos(2.0*pl); labS+=w*std::sin(2.0*pl);
            rotC+=w*std::cos(2.0*pr); rotS+=w*std::sin(2.0*pr);
        }
        const int j=i+1;
        if(j==steps/8||j==steps/4||j==steps/2||j==3*steps/4||j==steps)
            std::printf("%10.3f %13.6e %13.6e %13.6e\n",
                double(j)/steps,
                std::sqrt(labC*labC+labS*labS)/c0,
                std::sqrt(rotC*rotC+rotS*rotS)/c0,c0);
    }
}
