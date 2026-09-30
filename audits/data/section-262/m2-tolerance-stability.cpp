// Audit 262: is the EVEN-l m=2 moment stable against integrator
// tolerance past the point where the history grid saturates?
//
// This is the risk the whole design rests on and it was never
// measured.  256d found the ODD-l channel -- the net momentum --
// capped by the history grid: node count saturates below tolerance
// 1e-10, and past that point tightening the integrator makes the
// momentum WORSE, collapsing fivefold and swinging the direction
// through 0.95 rad, while the energy stays put to six digits.  Odd l
// was excluded from the pattern for exactly that reason (state.hpp).
// But even-l m=2 was ASSUMED safe because it is about a hundred
// times larger.  Assumed, not checked.
//
// So both channels go in one table at the same tolerances: the
// momentum as the known-bad control, the energy and c2/c0 as the
// known-good controls, and the m=2 amplitude and phase as the
// quantities in question.  If m=2 tracks the momentum, stage 1 is
// dead.  If it tracks the energy, it stands.
//
// Phase matters as much as amplitude here: a wandering m=2 phase
// would rotate the emitted azimuth just as a wandering momentum
// direction rotates the recoil.  It is measured against the ORBITAL
// frame, so it is relative to the periapsis.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const double ecc=std::getenv("ECC")?std::atof(std::getenv("ECC")):0.5;
    const double ms=std::getenv("MOM")?1.0:1.0e-6;
    const OsculatingElements el{-kk/(2.0*a),std::sqrt(kk*a*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const Vec3 mz=Vec3{0.0,0.0,1.0};
    std::printf("# a=0.25 a_pair, e=%.2f, momenty %s, cwierc okresu, nd=302\n",
        ecc,ms>0.5?"ON":"off");
    std::printf("%8s %6s %6s %13s %9s %11s %10s %13s %9s\n",
        "tol","wezly","span","energia [W]","c2/c0","m2 ampl","m2 faza",
        "|ped| [N]","kat ped");
    Vec3 refMomentum{};
    for(double tol:{1.0e-8,3.0e-9,1.0e-9,1.0e-10,1.0e-11,1.0e-12}){
        const State start=osculatingPeriapsisState(el,kk,
            mz*(firstMagneticMoment*ms),mz*(secondMagneticMoment*ms),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<50&&ok;++i) ok=engine.advance(st,0.25*period/50);
        if(!ok){ std::printf("%8.0e  nieudane\n",tol); continue; }
        const StateHistory& h=engine.history();
        double smin=1e300,smax=0.0;
        for(std::size_t i=1;i<h.size();++i){
            const double sp=h[i].time-h[i-1].time;
            if(sp>0.0){ smin=std::min(smin,sp); smax=std::max(smax,sp); }
        }
        FarFieldSampling smp;
        smp.directionCount=302;
        smp.controlRadius=1.0e7*bohrRadius;
        smp.radiationFieldOnly=true;
        const FieldFluxRates r=electromagneticFieldFluxRates(st,h,smp);
        const AngularPatternMoments& m=r.pattern;
        const double c2=5.0*m.a2[2]/m.a00;
        const double m2=15.0*std::sqrt(m.a2[0]*m.a2[0]+m.a2[4]*m.a2[4])
            /m.a00;
        const double m2phase=0.5*std::atan2(m.a2[0],m.a2[4]);
        const double pn=r.momentum.norm();
        const Vec3 pdir=pn>0.0?r.momentum*(1.0/pn):Vec3{};
        if(refMomentum.squaredNorm()==0.0&&pn>0.0) refMomentum=pdir;
        const double pang=(pn>0.0&&refMomentum.squaredNorm()>0.0)
            ?std::acos(std::clamp(dot(pdir,refMomentum),-1.0,1.0)):0.0;
        std::printf("%8.0e %6zu %6.2f %13.6e %9.5f %11.6f %10.6f %13.4e %9.2e\n",
            tol,h.size(),smax>0.0?smax/smin:0.0,r.energy,c2,m2,m2phase,
            pn,pang);
    }
}
