// Audit 256: why does the radiated MOMENTUM get worse as the
// integrator tolerance tightens, while the energy does not?
// Hypothesis: at tight tolerance the step shrinks, the retention
// window (4 sep/c) then holds more than maximumHistoryNodes=256
// nodes, decimation halves it and leaves UNEVEN segment spans.  An
// uneven grid degrades the Hermite second derivative, hence the
// acceleration field -- and the momentum, being a ~1e-06
// cancellation residue in the angular pattern, is destroyed while
// the energy, which has no cancellation, survives.
// Printed beside the flux: node count and the span unevenness.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A, ecc=0.5;
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const Vec3 mz=Vec3{0.0,0.0,1.0};
    std::printf("%10s %7s %10s %14s %14s %12s %12s\n",
        "tolerancja","wezly","span max/min","energia [W]","|ped| [N]",
        "|p|c/E","kat do 1e-08");
    Vec3 ref{};
    for(double tol:{1.0e-8,3.0e-9,1.0e-9,1.0e-10,1.0e-11,1.0e-12}){
        const State start=osculatingPeriapsisState(el,k,
            mz*firstMagneticMoment,mz*secondMagneticMoment,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<50&&ok;++i) ok=engine.advance(st,0.25*period/50);
        if(!ok){ std::printf("%10.0e  nieudane\n",tol); continue; }
        const StateHistory& h=engine.history();
        double smin=1e300,smax=0.0;
        for(std::size_t i=1;i<h.size();++i){
            const double s=h[i].time-h[i-1].time;
            if(s>0.0){ smin=std::min(smin,s); smax=std::max(smax,s); }
        }
        FarFieldSampling smp;
        smp.directionCount=302;
        smp.controlRadius=1.0e7*bohrRadius;
        smp.radiationFieldOnly=true;
        const FieldFluxRates r=
            electromagneticFieldFluxRates(st,h,smp);
        const double pn=r.momentum.norm();
        const Vec3 dir=pn>0.0?r.momentum*(1.0/pn):Vec3{};
        if(ref.squaredNorm()==0.0&&pn>0.0) ref=dir;
        const double ang=(pn>0.0&&ref.squaredNorm()>0.0)
            ?std::acos(std::clamp(dot(dir,ref),-1.0,1.0)):0.0;
        std::printf("%10.0e %7zu %10.2f %14.6e %14.6e %12.4e %12.4e\n",
            tol,h.size(),smax>0.0?smax/smin:0.0,r.energy,pn,
            pn*c/std::max(r.energy,1e-300),ang);
    }
}
