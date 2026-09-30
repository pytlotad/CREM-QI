// Audit 254: how far is the interpolant's second derivative from the
// dynamical acceleration, and does the gap modulate with orbital
// phase?
//
// 253b established that the acceleration entering the radiated field
// is interpolatedCharge's d2/dt2 of a cubic Hermite in position and
// velocity, never State::firstAcceleration.  farZoneChargeField
// (electrodynamics.hpp:447) builds the 1/R term from exactly that
// object, and its flux is what feeds orbitalRadiatedEnergy and hence
// the photon's firing time.  So the question of "in what field
// configuration the photon is generated" runs through this gap.
//
// Measured at every interior history node, where the integrator's
// own acceleration is stored and is therefore an independent
// reference:
//   err  = |a_interp(t_i) - a_stored(t_i)| / |a_stored(t_i)|
//   jump = |a_interp(t_i+) - a_interp(t_i-)| / |a_stored(t_i)|
// The jump is the cleanest statement: the Hermite second derivative
// is LINEAR in the segment parameter, so it is discontinuous across
// every node, and the radiated power inherits that sawtooth.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const double ecc=std::getenv("ECC")?std::atof(std::getenv("ECC")):0.5;
    std::printf("# a=0.25 a_pair, mimosrod e=%.2f, jeden pelny obieg\n",ecc);
    std::printf("%9s %7s %10s %10s %10s %10s %10s %10s\n",
        "tolerancja","wezly","err med","dfaza/wez","blad dlug","blad kata",
        "dlug/dfaza2","kat/dfaza");
    for(double tol:{1.0e-8,1.0e-10,1.0e-12}){
        const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a*(1.0-ecc*ecc))};
        const double period=osculatingPeriod(el.specificEnergy,k);
        const Vec3 mz=Vec3{0.0,0.0,1.0};
        const State start=osculatingPeriapsisState(el,k,
            mz*firstMagneticMoment,mz*secondMagneticMoment,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start; bool ok=true;
        for(int i=0;i<400&&ok;++i) ok=engine.advance(s,period/400);
        if(!ok){ std::printf("%9.0e  evolution failed\n",tol); continue; }
        const StateHistory& h=engine.history();
        std::vector<double> errs,jumps,spans,mags,angs; double ePeri=0,eApo=0,rPeri=1e30,rApo=0;
        for(std::size_t i=1;i+1<h.size();++i){
            const Vec3 aStored=h[i].firstAcceleration;
            const double an=aStored.norm();
            if(!(an>0.0)) continue;
            const double t=h[i].time;
            const double dl=t-h[i-1].time, dr=h[i+1].time-t;
            if(!(dl>0.0)||!(dr>0.0)) continue;
            const Vec3 aL=historicalCharge(h,s,true,t-0.25*dl).acceleration;
            const Vec3 aR=historicalCharge(h,s,true,t+0.25*dr).acceleration;
            const Vec3 aAt=historicalCharge(h,s,true,t).acceleration;
            const double err=(aAt-aStored).norm()/an;
            const double jump=(aR-aL).norm()/an;
            errs.push_back(err); jumps.push_back(jump);
            spans.push_back(0.5*(dl+dr)/period);
            // Split the error into length and direction.
            mags.push_back(std::abs(aAt.norm()-an)/an);
            angs.push_back(std::acos(std::clamp(
                dot(aAt,aStored)/(aAt.norm()*an),-1.0,1.0)));
            const double r=(h[i].firstPosition-h[i].secondPosition).norm();
            if(r<rPeri){ rPeri=r; ePeri=err; }
            if(r>rApo){ rApo=r; eApo=err; }
        }
        if(errs.empty()){ std::printf("%9.0e  brak wezlow\n",tol); continue; }
        std::sort(errs.begin(),errs.end()); std::sort(jumps.begin(),jumps.end());
        std::sort(spans.begin(),spans.end());
        std::sort(mags.begin(),mags.end()); std::sort(angs.begin(),angs.end());
        const double dphi=2.0*pi*spans[spans.size()/2];
        std::printf("%9.0e %7zu %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
            tol,errs.size(),errs[errs.size()/2],dphi,
            mags[mags.size()/2],angs[angs.size()/2],
            mags[mags.size()/2]/(dphi*dphi),angs[angs.size()/2]/dphi);
    }
}
