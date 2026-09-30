// Audit 273: is 272e right?  Does the stored acceleration carry noise
// that the quintic injects and the cubic ignores?
//
// 272e guessed that the quintic is worse because it reads
// State::firstAcceleration directly while the cubic derives the
// acceleration from positions and velocities.  But 254's own data
// already argues the other way: it found the stored-versus-cubic
// difference UNIFORM across nodes (median equal to maximum to three
// digits) and identified it as a coherent phase lead -- angle error
// equal to the phase per node to 1.000, length error to its square to
// 1.007.  A coherent lead is smooth.  Noise is not.
//
// The discriminator is roughness at the node scale, and it is the
// right one because of an asymmetry between the two interpolants: the
// cubic never reads the stored acceleration at all, so any roughness
// in it is invisible to the cubic, while the quintic reproduces it
// EXACTLY at every node.  So if the stored sequence is rough, the
// quintic's field is rough and 272e stands; if it is as smooth as the
// cubic's, 272e is wrong and the cause is elsewhere.
//
// Roughness of a sequence at an interior node: how far the value sits
// from the straight line through its two neighbours, in units of the
// value.  For a smooth function this is of order |f''| h^2 / (2|f|);
// for noise it is the noise itself.  Position is carried as a control,
// because it must come out very smooth whatever the accelerations do.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A, ecc=0.5;
    const OsculatingElements el{-kk/(2.0*a),std::sqrt(kk*a*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const Vec3 mz=Vec3{0.0,0.0,1.0};
    std::printf("%10s %7s %13s %13s %13s %10s\n","tolerancja","wezly",
        "szorst. ZAPIS","szorst. KUBIK","szorst. POZYC","iloraz");
    for(double tol:{1.0e-8,1.0e-9,1.0e-10,1.0e-11}){
        const State start=osculatingPeriapsisState(el,kk,
            mz*firstMagneticMoment,mz*secondMagneticMoment,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<50&&ok;++i) ok=engine.advance(st,0.25*period/50);
        if(!ok){ std::printf("%10.0e nieudane\n",tol); continue; }
        const StateHistory& h=engine.history();
        const auto roughness=[&](const std::vector<Vec3>& f,
                                 const std::vector<double>& t){
            std::vector<double> r;
            for(std::size_t i=1;i+1<f.size();++i){
                const double span=t[i+1]-t[i-1];
                if(!(span>0.0)) continue;
                const double w=(t[i]-t[i-1])/span;
                const Vec3 line=f[i-1]*(1.0-w)+f[i+1]*w;
                const double n=f[i].norm();
                if(n>0.0) r.push_back((f[i]-line).norm()/n);
            }
            if(r.empty()) return 0.0;
            std::sort(r.begin(),r.end());
            return r[r.size()/2];
        };
        std::vector<Vec3> stored,cubic,position;
        std::vector<double> times;
        for(std::size_t i=0;i<h.size();++i){
            stored.push_back(h[i].firstAcceleration);
            position.push_back(h[i].firstPosition);
            times.push_back(h[i].time);
            cubic.push_back(
                historicalCharge(h,st,true,h[i].time).acceleration);
        }
        const double rs=roughness(stored,times);
        const double rc=roughness(cubic,times);
        const double rp=roughness(position,times);
        std::printf("%10.0e %7zu %13.4e %13.4e %13.4e %10.3f\n",
            tol,h.size(),rs,rc,rp,rc>0.0?rs/rc:0.0);
    }
    std::printf("# iloraz >> 1 potwierdza 272e; iloraz ~ 1 je obala\n");
}
