// Audit 256: does the emission DIRECTION move with tolerance?
//
// 255b closed the amplitude channel: the interpolant does not shift
// harmonic amplitudes at all.  254e named the channel where it must
// still show, because there the error is the leading term and not a
// correction: the angle.  The interpolated acceleration is the true
// one rotated forward by one node spacing (254b), first order, and
// the angle enters the angular distribution of the radiation and so
// the photon's recoil direction.
//
// The recoil is the momentum flux through the far-field control
// sphere, electromagneticFieldFluxRates::momentum, built from
// farZoneChargeField, which takes its acceleration from
// historicalCharge.  For a pure E1 source the angular pattern is
// symmetric and the net momentum vanishes, so what is measured here
// is an interference quantity (E1-E2, E1-M1) -- which makes it MORE
// exposed to the interpolant than the power is, not less.
//
// Two controls, not one: tolerance moves the interpolant, direction
// count moves the sphere quadrature.  Without the second column a
// moving direction could be quadrature error alone.  All rows are at
// the same orbital phase, one quarter period, reached by the same
// fixed requested steps.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const double ecc=0.5;
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const Vec3 mz=Vec3{0.0,0.0,1.0};
    std::printf("# a=0.25 a_pair, e=%.2f, momenty WLACZONE (potrzebne dla "
        "interferencji E1-M1), cwierc okresu\n",ecc);
    std::printf("%10s %7s %14s %14s %12s %14s %14s\n",
        "tolerancja","kierunk","energia [W]","|ped| [N]","|p|c/E",
        "kat do odn.","beta");
    std::vector<Vec3> dirs; std::vector<double> tols,counts;
    Vec3 reference{};
    for(double tol:{1.0e-8,1.0e-10,1.0e-12}){
        for(int nd:{50,194,302}){
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
            if(!ok){ std::printf("%10.0e %7d  nieudane\n",tol,nd); continue; }
            FarFieldSampling smp; smp.directionCount=nd;
            const FieldFluxRates r=
                electromagneticFieldFluxRates(st,engine.history(),smp);
            const double pn=r.momentum.norm();
            const Vec3 dir=pn>0.0?r.momentum*(1.0/pn):Vec3{};
            if(reference.squaredNorm()==0.0&&pn>0.0) reference=dir;
            const double ang=(pn>0.0&&reference.squaredNorm()>0.0)
                ?std::acos(std::clamp(dot(dir,reference),-1.0,1.0)):0.0;
            const double beta=(st.firstVelocity-st.secondVelocity).norm()/c;
            std::printf("%10.0e %7d %14.6e %14.6e %12.4e %14.6e %14.6e\n",
                tol,nd,r.energy,pn,pn*c/std::max(r.energy,1e-300),ang,beta);
            dirs.push_back(dir); tols.push_back(tol); counts.push_back(nd);
        }
    }
    // Pairwise angles: tolerance at fixed quadrature, quadrature at
    // fixed tolerance.  This is the split that decides the question.
    std::printf("\n# kat miedzy kierunkami [rad]\n");
    std::printf("# TOLERANCJA przy stalej kwadraturze:\n");
    for(std::size_t j=0;j<3;++j)
        for(std::size_t i=0;i+1<3;++i){
            const std::size_t p=i*3+j,q=(i+1)*3+j;
            if(p<dirs.size()&&q<dirs.size())
                std::printf("#   nd=%3.0f  %.0e -> %.0e : %.6e\n",
                    counts[p],tols[p],tols[q],
                    std::acos(std::clamp(dot(dirs[p],dirs[q]),-1.0,1.0)));
        }
    std::printf("# KWADRATURA przy stalej tolerancji:\n");
    for(std::size_t i=0;i<3;++i)
        for(std::size_t j=0;j+1<3;++j){
            const std::size_t p=i*3+j,q=i*3+j+1;
            if(p<dirs.size()&&q<dirs.size())
                std::printf("#   tol=%.0e  %3.0f -> %3.0f : %.6e\n",
                    tols[p],counts[p],counts[q],
                    std::acos(std::clamp(dot(dirs[p],dirs[q]),-1.0,1.0)));
        }
}
