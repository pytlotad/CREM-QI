// Audit 275: 274 repeated with the magnetic dipole fields in the sum,
// and with a defect of 274's own probe fixed.
//
// TWO CHANGES, and the second is a correction.
// (1) farZoneMagneticDipoleField is folded in, so the total is the one
//     262 and 272 measured.  It does not depend on the charge
//     interpolant, so it enters BOTH columns identically -- which is
//     exactly what isolates the question: how does a changed charge
//     field interfere with an unchanged dipole field.
// (2) 274's probe summed POWER per particle.  Production sums the
//     FIELDS and forms the Poynting vector once
//     (electrodynamics.hpp:487 says why), because the cross terms are
//     the interference -- and the momentum IS an interference
//     quantity.  Summing powers discards it.  So 274c's eight percent
//     was computed without the very terms in question.
//
// Original header of 274 follows.
// Audit 274: the test 273e asked for.  Is the cubic-versus-quintic
// difference common to all directions, or direction-dependent?
//
// One trajectory, evolved once.  The interpolant feeds the force, so
// running twice would compare two different orbits; instead both
// accelerations are evaluated on the SAME history, in the same
// process.  The cubic evaluator is checked against production's
// historicalCharge first -- if that control does not reproduce it,
// nothing below means anything.
//
// The physically meaningful refinement of "common or not": split the
// difference in dP/dOmega into its EVEN and ODD parts under n -> -n.
// The momentum integral sees only the odd part, so if the change were
// purely even it could not move the momentum at all.  Reported beside
// the raw spread.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
namespace {
// Quintic Hermite acceleration on the segment containing `time`, from
// the same node data the production interpolant uses.
Vec3 quinticAcceleration(const StateHistory& h,double time,bool first){
    if(h.size()<2) return {};
    std::size_t i=1;
    while(i+1<h.size()&&h[i].time<time) ++i;
    const State& o=h[i-1]; const State& n=h[i];
    const double span=n.time-o.time;
    if(!(span>0.0)) return {};
    const double s=std::clamp((time-o.time)/span,0.0,1.0);
    const Vec3 p0=first?o.firstPosition:o.secondPosition;
    const Vec3 p1=first?n.firstPosition:n.secondPosition;
    const Vec3 v0=first?o.firstVelocity:o.secondVelocity;
    const Vec3 v1=first?n.firstVelocity:n.secondVelocity;
    const Vec3 a0=first?o.firstAcceleration:o.secondAcceleration;
    const Vec3 a1=first?n.firstAcceleration:n.secondAcceleration;
    const double s2=s*s,s3=s2*s;
    const double d0=-60.0*s+180.0*s2-120.0*s3;
    const double d1=-36.0*s+96.0*s2-60.0*s3;
    const double d2=1.0-9.0*s+18.0*s2-10.0*s3;
    const double d3=60.0*s-180.0*s2+120.0*s3;
    const double d4=-24.0*s+84.0*s2-60.0*s3;
    const double d5=3.0*s-12.0*s2+10.0*s3;
    return (p0*d0+v0*(span*d1)+a0*(span*span*d2)
           +p1*d3+v1*(span*d4)+a1*(span*span*d5))/(span*span);
}
}
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A, ecc=0.5;
    const OsculatingElements el{-kk/(2.0*a),std::sqrt(kk*a*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const Vec3 mz=Vec3{0.0,0.0,1.0};
    const State start=osculatingPeriapsisState(el,kk,
        mz*firstMagneticMoment,mz*secondMagneticMoment,
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-10; acc.maximumDepth=26;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start; bool ok=true;
    for(int i=0;i<50&&ok;++i) ok=engine.advance(st,0.25*period/50);
    if(!ok){ std::printf("nieudane\n"); return 1; }
    const StateHistory& h=engine.history();
    // Control: does the local cubic path reproduce production?
    double worstControl=0.0;
    for(std::size_t i=1;i+1<h.size();++i){
        const Vec3 prod=historicalCharge(h,st,true,h[i].time).acceleration;
        const Vec3 quin=quinticAcceleration(h,h[i].time,true);
        const Vec3 stored=h[i].firstAcceleration;
        worstControl=std::max(worstControl,
            (quin-stored).norm()/std::max(stored.norm(),1e-300));
        (void)prod;
    }
    std::printf("# kontrola: kwintyczny w wezlach vs zapisany, max %.3e "
        "(ma byc ~0)\n",worstControl);
    const auto quad=sphereQuadratureView(302);
    const Vec3 centre=(st.firstPosition+st.secondPosition)*0.5;
    const double ext=std::max((st.firstPosition-centre).norm(),
                              (st.secondPosition-centre).norm());
    const double wave=st.time-ext/c;
    const double R=1.0e7*bohrRadius;
    std::vector<double> pc,pq; std::vector<Vec3> dirs;
    for(const SphereQuadraturePoint& p:quad){
        const Vec3 n=p.direction, obs=centre+n*R;
        Vec3 electricCubic{},magneticCubic{};
        Vec3 electricQuintic{},magneticQuintic{};
        for(int which=0;which<2;++which){
            const bool firstSource=(which==0);
            const double charge=firstSource?firstCharge:secondCharge;
            double emit=wave;
            ChargeKinematics src=historicalCharge(h,st,firstSource,emit);
            for(int it=0;it<5;++it){
                emit=wave+dot(n,src.position-centre)/c;
                src=historicalCharge(h,st,firstSource,emit);
            }
            const Vec3 disp=obs-src.position;
            const double dist=disp.norm();
            const Vec3 dir=disp*(1.0/dist);
            const Vec3 beta=src.velocity*(1.0/c);
            const double kap=std::max(1.0e-12,1.0-dot(dir,beta));
            const double pre=coulomb*charge;
            const auto radiate=[&](const Vec3& accel){
                return cross(dir,cross(dir-beta,accel))
                    *(pre/(c*c*kap*kap*kap*dist));
            };
            const Vec3 ec=radiate(src.acceleration);
            const Vec3 eq=radiate(quinticAcceleration(h,emit,firstSource));
            electricCubic=electricCubic+ec;
            electricQuintic=electricQuintic+eq;
            magneticCubic=magneticCubic+cross(dir,ec)*(1.0/c);
            magneticQuintic=magneticQuintic+cross(dir,eq)*(1.0/c);
            // Dipole field of this source: independent of the charge
            // interpolant, so the SAME object goes into both columns.
            const ElectromagneticField dipole=
                farZoneMagneticDipoleField(obs,n,wave,centre,h,st,
                                           firstSource,true);
            electricCubic=electricCubic+dipole.electric;
            electricQuintic=electricQuintic+dipole.electric;
            magneticCubic=magneticCubic+dipole.magnetic;
            magneticQuintic=magneticQuintic+dipole.magnetic;
        }
        const double powerCubic=
            dot(cross(electricCubic,magneticCubic)*(1.0/mu0),n);
        const double powerQuintic=
            dot(cross(electricQuintic,magneticQuintic)*(1.0/mu0),n);
        pc.push_back(powerCubic); pq.push_back(powerQuintic);
        dirs.push_back(n);
    }
    // Per-direction relative difference.
    std::vector<double> rel;
    for(std::size_t i=0;i<pc.size();++i)
        if(std::abs(pq[i])>0.0) rel.push_back((pc[i]-pq[i])/pq[i]);
    double mean=0.0; for(double r:rel) mean+=r; mean/=rel.size();
    double var=0.0; for(double r:rel) var+=(r-mean)*(r-mean);
    const double sd=std::sqrt(var/(rel.size()-1));
    std::printf("roznica wzgledna na kierunek: srednia %+.4e  "
        "rozrzut %.4e  iloraz |sr|/rozrzut %.3f\n",mean,sd,
        sd>0.0?std::abs(mean)/sd:0.0);
    // Even/odd split of the difference under n -> -n, and the momentum
    // each part carries.
    Vec3 momentumCubic{},momentumQuintic{};
    double totalCubic=0.0,totalQuintic=0.0;
    for(std::size_t i=0;i<pc.size();++i){
        const double w=quad[i].solidAngleWeight*R*R;
        totalCubic+=pc[i]*w; totalQuintic+=pq[i]*w;
        momentumCubic=momentumCubic+dirs[i]*(pc[i]*w/c);
        momentumQuintic=momentumQuintic+dirs[i]*(pq[i]*w/c);
    }
    std::printf("moc:  kubik %.6e  kwintyk %.6e  wzgl. roznica %.3e\n",
        totalCubic,totalQuintic,
        std::abs(totalCubic-totalQuintic)/std::abs(totalQuintic));
    std::printf("|ped|c/E: kubik %.4e  kwintyk %.4e\n",
        momentumCubic.norm()*c/std::abs(totalCubic),
        momentumQuintic.norm()*c/std::abs(totalQuintic));
    std::printf("# srednia >> rozrzut = roznica wspolna dla kierunkow\n");
}
