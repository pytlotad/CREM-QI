// Audit 266: the momentum-consistency test of 257f, and why it cannot
// pass as registered.
//
// 257f pre-registered: <n> E_gamma/c over many photons against the
// independently accumulated radiatedMomentum.  It was registered BEFORE
// 259a decided to carry only even l, and the two are in conflict.
//
// An even-l pattern obeys f(-n) = f(n), so the integral of n f(n) over
// the sphere is EXACTLY zero -- substituting n -> -n flips its sign and
// leaves it equal to itself.  So <n> from the sampled photon directions
// is zero in expectation by construction, while radiatedMomentum is not.
// The test compares zero with a nonzero number.
//
// This is not a defect of the wiring, and it is not new: the PRESCRIBED
// draw, (1+cos^2 theta) about an axis with a uniform azimuth, is equally
// even, so the model's discrete photon recoil has had zero mean all
// along while its continuous field momentum flux has not.
//
// What can be measured is the SIZE of that pre-existing gap, in the one
// measure both sides share: a dimensionless anisotropy.  For the field
// it is |p| c / E.  For the photons it is |<n>|, since each photon
// carries E_gamma/c along n.  This probe puts the two side by side, and
// checks that the photon side really is statistically zero -- which is
// the part that IS a test of the wiring: it says the pattern path adds
// no spurious systematic recoil.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <random>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    std::printf("%10s %14s %14s %14s %12s\n","mimosrod",
        "pole |p|c/E","fotony |<n>|","szum 1/sqrt(N)","luka");
    for(double ecc:{0.0,0.25,0.50,0.80}){
        const OsculatingElements el{-kk/(2.0*a),
            std::sqrt(kk*a*(1.0-ecc*ecc))};
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
        if(!ok){ std::printf("%10.2f nieudane\n",ecc); continue; }
        FarFieldSampling smp;
        smp.directionCount=302;
        smp.controlRadius=1.0e7*bohrRadius;
        smp.radiationFieldOnly=true;
        const FieldFluxRates r=
            electromagneticFieldFluxRates(st,engine.history(),smp);
        const double fieldAnisotropy=
            r.momentum.norm()*c/std::max(r.energy,1.0e-300);
        // The same pattern, sampled through the production helper.
        const Vec3 rel=st.firstPosition-st.secondPosition;
        const Vec3 Lv=cross(rel,st.firstVelocity-st.secondVelocity);
        const Vec3 axis=Lv*(1.0/Lv.norm());
        // Periapsis direction of this orbit, from the eccentricity vector.
        const double strength=kk;
        const double distance=rel.norm();
        const Vec3 vrel=st.firstVelocity-st.secondVelocity;
        const Vec3 eccVector=(rel*(vrel.squaredNorm()-strength/distance)
            -vrel*dot(rel,vrel))*(1.0/strength);
        const Vec3 peri=eccVector.norm()>1.0e-6
            ?eccVector*(1.0/eccVector.norm()):cross(axis,Vec3{1,0,0});
        std::mt19937_64 rng(31337);
        std::uniform_real_distribution<double> u(0.0,1.0);
        const auto dr=[&]{ return u(rng); };
        const int n=4000000;
        Vec3 sum{}; long long got=0;
        for(int i=0;i<n;++i){
            const Vec3 d=drawLabDirectionFromOrbitalPattern(
                r.pattern,axis,peri,dr);
            if(d.squaredNorm()==0.0) continue;
            ++got;
            sum=sum+d;
        }
        const double photonAnisotropy=sum.norm()/got;
        const double noise=1.0/std::sqrt(double(got));
        std::printf("%10.2f %14.4e %14.4e %14.4e %12.4e\n",ecc,
            fieldAnisotropy,photonAnisotropy,noise,fieldAnisotropy);
    }
    std::printf("# \"luka\" = anizotropia pola, bo strona fotonowa jest "
        "zerem z konstrukcji\n");
}
