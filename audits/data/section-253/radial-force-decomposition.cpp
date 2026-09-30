// Audit 253: does any acceleration-related term shift the force
// balance in positronium from the Compton radius out to the Bohr
// radius?  The claim under test is that a changing magnetic field
// induces an electric repulsion between the charges that the model
// might have omitted.
//
// The acceleration field is NOT omitted -- retarded_charge_kinematics
// .hpp:505 carries the full 1/R term with its kappa^3 factors, and B
// = n x E / c exactly, so every Faraday-induced contribution is in
// the retarded solution by construction.  What is measured here is
// its SIZE: the radial force on one particle decomposed into the
// velocity (Coulomb) field, the acceleration field, the magnetic
// v x B term and the moment sector, each as a fraction of the
// Coulomb magnitude, on circular orbits from the Bohr radius down
// through the Compton barrier.
//
// The acceleration field CANNOT be isolated by zeroing the stored
// accelerations: interpolatedCharge (retarded_charge_kinematics.hpp
// :135) takes the retarded acceleration as the second derivative of
// a cubic Hermite interpolant of position and velocity and never
// reads State::firstAcceleration.  A first attempt that way returned
// a difference of exactly 0.000e+00, which is what a no-op control
// looks like.  The split below reproduces the model's own formula
// term by term at the same converged retarded point, and the
// "suma-model" column checks that the two terms add back up to the
// model's own lienardWiechertField call.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double rStar=comptonBarrierRadius;
    std::printf("# a_pair=%.6e m, r*=%.6e m, r*/a_pair=%.4e\n",
        A,rStar,rStar/A);
    std::printf("# znak dodatni = ODPYCHANIE.  Wszystko w jednostkach "
        "|F_Coulomb| na danym promieniu.\n");
    std::printf("%10s %9s %13s %13s %13s %13s %13s %13s %13s %9s\n",
        "r/a_pair","beta","predkosciowe","przyspiesz.","magnetyczne",
        "m->ladunek","m->m","Darwin","suma/Coul","rozjazd");
    for(double fa:{1.0,0.1,1.0e-2,3.0e-3,1.827e-3,1.0e-3}){
        const double a=fa*A;
        const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        const Vec3 mz=Vec3{0.0,0.0,1.0};
        const State start=osculatingPeriapsisState(el,k,
            mz*firstMagneticMoment,mz*secondMagneticMoment,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start; bool ok=true;
        for(int i=0;i<80&&ok;++i) ok=engine.advance(s,0.25*period/80);
        if(!ok){ std::printf("%10.4g  evolution failed\n",fa); continue; }
        const StateHistory& h=engine.history();
        const Vec3 sep=s.firstPosition-s.secondPosition;
        const double r=sep.norm();
        const Vec3 rhat=sep*(1.0/r);
        const double beta=(s.firstVelocity-s.secondVelocity).norm()/c;
        const double coulombScale=std::abs(coulomb*firstCharge*secondCharge)
            /(r*r);
        // Full charge field at particle 1 from particle 2.
        const ElectromagneticField full=lienardWiechertField(
            s.firstPosition,s.time,h,s,false,secondCharge);
        // Replicate lienardWiechertField term by term at the same root.
        ChargeKinematics src=historicalCharge(h,s,false,s.time);
        double tr=s.time-(s.firstPosition-src.position).norm()/c;
        for(int it=0;it<16;++it){
            src=historicalCharge(h,s,false,tr);
            const Vec3 d=s.firstPosition-src.position;
            const double dn=d.norm();
            const Vec3 dir=dn>0?d*(1.0/dn):Vec3{};
            const double res=tr+dn/c-s.time;
            tr-=res/std::max(1.0e-8,1.0-dot(dir,src.velocity*(1.0/c)));
        }
        src=historicalCharge(h,s,false,tr);
        const Vec3 disp=s.firstPosition-src.position;
        const double dist=disp.norm();
        const Vec3 dir=disp*(1.0/dist);
        const double fd=std::max(dist,nuclearCutoff);
        const Vec3 bv=src.velocity*(1.0/c);
        const double b2=bv.squaredNorm();
        const double kap=std::max(1.0e-8,1.0-dot(dir,bv));
        double ps=1.0;
        if(const double fl=separationFloor(); fl>0.0){
            const double rd=fd*kap/std::sqrt(std::max(1.0-b2,1.0e-300));
            const double ra=rd/std::sqrt(rd*rd+fl*fl); ps=ra*ra*ra;
        }
        const double pref=coulomb*secondCharge*ps;
        const Vec3 eVelT=(dir-bv)*((1.0-b2)/(kap*kap*kap*fd*fd))*pref;
        const Vec3 eAcc=cross(dir,cross(dir-bv,src.acceleration))
            *(1.0/(c*c*kap*kap*kap*fd))*pref;
        const Vec3 bVelT=cross(dir,eVelT)*(1.0/c);
        const Vec3 bAcc=cross(dir,eAcc)*(1.0/c);
        const ElectromagneticField novel{eVelT,bVelT};
        const double closure=(full.electric-eVelT-eAcc).norm()
            /std::max(full.electric.norm(),1.0e-300);
        // Moment sector at the same point.
        const ElectromagneticField mom=retardedMagneticDipoleField(
            s.firstPosition,s.time,h,s,false,1.0e-5);
        const double fVel=dot(novel.electric*firstCharge,rhat);
        const double fAcc=dot(eAcc*firstCharge,rhat);
        const double fMag=dot(cross(s.firstVelocity,
            novel.magnetic+bAcc)*firstCharge,rhat);
        const double fMom=dot((mom.electric
            +cross(s.firstVelocity,mom.magnetic))*firstCharge,rhat);
        // The moment-ON-moment force: this is the term that sets r*.
        const double fDD=dot(pairDipoleForce(sep,
            s.firstDipole,s.secondDipole),rhat);
        const double fDar=dot(darwinForceOnFirst(s.firstVelocity,
            s.secondVelocity,s.secondAcceleration,sep,
            firstCharge*secondCharge),rhat);
        std::printf("%10.4g %9.3e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %9.1e\n",
            fa,beta,fVel/coulombScale,fAcc/coulombScale,
            fMag/coulombScale,fMom/coulombScale,fDD/coulombScale,
            fDar/coulombScale,
            (fVel+fAcc+fMag+fMom+fDD+fDar)/coulombScale,closure);
    }
}
