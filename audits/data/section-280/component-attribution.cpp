// Audit 280: unambiguous attribution.  Force, torque and power for
// EVERY constituent, against the field flux, all from one step.
//
// 279's "mutual torque -0.5368" was four physically different things
// in one number.  retardedExternalForces builds chargeCharge as a
// SINGLE Lorentz force from the sum of the charge field and the
// retarded dipole field, and the charge field is itself the sum of a
// VELOCITY term (1/R^2, conservative, near-field) and an
// ACCELERATION term (1/R, radiative).  covariantDipoleGradientForce
// is added separately.  Only the acceleration term has any business
// being set against the radiated flux; if the velocity term carries
// most of that -0.5368, the 53.68 + 50.00 = 103.68 arithmetic was
// comparing a near-field exchange with a radiative one.
//
// The field split replicates lienardWiechertField term by term, the
// same replication audit 259b checked against production to 3e-16.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
namespace {
struct Split { ElectromagneticField velocity, acceleration; };
// Velocity and acceleration parts of the partner's Lienard-Wiechert
// field at `observer`, separated at the source.
Split splitLienardWiechert(const Vec3& observer,double time,
        const StateHistory& h,const State& s,bool sourceIsFirst,
        double sourceCharge){
    ChargeKinematics src=historicalCharge(h,s,sourceIsFirst,time);
    double tr=time-(observer-src.position).norm()/c;
    for(int it=0;it<16;++it){
        src=historicalCharge(h,s,sourceIsFirst,tr);
        const Vec3 d=observer-src.position;
        const double dn=d.norm();
        const Vec3 dir=dn>0?d*(1.0/dn):Vec3{};
        const double res=tr+dn/c-time;
        tr-=res/std::max(1.0e-8,1.0-dot(dir,src.velocity*(1.0/c)));
    }
    src=historicalCharge(h,s,sourceIsFirst,tr);
    const Vec3 disp=observer-src.position;
    const double dist=disp.norm();
    if(!(dist>0.0)) return {};
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
    const double pre=coulomb*sourceCharge*ps;
    const Vec3 ev=(dir-bv)*((1.0-b2)/(kap*kap*kap*fd*fd))*pre;
    const Vec3 ea=cross(dir,cross(dir-bv,src.acceleration))
        *(1.0/(c*c*kap*kap*kap*fd))*pre;
    return {{ev,cross(dir,ev)*(1.0/c)},{ea,cross(dir,ea)*(1.0/c)}};
}
}
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const OsculatingElements el{-kk/(2.0*A),std::sqrt(kk*A)};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
    const State start=osculatingPeriapsisState(el,kk,m1,
        m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    FarFieldSampling sampling;
    sampling.directionCount=194;
    sampling.controlRadius=1.0e6*bohrRadius;
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
    acc.reactionModel=
        ChargeRadiationReactionModel::individualLandauLifshitz;
    acc.computeOutwardFlux=true;
    acc.farFieldSampling=sampling;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start; bool ok=true;
    for(int i=0;i<200&&ok;++i) ok=engine.advance(st,period/200);
    if(!ok){ std::printf("nieudane\n"); return 1; }
    enum { VEL,ACC,DIP,GRAD,SELF,COUNT };
    const char* names[COUNT]={
        "pole predkosciowe ladunku","pole przyspieszeniowe ladunku",
        "pole dipolowe (Lorentz)","gradient dipolowy (tensor)",
        "reakcja wlasna (LL)"};
    Vec3 torque[COUNT]{}; Vec3 force[COUNT]{}; double power[COUNT]{};
    Vec3 fluxAngular{},fluxMomentum{}; double fluxEnergy=0.0;
    int n=0;
    for(int i=0;i<100&&ok;++i){
        ok=engine.advance(st,period/200);
        if(!ok) break;
        const StateHistory& h=engine.history();
        const Split atFirst=splitLienardWiechert(
            st.firstPosition,st.time,h,st,false,secondCharge);
        const Split atSecond=splitLienardWiechert(
            st.secondPosition,st.time,h,st,true,firstCharge);
        const ElectromagneticField dipAtFirst=
            retardedMagneticDipoleField(st.firstPosition,st.time,h,st,false);
        const ElectromagneticField dipAtSecond=
            retardedMagneticDipoleField(st.secondPosition,st.time,h,st,true);
        const MutualForces mutual=retardedExternalForces(st,h);
        const ParticleMultipoleRadiation radiation=
            particleMultipoleRadiation(st,mutual,h,true,
                ChargeRadiationReactionModel::individualLandauLifshitz,
                true,sampling);
        const auto add=[&](int slot,const Vec3& f1,const Vec3& f2){
            torque[slot]=torque[slot]+cross(st.firstPosition,f1)
                                     +cross(st.secondPosition,f2);
            force[slot]=force[slot]+f1+f2;
            power[slot]+=dot(f1,st.firstVelocity)+dot(f2,st.secondVelocity);
        };
        add(VEL,lorentzForce(firstCharge,st.firstVelocity,atFirst.velocity),
                lorentzForce(secondCharge,st.secondVelocity,
                             atSecond.velocity));
        add(ACC,lorentzForce(firstCharge,st.firstVelocity,
                             atFirst.acceleration),
                lorentzForce(secondCharge,st.secondVelocity,
                             atSecond.acceleration));
        add(DIP,lorentzForce(firstCharge,st.firstVelocity,dipAtFirst),
                lorentzForce(secondCharge,st.secondVelocity,dipAtSecond));
        add(GRAD,covariantDipoleGradientForce(st,h,true),
                 covariantDipoleGradientForce(st,h,false));
        add(SELF,radiation.chargeReaction.first,
                 radiation.chargeReaction.second);
        fluxAngular=fluxAngular+radiation.outwardFlux.angularMomentum;
        fluxMomentum=fluxMomentum+radiation.outwardFlux.momentum;
        fluxEnergy+=radiation.outwardFlux.energy;
        ++n;
    }
    if(n<1){ std::printf("brak probek\n"); return 1; }
    const double inv=1.0/n;
    std::printf("# pol obiegu, %d probek, wszystko z jednego stanu\n",n);
    std::printf("%30s %12s %12s %12s\n","czlon","moment z","moc [W]",
        "|sila| [N]");
    Vec3 torqueSum{},forceSum{}; double powerSum=0.0;
    for(int k=0;k<COUNT;++k){
        std::printf("%30s %12.4e %12.4e %12.4e\n",names[k],
            torque[k].z*inv,power[k]*inv,force[k].norm()*inv);
        torqueSum=torqueSum+torque[k]; forceSum=forceSum+force[k];
        powerSum+=power[k];
    }
    std::printf("%30s %12.4e %12.4e %12.4e\n","SUMA czlonow",
        torqueSum.z*inv,powerSum*inv,forceSum.norm()*inv);
    std::printf("%30s %12.4e %12.4e %12.4e\n","strumien pola (ze znakiem -)",
        -fluxAngular.z*inv,-fluxEnergy*inv,fluxMomentum.norm()*inv);
    std::printf("%30s %12.4e %12.4e\n","RESIDUUM suma+strumien",
        (torqueSum.z+fluxAngular.z)*inv,(powerSum+fluxEnergy)*inv);
    std::printf("\n# udzialy w strumieniu momentu pedu:\n");
    for(int k=0;k<COUNT;++k)
        std::printf("%30s %9.4f\n",names[k],
            fluxAngular.z!=0.0?torque[k].z/fluxAngular.z:0.0);
}
