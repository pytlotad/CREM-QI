// Audit 282: the third column -- LINEAR MOMENTUM -- and a bound on
// the quadrature error of audit 281.
//
// For a MUTUAL force Newton's third law makes F_1 + F_2 = 0, so the
// sum over both particles is a third-law residual, not a net force.
// Audit 241 found a momentum imbalance and 250b found it collapsing
// with the moments off; here it is resolved per component, which 241
// could not do.  The self reaction is not a mutual force and has no
// reason to cancel.
//
// 281d left the rectangle rule unbounded, and its error could be
// comparable to TEST 4's 0.47 percent residual -- which matters,
// because those numbers are now in README.  Sample count is doubled
// here and the closures recomputed, so the difference bounds it.
//
// Audit 280 header follows.
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
#ifndef SAMPLES
#define SAMPLES 100
#endif
// Okno STALE (pol obiegu) niezaleznie od SAMPLES: krok to
// period/(2*SAMPLES).  Pierwsza wersja mnozyla okno razem z
// liczba probek i nie mierzyla bledu kwadratury wcale.
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
    Vec3 forceFirst[COUNT]{},forceSecond[COUNT]{};
    Vec3 fluxAngular{},fluxMomentum{}; double fluxEnergy=0.0;
    // Conservative bookkeeping, so the power column can close.  The
    // velocity field's work is almost entirely -dU/dt and does not
    // radiate; without separating it the power residual is 55 times the
    // flux (audit 280f).
    double kineticFirst=0.0,kineticLast=0.0;
    double mechanicalFirst=0.0,mechanicalLast=0.0;
    double radiatedFirst=0.0,radiatedLast=0.0;
    double timeFirst=0.0,timeLast=0.0;
    int n=0;
    for(int i=0;i<SAMPLES&&ok;++i){
        ok=engine.advance(st,period/(2.0*SAMPLES));
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
            forceFirst[slot]=forceFirst[slot]+f1;
            forceSecond[slot]=forceSecond[slot]+f2;
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
        const double weight=retardedDipoleSectorWeight(st);
        const double kinetic=kineticEnergy(st.firstVelocity,firstMass)
            +kineticEnergy(st.secondVelocity,secondMass);
        const double mechanical=conservativeParticleEnergy(st,weight);
        if(n==0){ kineticFirst=kinetic; mechanicalFirst=mechanical;
                  radiatedFirst=st.radiatedEnergy; timeFirst=st.time; }
        kineticLast=kinetic; mechanicalLast=mechanical;
        radiatedLast=st.radiatedEnergy; timeLast=st.time;
        ++n;
    }
    if(n<1){ std::printf("brak probek\n"); return 1; }
    const double inv=1.0/n;
    std::printf("# %d probek, wszystko z jednego stanu\n",n);
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
    // ===== ENERGIA: trzy niezalezne domkniecia =====
    const double span=timeLast-timeFirst;
    if(!(span>0.0)) return 0;
    const double dt=span/(n-1);
    double workSum=0.0;
    for(int k=0;k<COUNT;++k) workSum+=power[k]*dt;
    const double deltaKinetic=kineticLast-kineticFirst;
    const double deltaMechanical=mechanicalLast-mechanicalFirst;
    const double radiated=radiatedLast-radiatedFirst;
    const double workRadiative=(power[ACC]+power[SELF])*dt;
    const double workVelocity=power[VEL]*dt;
    const double deltaPotential=deltaMechanical-deltaKinetic;
    std::printf("\n# ENERGIA, okno %.4e s\n",span);
    std::printf("%42s %14.6e\n","suma prac czlonow",workSum);
    std::printf("%42s %14.6e\n","zmiana energii kinetycznej",deltaKinetic);
    std::printf("%42s %14.6e\n","  TEST 1 wyliczenie sil (ma byc 0)",
        workSum-deltaKinetic);
    std::printf("%42s %14.6e\n","zmiana energii mechanicznej",
        deltaMechanical);
    std::printf("%42s %14.6e\n","wypromieniowana energia",radiated);
    std::printf("%42s %14.6e\n","  TEST 2 ksiega energii (ma byc 0)",
        deltaMechanical+radiated);
    std::printf("%42s %14.6e\n","praca czlonow RADIACYJNYCH",
        workRadiative);
    std::printf("%42s %14.6e\n","  TEST 3 radiacyjne + wypromien. (0)",
        workRadiative+radiated);
    std::printf("%42s %14.6e\n","praca czlonu PREDKOSCIOWEGO",
        workVelocity);
    std::printf("%42s %14.6e\n","zmiana energii potencjalnej",
        deltaPotential);
    std::printf("%42s %14.6e\n","  TEST 4 predkosciowe + dU (ma byc 0)",
        workVelocity+deltaPotential);
    // ===== PED: test trzeciej zasady per czlon =====
    std::printf("\n# PED: reszta trzeciej zasady |F1+F2| oraz |F1|\n");
    std::printf("%30s %13s %13s %11s\n","czlon","|F1+F2| [N]",
        "|F1| [N]","reszta/|F1|");
    Vec3 thirdLawSum{};
    for(int k=0;k<COUNT;++k){
        const Vec3 residualForce=(forceFirst[k]+forceSecond[k])*inv;
        const double scaleForce=forceFirst[k].norm()*inv;
        thirdLawSum=thirdLawSum+residualForce;
        std::printf("%30s %13.4e %13.4e %11.3e\n",names[k],
            residualForce.norm(),scaleForce,
            scaleForce>0.0?residualForce.norm()/scaleForce:0.0);
    }
    std::printf("%30s %13.4e\n","SUMA reszt",thirdLawSum.norm());
    std::printf("%30s %13.4e\n","strumien pedu pola",
        fluxMomentum.norm()*inv);
    std::printf("%30s %13.4e\n","  iloraz reszty do strumienia",
        fluxMomentum.norm()>0.0
            ?thirdLawSum.norm()/(fluxMomentum.norm()*inv):0.0);
    std::printf("\n# udzialy w wypromieniowanej energii:\n");
    for(int k=0;k<COUNT;++k)
        std::printf("%30s %12.4e\n",names[k],
            radiated!=0.0?-power[k]*dt/radiated:0.0);
}
