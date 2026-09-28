// Audit 251: is Ampere-Maxwell satisfied by the PRODUCTION field?
//
// The validation backend carries a Yee solver with an explicit
// Ampere-Maxwell update (maxwell_validation_backend.hpp:438, 1162),
// displacement current included, but that is validation-only.  The
// production path evaluates Lienard-Wiechert fields, which are exact
// solutions of Maxwell's equations by construction -- so the
// question is whether anything in production breaks that: the
// Plummer softening (a softened 1/r^2 is not a vacuum solution), the
// Darwin term, or the magnetization term added by hand to the
// two-pole dipole construction.
//
// Measured in vacuum, away from both charges, where the law reduces
// to curl B = (1/c^2) dE/dt.  Faraday, curl E = -dB/dt, is measured
// beside it as a control: both should hold for a genuine retarded
// solution, and a regulator that breaks one need not break the other.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const double mscale=std::getenv("NOMOM")?1.0e-6:1.0;
    const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}*(firstMagneticMoment*mscale);
    const State start=osculatingPeriapsisState(
        el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-12; acc.maximumDepth=26;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State s=start; bool ok=true;
    for(int i=0;i<200&&ok;++i) ok=engine.advance(s,0.25*period/200);
    if(!ok){ std::printf("evolution failed\n"); return 1; }
    const StateHistory& h=engine.history();
    const double sep=(s.firstPosition-s.secondPosition).norm();
    const Vec3 centre=(s.firstPosition+s.secondPosition)*0.5;
    std::printf("# mscale=%.0e separacja %.4e m = %.4f r*, moment radius %.4f r*\n",
        mscale,sep,sep/comptonBarrierRadius,
        magneticDipoleRadius()/comptonBarrierRadius);
    std::printf("%10s %14s %14s %14s %14s\n",
        "d/sep","Ampere rel","Faraday rel","divB rel","divE rel");
    // Both particles' fields summed: the total vacuum field.
    const auto field=[&](const Vec3& x,double t){
        const ElectromagneticField f1=
            fieldFromOtherParticleAt(x,t,s,h,false);
        const ElectromagneticField f2=
            fieldFromOtherParticleAt(x,t,s,h,true);
        return ElectromagneticField{f1.electric+f2.electric,
                                    f1.magnetic+f2.magnetic};
    };
    for(double f:{0.5,1.0,2.0,5.0,20.0,100.0}){
        const Vec3 x=centre+Vec3{0.3,0.5,0.8}*(f*sep/0.9899494937);
        const double hs=1.0e-4*f*sep;
        const double ht=hs/c;
        Vec3 curlB{},curlE{}; double divB=0.0,divE=0.0,sB=0.0,sE=0.0;
        for(int ax=0;ax<3;++ax){
            Vec3 d{}; (ax==0?d.x:ax==1?d.y:d.z)=hs;
            const ElectromagneticField p=field(x+d,s.time);
            const ElectromagneticField m=field(x-d,s.time);
            const Vec3 dB=(p.magnetic-m.magnetic)*(1.0/(2.0*hs));
            const Vec3 dE=(p.electric-m.electric)*(1.0/(2.0*hs));
            // curl_i = eps_ijk d_j F_k ; accumulate by axis
            const int j=ax,k2=(ax+1)%3,l=(ax+2)%3;
            (void)j;
            double* cb[3]={&curlB.x,&curlB.y,&curlB.z};
            double* ce[3]={&curlE.x,&curlE.y,&curlE.z};
            const double dBk=(k2==0?dB.x:k2==1?dB.y:dB.z);
            const double dBl=(l==0?dB.x:l==1?dB.y:dB.z);
            const double dEk=(k2==0?dE.x:k2==1?dE.y:dE.z);
            const double dEl=(l==0?dE.x:l==1?dE.y:dE.z);
            const double dBj=(j==0?dB.x:j==1?dB.y:dB.z);
            const double dEj=(j==0?dE.x:j==1?dE.y:dE.z);
            divB+=dBj; divE+=dEj;
            sB+=std::abs(dBj); sE+=std::abs(dEj);
            *cb[l]+=dBk; *cb[k2]-=dBl;
            *ce[l]+=dEk; *ce[k2]-=dEl;
        }
        const ElectromagneticField tp=field(x,s.time+ht);
        const ElectromagneticField tm=field(x,s.time-ht);
        const Vec3 dEdt=(tp.electric-tm.electric)*(1.0/(2.0*ht));
        const Vec3 dBdt=(tp.magnetic-tm.magnetic)*(1.0/(2.0*ht));
        const Vec3 amp=curlB-dEdt*(1.0/(c*c));
        const Vec3 far=curlE+dBdt;
        const double sA=std::max(curlB.norm(),dEdt.norm()/(c*c));
        const double sF=std::max(curlE.norm(),dBdt.norm());
        std::printf("%10.1f %14.6e %14.6e %14.6e %14.6e\n",
            f,sA>0?amp.norm()/sA:0.0,sF>0?far.norm()/sF:0.0,
            sB>0?std::abs(divB)/sB:0.0,sE>0?std::abs(divE)/sE:0.0);
    }
}
