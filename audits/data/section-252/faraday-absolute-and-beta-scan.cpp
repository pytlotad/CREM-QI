// Audit 252: why does the Faraday residual floor at 1.85e-04 when the
// Ampere residual in the SAME evaluation reaches 6.4e-07?
//
// First hypothesis: it is not a separate defect.  In the near zone
// curl E and dB/dt are both suppressed by beta^2 relative to curl B
// and dE/dt/c^2, so the same absolute noise divided by a beta^2
// smaller signal shows up as a relative floor beta^-2 higher.  The
// test is to print the ABSOLUTE residuals in matched units -- the
// Faraday residual divided by c has the units of the Ampere
// residual -- and the signal magnitudes beside them.  If the two
// absolute residuals agree, there is one noise source, not two.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double mscale=std::getenv("NOMOM")?1.0e-6:1.0;
    std::printf("# mscale=%.0e  (residuum Faradaya dzielone przez c, "
        "zeby bylo w jednostkach residuum Ampere'a)\n",mscale);
    std::printf("%8s %9s %11s %11s %11s %11s %11s %11s\n",
        "a/a_pair","beta","|curlB|","|rotB-dE|","|curlE|/c","|rotE+dB|/c",
        "Amp rel","Far rel");
    for(double fa:{1.0,0.5,0.25,0.125,0.0625}){
        const double a=fa*A;
        const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}
            *(firstMagneticMoment*mscale);
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
        if(!ok){ std::printf("%8.4f  evolution failed\n",fa); continue; }
        const StateHistory& h=engine.history();
        const double sep=(s.firstPosition-s.secondPosition).norm();
        const Vec3 centre=(s.firstPosition+s.secondPosition)*0.5;
        const double beta=(s.firstVelocity-s.secondVelocity).norm()/c;
        const auto field=[&](const Vec3& x,double t){
            const ElectromagneticField f1=
                fieldFromOtherParticleAt(x,t,s,h,false);
            const ElectromagneticField f2=
                fieldFromOtherParticleAt(x,t,s,h,true);
            return ElectromagneticField{f1.electric+f2.electric,
                                        f1.magnetic+f2.magnetic};
        };
        const Vec3 x=centre+Vec3{0.3,0.5,0.8}*(sep/0.9899494937);
        const double hs=1.0e-4*sep, ht=hs/c;
        Vec3 curlB{},curlE{};
        for(int ax=0;ax<3;++ax){
            Vec3 d{}; (ax==0?d.x:ax==1?d.y:d.z)=hs;
            const ElectromagneticField p=field(x+d,s.time);
            const ElectromagneticField m=field(x-d,s.time);
            const Vec3 dB=(p.magnetic-m.magnetic)*(1.0/(2.0*hs));
            const Vec3 dE=(p.electric-m.electric)*(1.0/(2.0*hs));
            const int k2=(ax+1)%3,l=(ax+2)%3;
            double* cb[3]={&curlB.x,&curlB.y,&curlB.z};
            double* ce[3]={&curlE.x,&curlE.y,&curlE.z};
            const double dBk=(k2==0?dB.x:k2==1?dB.y:dB.z);
            const double dBl=(l==0?dB.x:l==1?dB.y:dB.z);
            const double dEk=(k2==0?dE.x:k2==1?dE.y:dE.z);
            const double dEl=(l==0?dE.x:l==1?dE.y:dE.z);
            *cb[l]+=dBk; *cb[k2]-=dBl;
            *ce[l]+=dEk; *ce[k2]-=dEl;
        }
        const ElectromagneticField tp=field(x,s.time+ht);
        const ElectromagneticField tm=field(x,s.time-ht);
        const Vec3 dEdt=(tp.electric-tm.electric)*(1.0/(2.0*ht));
        const Vec3 dBdt=(tp.magnetic-tm.magnetic)*(1.0/(2.0*ht));
        const double ampAbs=(curlB-dEdt*(1.0/(c*c))).norm();
        const double farAbs=(curlE+dBdt).norm()/c;
        const double sA=std::max(curlB.norm(),dEdt.norm()/(c*c));
        const double sF=std::max(curlE.norm(),dBdt.norm());
        std::printf("%8.4f %9.3e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
            fa,beta,curlB.norm(),ampAbs,curlE.norm()/c,farAbs,
            ampAbs/sA,farAbs*c/sF);
    }
}
