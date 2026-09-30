// Audit 259: does the new production moment output agree with the
// probe that motivated it?  257h flagged this as the first step of
// the realization, not of the proposal.
//
// Three checks, all of them falsifiable:
//   a00 == energia   -- free invariant: a00's basis function is 1.
//   momenty: prod vs sonda -- the same six numbers, reconstructed
//                      independently by replicating the loop.
//   c2/c0 wzgledem L -- must come out 0.500, which is what production
//                      ASSUMES today (pdf (3/8)(1+mu^2)).  Here it is
//                      a prediction of the accumulated moments, so a
//                      miss is a wiring error, not a physics result.
// The l=2 part is nQn with Q symmetric traceless; the P2 coefficient
// about an axis is exactly a.Q.a, so c2/c0 = a.Q.a/(a00/4pi).
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    std::printf("%22s %12s %12s %10s\n",
        "konfiguracja","a00/energia","moment max","c2/c0 (L)");
    struct Case{const char* n;double e;double ms;};
    for(const Case& cs:{Case{"kolowa, momenty off",0.0,1.0e-6},
                        Case{"e=0.50, momenty off",0.50,1.0e-6},
                        Case{"e=0.80, momenty off",0.80,1.0e-6},
                        Case{"e=0.50, momenty ON",0.50,1.0}}){
        const OsculatingElements el{-kk/(2.0*a),
            std::sqrt(kk*a*(1.0-cs.e*cs.e))};
        const double period=osculatingPeriod(el.specificEnergy,kk);
        const Vec3 mz=Vec3{0.0,0.0,1.0};
        const State start=osculatingPeriapsisState(el,kk,
            mz*(firstMagneticMoment*cs.ms),mz*(secondMagneticMoment*cs.ms),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<50&&ok;++i) ok=engine.advance(st,0.25*period/50);
        if(!ok){ std::printf("%22s nieudane\n",cs.n); continue; }
        const StateHistory& h=engine.history();
        FarFieldSampling smp;
        smp.directionCount=302;
        smp.controlRadius=1.0e7*bohrRadius;
        smp.radiationFieldOnly=true;
        const FieldFluxRates r=electromagneticFieldFluxRates(st,h,smp);
        // Independent replication of the same six projections.
        const auto quad=sphereQuadratureView(302);
        const Vec3 centre=(st.firstPosition+st.secondPosition)*0.5;
        const double ext=std::max((st.firstPosition-centre).norm(),
                                  (st.secondPosition-centre).norm());
        const double wave=st.time-ext/c;
        const Vec3 comV=(st.firstVelocity*firstMass
            +st.secondVelocity*secondMass)/(firstMass+secondMass);
        double M[6]={0,0,0,0,0,0};
        for(const SphereQuadraturePoint& p:quad){
            const Vec3 n=p.direction, obs=centre+n*smp.controlRadius;
            ElectromagneticField f=farZoneChargeField(
                obs,n,wave,centre,h,st,true,firstCharge,true);
            const ElectromagneticField f2=farZoneChargeField(
                obs,n,wave,centre,h,st,false,secondCharge,true);
            const ElectromagneticField d1=farZoneMagneticDipoleField(
                obs,n,wave,centre,h,st,true,true);
            const ElectromagneticField d2=farZoneMagneticDipoleField(
                obs,n,wave,centre,h,st,false,true);
            f.electric=f.electric+f2.electric+d1.electric+d2.electric;
            f.magnetic=f.magnetic+f2.magnetic+d1.magnetic+d2.magnetic;
            const double pw=dot(cross(f.electric,f.magnetic)*(1.0/mu0),n)
                *smp.controlRadius*smp.controlRadius*p.solidAngleWeight
                *(1.0-dot(n,comV)/c);
            M[0]+=pw; M[1]+=pw*n.x*n.y; M[2]+=pw*n.y*n.z;
            M[3]+=pw*0.5*(3.0*n.z*n.z-1.0); M[4]+=pw*n.x*n.z;
            M[5]+=pw*0.5*(n.x*n.x-n.y*n.y);
        }
        const double prod[6]={r.pattern.a00,r.pattern.a2[0],r.pattern.a2[1],
            r.pattern.a2[2],r.pattern.a2[3],r.pattern.a2[4]};
        double worst=0.0;
        for(int i=0;i<6;++i){
            const double sc=std::max(std::abs(prod[0]),1e-300);
            worst=std::max(worst,std::abs(prod[i]-M[i])/sc);
        }
        // c2/c0 about the pair's angular momentum.
        const Vec3 rel=st.firstPosition-st.secondPosition;
        const Vec3 Lv=cross(rel,st.firstVelocity-st.secondVelocity);
        const Vec3 ax=Lv*(1.0/Lv.norm());
        const double f4pi=1.0/(4.0*pi);
        const double Cxy=prod[1]*15.0*f4pi, Cyz=prod[2]*15.0*f4pi;
        const double Cm0=prod[3]*5.0*f4pi,  Cxz=prod[4]*15.0*f4pi;
        const double Cm2=prod[5]*15.0*f4pi;
        const double Qxx=-0.5*Cm0+0.5*Cm2, Qyy=-0.5*Cm0-0.5*Cm2, Qzz=Cm0;
        const double aQa=Qxx*ax.x*ax.x+Qyy*ax.y*ax.y+Qzz*ax.z*ax.z
            +Cxy*ax.x*ax.y+Cyz*ax.y*ax.z+Cxz*ax.x*ax.z;
        std::printf("%22s %12.4e %12.4e %10.5f\n",cs.n,
            prod[0]/r.energy,worst,aQa/(prod[0]*f4pi));
    }
}
