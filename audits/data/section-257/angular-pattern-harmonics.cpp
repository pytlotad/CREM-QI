// Audit 257 feasibility: what does the model's OWN computed angular
// power pattern look like, and how far is it from the (1+cos^2 theta)
// that production assumes?
//
// 256a established that the production photon direction is drawn
// from CDF(mu)=(mu^3+3mu+4)/8, i.e. pdf=(3/8)(1+mu^2): the rotating
// electric dipole pattern about the orbital normal.  In Legendre
// coefficients that is exactly c2/c0 = 1/2 and every other c_l = 0.
//
// The model already evaluates dP/dOmega at each Lebedev node inside
// electromagneticFieldFluxRates and throws the per-direction values
// away, keeping only the sum.  This probe reconstructs them the same
// way that function does and projects onto Legendre polynomials
// about the angular-momentum axis.
//
// FIRST VERSION WAS THE WRONG PROJECTION.  Legendre-in-mu about L
// averages over azimuth, and the azimuth average of an instantaneous
// dipole pattern about an IN-PLANE acceleration is identically
// (1+mu^2)/2 for any planar orbit, eccentric or not.  That version
// returned c2/c0 = 0.5 for every case including e=0.80 -- it was
// measuring that the acceleration lies in the orbital plane, which is
// trivially true, and destroying exactly the structure in question.
//
// The instantaneous pattern is NOT azimuthally symmetric about L: it
// is a dipole about the acceleration direction, which rotates within
// the plane.  So this version measures the AZIMUTHAL harmonics in the
// orbital plane, with the azimuth origin on the instantaneous
// separation direction.  For a pure dipole about an in-plane axis,
// sin^2(Theta) = 1 - (1-mu^2)(1+cos 2phi)/2, so the m=2 amplitude is
// 2pi/3 against a total of 8pi/3 -- a predicted 0.25 of the total.
// Production draws the azimuth UNIFORMLY (crem_collapse.hpp:4797),
// which discards that entire structure.
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
    std::printf("%22s %9s %9s %9s %9s %9s\n",
        "konfiguracja","c2/c0","m=1","m=2","faza m2","m=3");
    struct Case{const char* name;double ecc;double ms;};
    for(const Case& cs:{Case{"kolowa, momenty off",0.0,1.0e-6},
                        Case{"kolowa, momenty ON",0.0,1.0},
                        Case{"e=0.25, momenty off",0.25,1.0e-6},
                        Case{"e=0.50, momenty off",0.50,1.0e-6},
                        Case{"e=0.50, momenty ON",0.50,1.0},
                        Case{"e=0.80, momenty off",0.80,1.0e-6}}){
        const OsculatingElements el{-kk/(2.0*a),
            std::sqrt(kk*a*(1.0-cs.ecc*cs.ecc))};
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
        if(!ok){ std::printf("%22s nieudane\n",cs.name); continue; }
        const StateHistory& h=engine.history();
        // Angular-momentum axis of the pair, the axis production uses.
        const Vec3 rel=st.firstPosition-st.secondPosition;
        const Vec3 vrel=st.firstVelocity-st.secondVelocity;
        const Vec3 Lv=cross(rel,vrel);
        const Vec3 axis=Lv*(1.0/Lv.norm());
        const auto quad=sphereQuadratureView(302);
        const double R=1.0e7*bohrRadius;
        const Vec3 centre=(st.firstPosition+st.secondPosition)*0.5;
        const double srcExtent=std::max((st.firstPosition-centre).norm(),
                                        (st.secondPosition-centre).norm());
        const double wave=st.time-srcExtent/c;
        double cl[6]={0,0,0,0,0,0};
        double m1c=0,m1s=0,m2c=0,m2s=0,m3c=0,m3s=0;
        // In-plane basis: e1 on the separation direction, e2 = L x e1.
        const Vec3 e1=rel*(1.0/rel.norm());
        const Vec3 e2=cross(axis,e1);
        for(const SphereQuadraturePoint& p:quad){
            const Vec3 n=p.direction;
            const Vec3 obs=centre+n*R;
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
            const double flux=dot(cross(f.electric,f.magnetic)*(1.0/mu0),n);
            const double m=dot(n,axis);
            const double P[6]={1.0,m,0.5*(3*m*m-1.0),
                0.5*(5*m*m*m-3*m),
                0.125*(35*m*m*m*m-30*m*m+3.0),
                0.125*(63*m*m*m*m*m-70*m*m*m+15*m)};
            for(int l=0;l<6;++l)
                cl[l]+=p.solidAngleWeight*flux*P[l]*(2.0*l+1.0);
            const double phi=std::atan2(dot(n,e2),dot(n,e1));
            const double w=p.solidAngleWeight*flux;
            m1c+=w*std::cos(phi);      m1s+=w*std::sin(phi);
            m2c+=w*std::cos(2.0*phi);  m2s+=w*std::sin(2.0*phi);
            m3c+=w*std::cos(3.0*phi);  m3s+=w*std::sin(3.0*phi);
        }
        std::printf("%22s %9.5f %9.5f %9.5f %9.4f %9.5f\n",cs.name,
            cl[2]/cl[0],
            std::sqrt(m1c*m1c+m1s*m1s)/cl[0],
            std::sqrt(m2c*m2c+m2s*m2s)/cl[0],
            0.5*std::atan2(m2s,m2c),
            std::sqrt(m3c*m3c+m3s*m3s)/cl[0]);
    }
}
