#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
static Vec3 v3(const char* s){ Vec3 v; std::sscanf(s,"%lf,%lf,%lf",&v.x,&v.y,&v.z); return v; }
int main(int argc,char**argv){
    // args: m1 m2 Lhat Lafter_hbar Egamma_J a_before
    const Vec3 m1=v3(argv[1]), m2=v3(argv[2]); Vec3 Lh=v3(argv[3]); Lh=Lh*(1/Lh.norm());
    const double Lafter=std::atof(argv[4])*hbar, Eg=std::atof(argv[5]), abefore=std::atof(argv[6]);
    const double mu=firstMass*secondMass/(firstMass+secondMass), K=pairCoulombStrength/mu;
    const double epsBefore=-K/(2*abefore), eps=epsBefore-Eg/mu;
    OsculatingElements el{eps,Lafter/mu};
    const double rp=osculatingPeriapsis(el,K), ra=osculatingApoapsis(el,K);
    const double dap=dipoleAwarePeriapsis(el,K,m1,m2,Lh*Lafter,mu);
    std::printf("after photon: a=%.4e m, L=%.3e hbar, Kepler peri=%.4e r*, apo=%.4e r*, dipoleAwarePeriapsis=%.4e r*\n",
        -K/(2*eps),Lafter/hbar,rp/comptonBarrierRadius,ra/comptonBarrierRadius,dap/comptonBarrierRadius);
    const double totalM=firstMass+secondMass;
    for(double x: {1e-4,1e-3,1e-2,0.1,0.3,0.5,0.8,1.0,1.2,1.5,2.0,3.0,10.0,100.0,300.0,360.0,364.0,365.0}){
        const double r=x*comptonBarrierRadius; if(r>ra) continue;
        const double dip=azimuthAveragedDipoleEnergy(r,m1,m2,Lh*Lafter)/mu;
        const Vec3 seed=std::abs(Lh.x)<0.9?Vec3{1,0,0}:Vec3{0,1,0}; Vec3 rh=cross(seed,Lh); rh=rh*(1/rh.norm()); const Vec3 th=cross(Lh,rh);
        State tp{}; const double vt=(Lafter/mu)/r;
        tp.firstPosition=rh*(r*secondMass/totalM); tp.secondPosition=rh*(-r*firstMass/totalM);
        tp.firstVelocity=th*(vt*secondMass/totalM); tp.secondVelocity=th*(-vt*firstMass/totalM);
        tp.firstProperDipole=m1; tp.secondProperDipole=m2; synchronizeCovariantDipoles(tp);
        const double so=chargeDipoleInteractionEnergy(tp)/mu;
        const double h=eps+K/r-dip-so-(Lafter/mu)*(Lafter/mu)/(2*r*r);
        std::printf("  r=%9.4g r*  h=%+.4e  (Coulomb+E %.3e, dipole %.3e, spin-orbit %.3e) [J/kg]\n",x,h,eps+K/r,dip,so);
    }
}
