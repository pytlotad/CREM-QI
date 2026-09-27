// Audit 243: where does BMT's uniform-field assumption fail?
//
// BMT is derived for a field uniform over the particle.  It carries
// no gradient term in the spin equation.  The natural dimensionless
// measure of the neglected term is the field's relative variation
// over the spin's own length scale, the reduced Compton wavelength:
//     eta = lambda_C |grad B| / |B|.
// For a pure dipole field B ~ 1/r^3 that is 3 lambda_C / r, which at
// r* = 193.3035 fm is 5.99 -- order unity and then some.  Measured
// here on the model's OWN field, which carries Plummer softening and
// retardation that the 1/r^3 estimate does not.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    const double lambdaC=hbar/(electron.mass*c);
    const double rstar=comptonBarrierRadius;
    std::printf("# lambda_C = %.6e m, r* = %.6e m, lambda_C/r* = %.4f\n",
                lambdaC,rstar,lambdaC/rstar);
    std::printf("%10s %12s %12s %12s %12s\n",
                "r/a_pair","r [fm]","|B| [T]","eta","3 lC/r");
    for(double f:{1.0,0.25,0.05,0.01,4.0e-3,2.0e-3,1.0e-3,
                  rstar/A,5.0e-4}){
        const double a=f*A;
        if(!(a>2.0*nuclearCutoff)) continue;
        const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}*firstMagneticMoment;
        const State start=osculatingPeriapsisState(
            el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=26;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start; bool ok=true;
        const int n=200;
        for(int i=0;i<n&&ok;++i) ok=engine.advance(s,0.25*period/n);
        if(!ok){ std::printf("%10.2e failed\n",f); continue; }
        const StateHistory& h=engine.history();
        const Vec3 p=s.firstPosition;
        const double step=1.0e-4*(s.firstPosition-s.secondPosition).norm();
        const auto Bat=[&](const Vec3& x){
            return fieldFromOtherParticleAt(x,s.time,s,h,true).magnetic; };
        const Vec3 B0=Bat(p);
        double grad2=0.0;
        for(int ax=0;ax<3;++ax){
            Vec3 d{}; (ax==0?d.x:ax==1?d.y:d.z)=step;
            const Vec3 dB=(Bat(p+d)-Bat(p-d))*(1.0/(2.0*step));
            grad2+=dB.squaredNorm();
        }
        const double gradB=std::sqrt(grad2);
        const double sep=(s.firstPosition-s.secondPosition).norm();
        const double eta=B0.norm()>0.0?lambdaC*gradB/B0.norm():0.0;
        std::printf("%10.2e %12.3f %12.4e %12.4e %12.4e\n",
                    f,sep*1e15,B0.norm(),eta,3.0*lambdaC/sep);
    }
    std::printf("# eta >= 1 oznacza, ze pominiety czlon gradientowy\n");
    std::printf("# jest tego samego rzedu co zachowany.\n");
}
