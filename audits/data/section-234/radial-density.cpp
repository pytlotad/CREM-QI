// Audit 234: <n(0)> read off the radial orbit the full-hbar photon
// leaves behind.  13a factorizes the singlet rate as
// Gamma = 4 (pi r_e^2 c) <n(0)> and 13b names <n(0)> as the one
// factor the model cannot supply.  228-230 showed the model DOES
// reach a radial orbit, for about three periods at the end of the
// waiting time, so the density can now be read from it directly:
// the fraction of a radial period spent inside R, divided by the
// volume of that ball.
//
// The classical expectation is the point.  For a radial Kepler orbit
// v ~ sqrt(2K/(mu r)) near the origin, so the time in a shell goes as
// r^(1/2) dr against a volume 4 pi r^2 dr, and the density diverges
// as r^(-3/2).  A mean over a ball of radius R diverges the same way.
// Whether the model's floors cure that is what is measured.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    // The post-photon orbit of audit 229c, seed 1: a = 0.4127 a_pair
    // and L = 3.41e-05 hbar.  Launched at apoapsis r = 2a, where a
    // radial orbit is momentarily at rest bar the tiny tangential L/r.
    const double a=0.4127*A;
    const double Lspec=3.40959940739912388e-05*hbar/mu;
    const double rApo=2.0*a;
    const double vTan=Lspec/rApo;
    State s;
    const double w1=secondMass/(firstMass+secondMass);
    const double w2=-firstMass/(firstMass+secondMass);
    s.firstPosition=Vec3{rApo,0,0}*w1;
    s.secondPosition=Vec3{rApo,0,0}*w2;
    s.firstVelocity=Vec3{0,vTan,0}*w1;
    s.secondVelocity=Vec3{0,vTan,0}*w2;
    const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
    s.firstDipole=m1;
    s.secondDipole=m1*(secondMagneticMoment/firstMagneticMoment);
    const double period=2.0*pi*std::sqrt(a*a*a/k);
    std::printf("# a = %.4f a_pair, r_apo = %.4e m, T_rad = %.4e s\n",
                a/A,rApo,period);
    std::printf("# L = 3.4096e-05 hbar, v_tan at apoapsis = %.4e m/s\n",
                vTan);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-9; acc.maximumDepth=20;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(s,acc);
    // Thresholds in metres, from a_pair down past the Compton barrier.
    std::vector<double> R;
    for(double e=-11.0;e>=-12.51;e-=0.25) R.push_back(std::pow(10.0,e));
    std::vector<double> tCross(R.size(),-1.0);
    double t=0.0; const double dt=period/20000.0;
    double rMin=1e300;
    for(long step=0;step<40000;++step){
        if(!engine.advance(s,dt)) break;
        t+=dt;
        const double r=(s.firstPosition-s.secondPosition).norm();
        rMin=std::min(rMin,r);
        for(std::size_t i=0;i<R.size();++i)
            if(tCross[i]<0.0&&r<R[i]) tCross[i]=t;
        if(r<3.0e-13) break;
    }
    std::printf("# deepest separation reached %.4e m after %.4e s\n",
                rMin,t);
    std::printf("\n%12s %14s %14s %14s %12s\n",
                "R [m]","t inside [s]","frac of T","n(R) [1/m^3]",
                "/(1/pi a^3)");
    const double quantum=1.0/(pi*A*A*A);
    for(std::size_t i=0;i<R.size();++i){
        if(tCross[i]<0.0) continue;
        const double tin=2.0*(t-tCross[i]);   // in and out
        const double frac=tin/period;
        const double vol=4.0/3.0*pi*R[i]*R[i]*R[i];
        std::printf("%12.3e %14.6e %14.6e %14.6e %12.4e\n",
                    R[i],tin,frac,frac/vol,frac/vol/quantum);
    }
    std::printf("\n# quantum |psi(0)|^2 = 1/(pi a_Ps^3) = %.6e 1/m^3\n",
                quantum);
    std::printf("# a classical radial orbit predicts n ~ R^(-3/2).\n");
}
