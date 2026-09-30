// Audit 255: is 221f's k=2 floor resolved, and does the interpolant
// change it?
//
// Two questions in one scan, following 254e.
//
// (1) TOLERANCE.  221f ran at relativeTolerance 1e-10 only and 221h
// recorded its minimum, 6.141e-06, as "where the scan landed, not a
// measured floor".  If the k=2 amplitudes move when the tolerance
// moves two decades, they are the integrator, not the orbit.
//
// (2) INTERPOLANT.  A fact found while writing this probe: 221f
// built its spectrum from st.firstAcceleration, the integrator's
// STORED acceleration, which never passes through historicalCharge.
// So 254's first-order interpolant error does not touch 221f's own
// numbers -- but it does touch what the model RADIATES, because
// farZoneChargeField takes its acceleration from historicalCharge
// (253b, 254a).  The second column repeats the same harmonic
// decomposition on the interpolated acceleration at the same
// instants, so the two columns differ by the interpolant alone.
//
// The interpolant is sampled at the present time, which is the END
// of its last segment -- the worst position, per 254b.  It is the
// same place 254 measured, and it is an upper bound on what the
// retarded evaluations mid-segment actually see.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double a=pairBohrRadius({electron,positron});
    const double b2=pairCoulombStrength/(mu*c*c*a);
    const int samples=256;
    std::printf("# a = a_pair, beta_rel^2 = %.4e, momenty wylaczone, "
        "%d probek na obieg\n",b2,samples);
    std::printf("%12s %10s %13s %13s %13s %9s\n",
        "v/v_Coulomb","tolerancja","e realized","k=2 zapisane",
        "k=2 interp","iloraz");
    for(double f:{-0.5,-0.25,-0.125,0.0,0.25,1.0}){
        const double s0=1.0+f*b2;
        for(double tol:{1.0e-8,1.0e-10,1.0e-12}){
            const double period=osculatingPeriod(-k/(2.0*a),k);
            const double v=std::sqrt(k/a)*s0;
            const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*1.0e-6);
            State start;
            const double w1=secondMass/(firstMass+secondMass);
            const double w2=-firstMass/(firstMass+secondMass);
            start.firstPosition=Vec3{a,0,0}*w1;
            start.secondPosition=Vec3{a,0,0}*w2;
            start.firstVelocity=Vec3{0,v,0}*w1;
            start.secondVelocity=Vec3{0,v,0}*w2;
            start.firstDipole=m1;
            start.secondDipole=m1*(secondMagneticMoment/firstMagneticMoment);
            ClassicalTrajectoryEngine::Accuracy acc;
            acc.relativeTolerance=tol; acc.maximumDepth=26;
            acc.reactionModel=ChargeRadiationReactionModel::disabled;
            acc.computeOutwardFlux=false;
            ClassicalTrajectoryEngine engine(start,acc);
            State st=start;
            double lo=1e300,hi=0.0; bool ok=true;
            std::vector<Vec3> stored,interp;
            stored.reserve(samples); interp.reserve(samples);
            for(int i=0;i<samples&&ok;++i){
                ok=engine.advance(st,period/samples);
                if(!ok) break;
                stored.push_back(st.firstAcceleration*firstCharge
                                +st.secondAcceleration*secondCharge);
                const StateHistory& h=engine.history();
                interp.push_back(
                    historicalCharge(h,st,true,st.time).acceleration*firstCharge
                   +historicalCharge(h,st,false,st.time).acceleration*secondCharge);
                const double r=(st.firstPosition-st.secondPosition).norm();
                lo=std::min(lo,r); hi=std::max(hi,r);
            }
            if(!ok){ std::printf("%12.9f %10.0e  nieudane\n",s0,tol); continue; }
            const auto ratio=[&](const std::vector<Vec3>& s){
                double amp[2]={0,0};
                for(int kk=1;kk<=2;++kk){
                    std::complex<double> cx{},cy{},cz{};
                    for(int i=0;i<samples;++i){
                        const double ph=-2.0*pi*kk*i/samples;
                        const std::complex<double> w{std::cos(ph),std::sin(ph)};
                        cx+=w*s[i].x; cy+=w*s[i].y; cz+=w*s[i].z;
                    }
                    amp[kk-1]=std::sqrt(std::norm(cx)+std::norm(cy)
                                       +std::norm(cz));
                }
                return amp[0]>0.0?amp[1]/amp[0]:0.0;
            };
            const double rs=ratio(stored),ri=ratio(interp);
            std::printf("%12.9f %10.0e %13.6e %13.6e %13.6e %9.4f\n",
                s0,tol,(hi-lo)/(hi+lo),rs,ri,rs>0.0?ri/rs:0.0);
        }
    }
}
