// Audit 238: the refusal probability as kinematics, with no
// trajectory at all.  237f item 1 asked for the e^2 >= 0 ceiling to
// be exposed so p(a, e) could be evaluated over the emission
// pattern.  Reading the inline computation (crem_collapse.hpp
// ~4850) it turns out to have a closed form, so it is reproduced
// here exactly and validated against the 0.4012 +- 0.0379 that 233a
// measured on a hundred trajectories.
//
//   angularAfter   = |L_before - hbar * helicity * nhat| / mu
//   invariantBefore= M c^2 + mu * eps
//   eps_min        = -A^2/(2 angularAfter^2)
//   invariantFloor = M c^2 + mu * eps_min
//   ceiling        = (B^2 - F^2)/(2B),  refuse when ceiling <= 0
//
// B > 0 always, so the refusal is B^2 <= F^2.  With eps_min the
// energy of the CIRCULAR orbit carrying angularAfter, that is the
// statement that the photon's vector subtraction has left more
// angular momentum than the remaining energy can hold -- which
// happens because subtracting hbar n from hbar z can REACH 2 hbar.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double M=firstMass+secondMass;
    const double A=pairCoulombStrength/mu;
    const double aPair=pairBohrRadius({electron,positron});
    // Exact refusal test, straight from the inline code.
    const auto refused=[&](double Lspec,double eps,double cosTheta,
                           double helicity){
        // L vector along +z, photon direction at polar cosTheta;
        // the azimuth drops out of the magnitude.
        const double sinTheta=std::sqrt(std::max(0.0,1.0-cosTheta*cosTheta));
        const double Lz=Lspec*mu-helicity*hbar*cosTheta;
        const double Lx=-helicity*hbar*sinTheta;
        const double angularAfter=std::sqrt(Lz*Lz+Lx*Lx)/mu;
        const double B=M*c*c+mu*eps;
        if(!(angularAfter>0.0&&B>0.0)) return true;
        const double epsMin=-A*A/(2.0*angularAfter*angularAfter);
        const double F=M*c*c+mu*epsMin;
        return !((B*B-F*F)/(2.0*B)>0.0);
    };
    // Emission pattern (1+cos^2) and the code's own helicity split.
    const auto refusalProbability=[&](double Lspec,double eps){
        const int N=200001; double num=0.0,den=0.0;
        for(int i=0;i<N;++i){
            const double ct=-1.0+2.0*i/(N-1.0);
            const double w=(1.0+ct*ct)*((i==0||i==N-1)?0.5:1.0);
            const double pPlus=std::clamp((1.0+ct)*(1.0+ct)
                /(2.0*(1.0+ct*ct)),0.0,1.0);
            const double r=pPlus*(refused(Lspec,eps,ct,+1.0)?1.0:0.0)
                          +(1.0-pPlus)*(refused(Lspec,eps,ct,-1.0)?1.0:0.0);
            num+=w*r; den+=w;
        }
        return num/den;
    };
    // The measured preparation: circular at a_pair, L = 1 hbar.
    const double epsPair=-A/(2.0*aPair);
    const double Lpair=std::sqrt(A*aPair);
    std::printf("preparacja 213b: L = %.6f hbar, eps = %.6e J/kg\n",
                Lpair*mu/hbar,epsPair);
    std::printf("p(kwadratura)       = %.6f\n",
                refusalProbability(Lpair,epsPair));
    std::printf("p(zamknieta, 27/64) = %.6f\n",27.0/64.0);
    std::printf("p(zmierzone, 233a)  = 0.4012 +- 0.0379\n");
    std::printf("\n# Granica: dla e=0 najmniejsze mozliwe L_po jest\n");
    std::printf("# |L - hbar| (foton wzdluz +z, skretnosc +), wiec kazdy\n");
    std::printf("# kierunek odmawia gdy hbar - L > L, czyli L < hbar/2,\n");
    std::printf("# czyli a < a_pair/4.  L(a) w hbar podane obok.\n");
    std::printf("\n# p(a, e): odmowa gdy L_po > L_kolowe(eps)\n");
    std::printf("%10s","a/a_pair");
    for(double e:{0.0,0.1,0.3,0.5,0.7,0.9}) std::printf(" %9.1f",e);
    std::printf("   <- mimosrod\n");
    for(double f:{1.0,0.5,0.30,0.26,0.2501,0.2499,0.24,0.1}){
        const double a=f*aPair;
        std::printf("%10.4f",f);
        for(double e:{0.0,0.1,0.3,0.5,0.7,0.9}){
            const double eps=-A/(2.0*a);
            const double L=std::sqrt(A*a*(1.0-e*e));
            std::printf(" %9.4f",refusalProbability(L,eps));
        }
        std::printf("   L=%.4f hbar\n",
                    std::sqrt(A*a)*mu/hbar);
    }
}
