// Audit 313b: czy wektory precesji obu momentow sa rowne (synchronizacja)?
// Stan startowy ziarna 40 z pierwszej linii ACTION (mu1, mu2, nhat, L).
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
    // argv: mu1x mu1y mu1z mu2x mu2y mu2z nx ny nz Lspecific
    double v[10]; for(int i=0;i<10;++i) v[i]=std::atof(argv[i+1]);
    const Vec3 m1{v[0],v[1],v[2]}, m2{v[3],v[4],v[5]}, n{v[6],v[7],v[8]};
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const Vec3 L=n*(v[9]*mu_);
    const double a=pairBohrRadius(activePair);
    auto show=[&](const char* tag,const Vec3& d1,const Vec3& d2){
        const auto r=orbitAveragedBmtAngularVelocities(a,L,d1,d2,mu_,0.0,orbitPlaneDirection(L,Vec3{1,0,0}));
        const Vec3 nh=L/L.norm();
        auto par=[&](const Vec3& w){return dot(w,nh);};
        auto perp=[&](const Vec3& w){return (w-nh*dot(w,nh)).norm();};
        std::printf("%-22s |W1|=%.4e |W2|=%.4e  W1.L=%+.4e W2.L=%+.4e  perp1=%.3e perp2=%.3e  cos(W1,W2)=%+.6f  |W1-W2|/|W1|=%.3e\n",
            tag,r.first.norm(),r.second.norm(),par(r.first),par(r.second),perp(r.first),perp(r.second),
            dot(r.first,r.second)/(r.first.norm()*r.second.norm()),(r.first-r.second).norm()/r.first.norm());
    };
    show("stan ziarna 40",m1,m2);
    show("oba momenty wzdluz L",n*m1.norm(),n*m2.norm());
    show("momenty bez wzajemnego",m1,m2*1e-12);  // usuwa pole dipola partnera
    std::printf("(rad/s; precesja spin-orbita oczekiwana ~alpha^2 w_orb = %.3e)\n",
        7.297e-3*7.297e-3*std::sqrt(pairCoulombStrength/(mu_*a*a*a)));
}
