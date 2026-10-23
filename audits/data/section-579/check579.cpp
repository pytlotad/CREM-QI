#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <vector>
#include <cmath>
struct CaptureChannel { int n=0, l=0; double probability=0.0; };
inline bool captureModeRule() {
    static const bool on=std::getenv("CREM_CAPTURE_TRAJECTORY_HAZARD")==nullptr;
    return on;
}
inline std::vector<CaptureChannel> captureModeRuleChannels(double kineticEnergy,
        double angularMomentum,int maximumLevel=16) {
    std::vector<CaptureChannel> out;
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength;
    const double energyUnit=mu*k*k/(hbar*hbar);
    const double K=kineticEnergy/energyUnit;
    if(!(K>0.0)||!(angularMomentum>=0.0)) return out;
    const int lin=static_cast<int>(std::floor(angularMomentum/hbar));
    const double a=fineStructureConstant;
    const int nodes=200001; const double span=20.0;
    for(int lp: {lin-1,lin+1}) {
        if(lp<0||lp+1>maximumLevel) continue;
        const double lg=std::max(lin,lp);
        const double ah=1.0/(2.0*K), ecc=std::sqrt(1.0+2.0*K*lg*lg), T=std::sqrt(ah*ah*ah);
        std::vector<double> t(nodes), ax(nodes), ay(nodes), w(nodes);
        const double dF=2.0*span/(nodes-1);
        for(int i=0;i<nodes;++i) {
            const double F=-span+i*dF;
            const double x=ah*(ecc-std::cosh(F)), y=ah*std::sqrt(ecc*ecc-1.0)*std::sinh(F);
            const double r=std::sqrt(x*x+y*y), r3=r*r*r;
            t[i]=T*(ecc*std::sinh(F)-F); ax[i]=-x/r3; ay[i]=-y/r3;
            w[i]=T*(ecc*std::cosh(F)-1.0)*dF*((i==0||i==nodes-1)?0.5:1.0);
        }
        const double sign=lp==lin-1?1.0:-1.0;   // + co-rotating (Delta l = -1)
        for(int n=lp+1;n<=maximumLevel;++n) {
            const double omega=K+0.5/(n*n);
            double re=0.0, im=0.0;
            for(int i=0;i<nodes;++i) {
                const double ph=-omega*t[i], cs=std::cos(ph), sn=std::sin(ph);
                // (ax + i sign ay)(cs + i sn)
                re+=w[i]*(ax[i]*cs-sign*ay[i]*sn);
                im+=w[i]*(ax[i]*sn+sign*ay[i]*cs);
            }
            const double dEdw=(2.0*a*a*a/(3.0*pi))*(re*re+im*im)/2.0;
            const double weight=2.0*lg/(2.0*lin+1.0);
            const double P=weight*dEdw/omega/(static_cast<double>(n)*n*n);
            if(P>0.0) out.push_back({n,lp,P});
        }
    }
    return out;
}


// Audit 579: check of captureModeRuleChannels (copied from crem_interactions.cpp) against pw577.py.
int main(){
  for(int pass=0;pass<2;++pass){
    if(pass==1) applyPairFromOption("proton,electron");
    const double KoverB=pass==0?0.18:0.05;
    const double B=pairBindingEnergy(activePair), K=KoverB*B;
    std::printf("%s K/B = %.2f\n",pass==0?"Ps":"H",KoverB);
    for(int lin=0;lin<=4;++lin){
      const auto ch=captureModeRuleChannels(K,(lin+0.5)*hbar,4);
      for(const auto& c: ch) std::printf("  %d -> (%d,%d): %.4e\n",lin,c.n,c.l,c.probability);
    }
  }
}
