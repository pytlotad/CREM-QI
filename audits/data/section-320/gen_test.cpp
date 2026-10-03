#include "modules/crem_collapse.hpp"
#include "modules/contact_annihilation.hpp"
#include <cstdio>
int main() {
    // (1) experiment-6 generator, fixed axis, both ways of reading the spin
    std::uint64_t st=splitMix64(42^0x5e1f7e57u);
    const int N=2000000; double p2=0, cpt=0, x=0, worst=0, worstE=0; int two=0;
    double hist[5]={0};
    for(int i=0;i<N;++i){
        auto p=drawContactAnnihilationPhotons(0.0,Vec3{0,0,1},st);
        if(p.count==2){++two;continue;}
        double c=p.normalSpinCosine; p2+=0.5*(3*c*c-1); cpt+=p.orderedSpinCorrelation; x+=p.fraction[0];
        hist[std::min(4,int(std::abs(c)*5))]+=1;
        auto cl=contactPhotonClosure(p); worst=std::max(worst,cl.momentum); worstE=std::max(worstE,cl.energy);
        // the normal recomputed from the photons themselves
    }
    std::printf("3g gen: 2g count %d, <P2> %.5f (-0.05), <S.(k1xk2)> %.5f (0), <x1> %.5f (0.66667), max|sum k|c/W %.2e, max|dE|/W %.2e\n",two,p2/N,cpt/N,x/N,worst,worstE);
    std::printf("|cos| bins (0.2 wide), measured/expected:");
    for(int b=0;b<5;++b){double a=b*0.2,c=a+0.2; double ex=((c-c*c*c/9)-(a-a*a*a/9))/(1-1.0/9); std::printf(" %.4f",hist[b]/N/ex);} std::printf("\n");
    // (2) 2 gamma branch and branching at w=1 / w=0.5
    int t1=0,t5=0; for(int i=0;i<N;++i){ if(drawContactAnnihilationPhotons(1.0,Vec3{},st).count==2)++t1; if(drawContactAnnihilationPhotons(0.5,Vec3{0,0,0.5},st).count==2)++t5;}
    double eps=orePowellSuppression();
    std::printf("P(2g): w=1 %.6f (1), w=0.5 %.6f (%.6f)\n",double(t1)/N,double(t5)/N,0.5/(0.5+0.5*eps));
    // (3) cascade sampler: per-photon marginal
    double m1=0,m2=0,mx=0; int hi=0;
    for(int i=0;i<N;++i){auto e=annihilationPhotonEnergiesFor(2.0,false,st); m1+=e[0]; m2+=e[0]*e[0]; mx+=e[1]; if(e[2]>0.9)++hi;}
    std::printf("cascade sampler (W=2 -> x): <x1> %.5f <x1^2> %.5f <x2> %.5f P(x3>0.9) %.5f  (exact 0.6667 0.5001 0.6667 0.1938)\n",m1/N,m2/N,mx/N,double(hi)/N);
}
