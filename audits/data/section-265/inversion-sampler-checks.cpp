// Audit 265: the fixed-consumption sampler.  Three new checks.
//
//   ZUZYCIE      -- exactly two uniforms per draw, always.  That is the
//                   whole point: it restores seed-for-seed comparison.
//   ZGODNOSC     -- the inversion and rejection samplers must produce the
//                   SAME distribution.  Compared through the second-moment
//                   tensor, the same estimator 263 used, on a pattern with
//                   all five moments nonzero.
//   ZAWIERANIE   -- the claim that the new sampler contains the old
//                   prescribed draw.  For c0 = 4/3, cm0 = 2/3 the polar
//                   marginal's cubic IS the depressed cubic
//                   mu^3+3mu+(4-8u)=0 that crem_collapse solves by
//                   Cardano, so for the same u the two must agree.  If
//                   they do not, the derivation is wrong.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double c0=1.0;
    const double cm0=0.40,cm2=0.30,cxy=0.25,cyz=0.20,cxz=0.15;
    AngularPatternMoments m;
    m.a00=4.0*pi*c0;
    m.a2[0]=cxy*(4.0*pi/15.0);  m.a2[1]=cyz*(4.0*pi/15.0);
    m.a2[2]=cm0*(4.0*pi/5.0);   m.a2[3]=cxz*(4.0*pi/15.0);
    m.a2[4]=cm2*(4.0*pi/15.0);

    // --- ZUZYCIE ---
    std::mt19937_64 rng(11111);
    std::uniform_real_distribution<double> u(0.0,1.0);
    int calls=0;
    const auto counted=[&]{ ++calls; return u(rng); };
    int worst=0,best=1000;
    for(int i=0;i<200000;++i){
        calls=0;
        const Vec3 d=drawDirectionFromPattern(m,counted);
        if(d.squaredNorm()==0.0){ std::printf("nieudane losowanie!\n"); return 1; }
        worst=std::max(worst,calls); best=std::min(best,calls);
    }
    std::printf("ZUZYCIE: min=%d max=%d  (oczekiwane 2 i 2)\n",best,worst);

    // --- ZGODNOSC dwoch samplerow ---
    const auto tensor=[&](bool inversion,int n,double out[3][3]){
        std::mt19937_64 r(2024);
        std::uniform_real_distribution<double> uu(0.0,1.0);
        const auto dr=[&]{ return uu(r); };
        double acc[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        long long got=0;
        for(int i=0;i<n;++i){
            const Vec3 d=inversion
                ?drawDirectionFromPatternByInversion(m,dr)
                :drawDirectionFromPatternByRejection(m,dr);
            if(d.squaredNorm()==0.0) continue;
            ++got;
            const double v[3]={d.x,d.y,d.z};
            for(int a=0;a<3;++a) for(int b=0;b<3;++b) acc[a][b]+=v[a]*v[b];
        }
        for(int a=0;a<3;++a) for(int b=0;b<3;++b)
            out[a][b]=(15.0*c0/2.0)*(acc[a][b]/got-(a==b?1.0/3.0:0.0));
    };
    double Qinv[3][3],Qrej[3][3];
    tensor(true,4000000,Qinv);
    tensor(false,4000000,Qrej);
    const double Qtrue[3][3]={
        {-0.5*cm0+0.5*cm2, 0.5*cxy,          0.5*cxz},
        { 0.5*cxy,         -0.5*cm0-0.5*cm2, 0.5*cyz},
        { 0.5*cxz,          0.5*cyz,         cm0}};
    double wInv=0,wRej=0,wPair=0;
    for(int a=0;a<3;++a) for(int b=0;b<3;++b){
        wInv=std::max(wInv,std::abs(Qinv[a][b]-Qtrue[a][b]));
        wRej=std::max(wRej,std::abs(Qrej[a][b]-Qtrue[a][b]));
        wPair=std::max(wPair,std::abs(Qinv[a][b]-Qrej[a][b]));
    }
    std::printf("ZGODNOSC: odwrotnosciowy vs analityczny %.3e,"
        " odrzuceniowy vs analityczny %.3e, miedzy soba %.3e\n",
        wInv,wRej,wPair);

    // --- ZAWIERANIE: czy dla c0=4/3, cm0=2/3 wychodzi korzen Cardana ---
    AngularPatternMoments pres;
    pres.a00=4.0*pi*(4.0/3.0);
    pres.a2[2]=(4.0*pi/5.0)*(2.0/3.0);
    const auto signedCbrt=[](double v){
        return std::copysign(std::cbrt(std::abs(v)),v); };
    double worstCardano=0.0;
    for(int i=1;i<1000;++i){
        const double uu=double(i)/1000.0;
        const double q=4.0-8.0*uu;
        const double disc=(q*q)/4.0+1.0;
        const double root=std::sqrt(disc);
        const double cardano=signedCbrt(-q/2.0+root)+signedCbrt(-q/2.0-root);
        double seq[2]={uu,0.5}; int k=0;
        const auto fixed=[&]{ return seq[k++]; };
        const Vec3 d=drawDirectionFromPatternByInversion(pres,fixed);
        worstCardano=std::max(worstCardano,std::abs(d.z-cardano));
    }
    std::printf("ZAWIERANIE: najgorsza roznica mu od korzenia Cardana"
        " na 999 wartosciach u: %.3e\n",worstCardano);
}
