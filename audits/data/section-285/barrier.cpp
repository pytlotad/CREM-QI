// Audit 285: frakcja czasu w barierze na orbicie radialnej.
//
// Pierwszy pomiar wspolnego modelu emisji: kanal kontaktowy
// (anihilacja) ma miec tempo proporcjonalne do czasu spedzanego w
// poblizu zera, a model ma twardy rdzen na r* = (g/2) lambda_C/2,
// wiec orbita nigdy nie dochodzi do r=0 -- dochodzi do bariery.
// Czyli <n(0)> jest w istocie <n(r*)>.
//
// Jest na to sprawdzalna liczba.  Tozsamosc alpha^-5 z rejestru daje
//   tau_2gamma = 2 lambda_C/(alpha^5 c) = 4 r*/(alpha^5 c),
// czyli tau_2gamma c/r* = 4/alpha^5.  Jesli tempo kontaktowe wynika z
// rezydencji w barierze, to zmierzona frakcja f i polokres T musza
// dawac 1/(f Gamma_in) = tau_2gamma z Gamma_in rzedu jednosci w
// jednostkach c/r*.  Jesli Gamma_in wyjdzie mikroskopijne, rezydencja
// SAMA tego nie wyjasnia i brakujacy czynnik trzeba bedzie nazwac.
//
// Orbita: spadek ze spoczynku z a_pair, czyli L=0 i r_max=a_pair --
// konfiguracja, ktora 228 nazwalo jedyna z kontaktem.  Reakcja
// wylaczona, zeby mierzyc czysta rezydencje zachowawcza.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double rStar=comptonBarrierRadius;
    // Kepler: polokres radialny dla r_max=A to polowa okresu elipsy o
    // polosi A/2.
    const double halfPeriod=pi*std::sqrt(std::pow(0.5*A,3)/kk);
    std::printf("# a_pair=%.6e m, r*=%.6e m, r*/a_pair=%.4e\n",
        A,rStar,rStar/A);
    std::printf("# polokres radialny Keplera=%.6e s = %.4e r*/c\n",
        halfPeriod,halfPeriod/(rStar/c));
    // Stan: oba ladunki w spoczynku, separacja a_pair.
    const double w1=secondMass/(firstMass+secondMass);
    const double w2=-firstMass/(firstMass+secondMass);
    State start;
    start.firstPosition=Vec3{A,0,0}*w1;
    start.secondPosition=Vec3{A,0,0}*w2;
    const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*1.0e-6);
    start.firstDipole=m1;
    start.secondDipole=m1*(secondMagneticMoment/firstMagneticMoment);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-9;
    acc.maximumDepth=std::getenv("CREM_MAX_DEPTH")
        ?std::atoi(std::getenv("CREM_MAX_DEPTH")):30;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State st=start;
    const double thresholds[]={0.5,1.0,2.0,5.0,10.0,100.0};
    double inside[6]={0,0,0,0,0,0};
    double minimumSeparation=1e300,elapsed=0.0;
    const double step=halfPeriod/200000.0;
    bool ok=true; long long steps=0;
    while(ok&&elapsed<1.5*halfPeriod){
        const double before=(st.firstPosition-st.secondPosition).norm();
        ok=engine.advance(st,step);
        if(!ok) break;
        const double after=(st.firstPosition-st.secondPosition).norm();
        const double mid=0.5*(before+after);
        elapsed+=step; ++steps;
        minimumSeparation=std::min(minimumSeparation,after);
        for(int k=0;k<6;++k)
            if(mid<thresholds[k]*rStar) inside[k]+=step;
    }
    std::printf("# krokow=%lld, czas=%.6e s (%.3f polokresu), %s\n",
        steps,elapsed,elapsed/halfPeriod,ok?"ok":"PRZERWANE");
    std::printf("# minimalna separacja=%.6e m = %.4f r*\n",
        minimumSeparation,minimumSeparation/rStar);
    std::printf("\n%10s %14s %14s\n","prog [r*]","czas w [s]","frakcja");
    for(int k=0;k<6;++k)
        std::printf("%10.1f %14.6e %14.6e\n",thresholds[k],inside[k],
            elapsed>0.0?inside[k]/elapsed:0.0);
    // Czego wymagaloby tau_2gamma = 4 r*/(alpha^5 c).
    const double alpha=pairCoulombStrength/(hbar*c);
    const double tau2g=4.0*rStar/(std::pow(alpha,5)*c);
    std::printf("\n# alpha=1/%.6f, tau_2gamma=4r*/(alpha^5 c)=%.6e s\n",
        1.0/alpha,tau2g);
    std::printf("# tau_2gamma c/r* = %.4e  (cel: 4/alpha^5)\n",
        tau2g*c/rStar);
    for(int k=0;k<6;++k){
        const double f=elapsed>0.0?inside[k]/elapsed:0.0;
        if(!(f>0.0)) continue;
        const double gammaIn=1.0/(f*tau2g);
        std::printf("# prog %5.1f r*: f=%.4e -> Gamma_in=%.4e 1/s"
            " = %.4e c/r*\n",thresholds[k],f,gammaIn,
            gammaIn*(rStar/c));
    }
}
