// Audit 289: brama (iii) -- czy sciezka ZACHOWAWCZA ma te sama
// asymetrie para/orto w pedzie, co retardowana?
//
// 286c zmierzylo, ze niezbilansowanie pedu w sektorze dipolowym jest
// SPECYFICZNE DLA PARA: oba czlony retardowane daja reszte trzeciej
// zasady 2.000 i 8.027 dla para, a DOKLADNIE 0 dla orto.  288d
// zapisalo to jako powod, zeby brame (iii) sprawdzic osobno.
//
// Ale przelacznik nie usuwa sektora dipolowego -- zamienia go.
// retardedExternalForces niesie sile Lorentza od retardowanego pola
// momentu plus covariantDipoleGradientForce; allExternalForces niesie
// CHWILOWE chargeDipoleForces (plus mutualForces i darwinForces).
// Pytanie jest wiec konkretne i tanie: czy sciezka zachowawcza ma te
// sama asymetrie.  Jesli NIE, przelacznik usuwa mechanizm asymetrii i
// brama (iii) jest zagrozona niezaleznie od dlugosci przebiegu.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const double a=0.25*A;
    const OsculatingElements el{-kk/(2.0*a),std::sqrt(kk*a)};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    std::printf("%6s %28s %13s %13s %11s\n","kanal","czlon",
        "|F1+F2|","|F1|","reszta/|F1|");
    for(int channel=0;channel<2;++channel){
        const bool para=(channel==0);
        const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
        const Vec3 m2=Vec3{0,0,para?1.0:-1.0}*secondMagneticMoment;
        const State start=osculatingPeriapsisState(el,kk,m1,m2,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        acc.useRetardedExternalForces=false;   // sciezka zachowawcza
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<100&&ok;++i) ok=engine.advance(st,period/200);
        if(!ok){ std::printf("%6s przerwane\n",para?"para":"orto");
                 continue; }
        Vec3 coulombFirst{},coulombSecond{};
        Vec3 darwinFirst{},darwinSecond{};
        Vec3 chargeDipoleFirst{},chargeDipoleSecond{};
        Vec3 totalFirst{},totalSecond{};
        int n=0;
        for(int i=0;i<100&&ok;++i){
            ok=engine.advance(st,period/200);
            if(!ok) break;
            const StateHistory& h=engine.history();
            const MutualForces coulomb=mutualForces(st);
            const MutualForces darwin=darwinForces(st);
            State localState{st};
            localState.firstAcceleration={};
            localState.secondAcceleration={};
            const StateHistory localHistory{localState};
            const MutualForces chargeDipole=
                chargeDipoleForces(st,localHistory);
            const MutualForces total=allExternalForces(st);
            coulombFirst=coulombFirst+coulomb.first;
            coulombSecond=coulombSecond+coulomb.second;
            darwinFirst=darwinFirst+darwin.first;
            darwinSecond=darwinSecond+darwin.second;
            chargeDipoleFirst=chargeDipoleFirst+chargeDipole.first;
            chargeDipoleSecond=chargeDipoleSecond+chargeDipole.second;
            totalFirst=totalFirst+total.first;
            totalSecond=totalSecond+total.second;
            (void)h;
            ++n;
        }
        if(n<1) continue;
        const double inv=1.0/n;
        const auto report=[&](const char* name,const Vec3& f1,
                              const Vec3& f2){
            const Vec3 residual=(f1+f2)*inv;
            const double scale=f1.norm()*inv;
            std::printf("%6s %28s %13.4e %13.4e %11.4f\n",
                para?"para":"orto",name,residual.norm(),scale,
                scale>0.0?residual.norm()/scale:0.0);
        };
        report("mutualForces (Coulomb)",coulombFirst,coulombSecond);
        report("darwinForces",darwinFirst,darwinSecond);
        report("chargeDipoleForces (chwil.)",chargeDipoleFirst,
               chargeDipoleSecond);
        report("allExternalForces (SUMA)",totalFirst,totalSecond);
    }
    std::printf("# 286c dla sciezki retardowanej: para 2.000 i 8.027,"
        " orto 0.000 i 0.000\n");
}
