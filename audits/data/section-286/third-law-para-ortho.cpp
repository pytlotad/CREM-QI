// Audit 286: dwa punkty zostawione w 282f.
//
// (1) ORTO.  Wyprowadzenie F_1 = F_2 z 282b dotyczy momentow
//     ROWNOLEGLYCH.  Dla orto m_2 = -m_1, a pole dipolowe jest
//     liniowe w m i parzyste w r, wiec B(r_21,m_1) = -B(r_12,m_2) i
//     F_2 = (-q_1)(-v_1) x (-B) = -F_1: trzecia zasada ZACHODZI i
//     iloraz |F1+F2|/|F1| ma wyjsc 0, nie 2.  Przewidywanie
//     jakosciowo odwrotne.
//
// (2) GRADIENT.  Dla F = grad(m.B): m.B(r) jest parzyste w r, wiec
//     gradient jest nieparzysty, wiec dla para F_2 = -F_1 i iloraz
//     ma byc 0.  Zmierzono 3.202.  Sprzecznosc jest informacja: albo
//     covariantDipoleGradientForce nie jest prostym grad(m.B), albo
//     retardacja lamie parzystosc -- ale wtedy nie moze byc
//     jednoczesnie zaniedbywalna dla czesci Lorentza, gdzie wyszlo
//     dokladnie 2.000.
//
// Kontrola: pairDipoleForce jest sila CHWILOWA i wzajemna z
// konstrukcji, wiec jej iloraz musi byc 0 w obu kanalach.  Jesli nie
// jest, blad jest w moim pomiarze, nie w modelu.
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
    std::printf("%6s %26s %13s %13s %11s\n","kanal","czlon",
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
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<100&&ok;++i) ok=engine.advance(st,period/200);
        if(!ok){ std::printf("%6s przerwane\n",para?"para":"orto");
                 continue; }
        Vec3 lorentzFirst{},lorentzSecond{};
        Vec3 gradientFirst{},gradientSecond{};
        Vec3 instantFirst{},instantSecond{};
        int n=0;
        for(int i=0;i<100&&ok;++i){
            ok=engine.advance(st,period/200);
            if(!ok) break;
            const StateHistory& h=engine.history();
            const ElectromagneticField dipAtFirst=
                retardedMagneticDipoleField(st.firstPosition,st.time,h,
                                            st,false);
            const ElectromagneticField dipAtSecond=
                retardedMagneticDipoleField(st.secondPosition,st.time,h,
                                            st,true);
            lorentzFirst=lorentzFirst
                +lorentzForce(firstCharge,st.firstVelocity,dipAtFirst);
            lorentzSecond=lorentzSecond
                +lorentzForce(secondCharge,st.secondVelocity,dipAtSecond);
            gradientFirst=gradientFirst
                +covariantDipoleGradientForce(st,h,true);
            gradientSecond=gradientSecond
                +covariantDipoleGradientForce(st,h,false);
            const Vec3 separationVector=st.firstPosition-st.secondPosition;
            instantFirst=instantFirst
                +pairDipoleForce(separationVector,st.firstDipole,
                                 st.secondDipole);
            instantSecond=instantSecond
                +pairDipoleForce(separationVector*-1.0,st.secondDipole,
                                 st.firstDipole);
            ++n;
        }
        if(n<1) continue;
        const double inv=1.0/n;
        const auto report=[&](const char* name,const Vec3& f1,
                              const Vec3& f2){
            const Vec3 residual=(f1+f2)*inv;
            const double scale=f1.norm()*inv;
            std::printf("%6s %26s %13.4e %13.4e %11.4f\n",
                para?"para":"orto",name,residual.norm(),scale,
                scale>0.0?residual.norm()/scale:0.0);
        };
        report("dipol Lorentz (retard.)",lorentzFirst,lorentzSecond);
        report("gradient kowariantny",gradientFirst,gradientSecond);
        report("pairDipoleForce (chwil.)",instantFirst,instantSecond);
    }
    std::printf("# przewidywania: Lorentz para 2.000, orto 0;"
        " chwilowa 0 w obu\n");
}
