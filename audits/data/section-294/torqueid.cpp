// Audit 294: dlaczego ta sama struktura dipolowa daje efekt 2.000 w
// PEDZIE i trzy rzedy mniejszy w MOMENCIE PEDU?  290d zapisalo to
// jako niewyjasnione.  Wyprowadzenie jest jednak krotkie.
//
// 289a zmierzylo F2 = +F1 dla czlonu dipolowego (reszta 2.0002).
// Stad dwie rozne konsekwencje dla dwoch wielkosci zachowywanych:
//   PED:          F1 + F2 = 2 F1,  czyli pelny efekt -- srodek masy
//                 przyspiesza, i to jest iloraz 2.000 z 286c/289a.
//   MOMENT PEDU:  tau = r1 x F1 + r2 x F2
//                     = (r1 + r2) x F1 + r2 x (F2 - F1).
//                 Dla ROWNYCH mas (e-/e+) zachodzi r1 + r2 = 2 r_cm,
//                 wiec
//                     tau = 2 r_cm x F1 + r2 x (F2 - F1).
//                 Para jest przygotowana z r_cm = 0 i zerowym pedem
//                 calkowitym, wiec r_cm rosnie DOPIERO w miare jak
//                 sila dipolowa ciagnie srodek masy -- stad tlumienie
//                 czynnikiem |r_cm|/|r_rel|, a nie brak efektu.
//
// ZAREJESTROWANE PRZED URUCHOMIENIEM:
//   R294a: TOZSAMOSC.  Zmierzony moment siły czlonu dipolowego ma sie
//          rownac sumie 2 r_cm x F1 + r2 x (F2-F1) do precyzji
//          maszynowej.  To jest algebra, nie fizyka -- jesli nie
//          zamknie sie, sonda jest zepsuta.
//   R294b: SKALA.  Iloraz |tau| / (|r_rel|/2 * |F1|), czyli moment
//          siły odniesiony do naiwnej skali, ma byc rzedu
//          4|r_cm|/|r_rel|, czyli maly -- to jest wyjasnienie trzech
//          rzedow z 290d.
//   R294d DODANE PO PIERWSZYM PRZEBIEGU, bo R294c padlo i powod jest
//          moim bledem normalizacji: sumowalem tau WEKTOROWO, a skale
//          naiwna SKALARNIE, wiec iloraz mieszal dwa rozne tlumienia.
//          Trzeba je rozdzielic:
//            <|tau|>  -- srednia MODULOW, czysta skala chwilowa,
//            |<tau>|  -- modul SREDNIEJ, czyli to, co zostaje po
//                        kasacji po orbicie.
//          Przewidywanie: <|tau|>/naiwne bedzie O(1) dla ORTO
//          (bo tam tau = -2 r2 x F1 i |r2| = |r_rel|/2) i rzedu
//          4|rcm|/|rrel| ~ 7e-04 dla PARA; natomiast |<tau>|/<|tau|>
//          bedzie MALE w obu kanalach i to ONO jest glownym zrodlem
//          trzech rzedow z 290d.
//   R294c: ORTO.  Tam F2 = -F1 (289a: reszta 0.0000), wiec czlon
//          2 r_cm x F1 ZNIKA z tozsamosci i caly moment siły siedzi w
//          r2 x (F2-F1) = -2 r2 x F1, ktory NIE jest tlumiony.
//          Przewidywanie: w orto iloraz z R294b bedzie O(1), czyli
//          RZEDY WIEKSZY niz w para.  Jesli tak, asymetria 290d jest
//          wyjasniona w calosci; jesli nie, wyprowadzenie jest zle.
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
    std::printf("# masy rowne: m1=%.6e m2=%.6e, iloraz=%.15f\n",
        firstMass,secondMass,firstMass/secondMass);
    std::printf("%3s %6s %11s %11s %11s %10s %11s %11s %11s\n","N",
        "kanal",
        "|<tau>|","<|tau|>","|r2 x dF|","zamkn.",
        "4|rcm|/|rrel|","<|tau|>/naiw","|<t>|/<|t|>");
    // 294e: kasacja przeciw dlugosci okna.  PREDYKCJA: czysta kasacja
    // oscylacji o zerowej sredniej da |<tau>|/<|tau|> ~ 1/N; czesc
    // SEKULARNA da podloge, ktora nie spada.
    for(int orbits : {1,2,4,8})
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
        Vec3 tauDirect{},tauCm{},tauMutual{};
        double cmOverRel=0.0,naive=0.0,tauAbs=0.0; int n=0;
        const int steps=200*orbits;
        for(int i=0;i<steps&&ok;++i){
            ok=engine.advance(st,period/200);
            if(!ok) break;
            State localState{st};
            localState.firstAcceleration={};
            localState.secondAcceleration={};
            const StateHistory localHistory{localState};
            const MutualForces f=chargeDipoleForces(st,localHistory);
            const Vec3& r1=st.firstPosition;
            const Vec3& r2=st.secondPosition;
            const Vec3 rcm=(r1*firstMass+r2*secondMass)
                /(firstMass+secondMass);
            const Vec3 rrel=r1-r2;
            const Vec3 tauNow=cross(r1,f.first)+cross(r2,f.second);
            tauDirect=tauDirect+tauNow;
            tauAbs+=tauNow.norm();
            tauCm=tauCm+cross(rcm*2.0,f.first);
            tauMutual=tauMutual+cross(r2,f.second-f.first);
            cmOverRel+=4.0*rcm.norm()/rrel.norm();
            naive+=0.5*rrel.norm()*f.first.norm();
            ++n;
        }
        if(n<1) continue;
        const double inv=1.0/n;
        const Vec3 sum=tauCm+tauMutual;
        const double closure=(tauDirect.norm()>0.0)
            ?(sum-tauDirect).norm()/tauDirect.norm():0.0;
        std::printf("%3d %6s %11.4e %11.4e %11.4e %10.2e %11.4e %11.4e"
            " %11.4e\n",orbits,
            para?"para":"orto",tauDirect.norm()*inv,tauAbs*inv,
            tauMutual.norm()*inv,closure,cmOverRel*inv,
            (naive>0.0)?tauAbs/naive:0.0,
            (tauAbs>0.0)?tauDirect.norm()/tauAbs:0.0);
    }
    std::printf("# R294a: zamkniecie ma byc ~1e-15 (algebra)\n");
    std::printf("# R294e: |<t>|/<|t|> ma spadac jak 1/N, jesli to"
        " czysta kasacja; podloga = czlon sekularny\n");
}
