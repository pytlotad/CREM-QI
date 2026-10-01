// Audit 290: przypisanie per czlon dla sciezki ZACHOWAWCZEJ.
//
// Tabela z 280-282 byla mierzona na sile retardowanej, ktora po 288
// jest domyslna tylko dla trybow CIAGLYCH; tryb stochastyczny ma
// teraz sciezke zachowawcza.  Tamta tabela opisuje wiec jedna z dwoch
// konfiguracji produkcyjnych, a drugiej brakuje -- 289 dalo dla niej
// tylko kolumne pedu.
//
// Przewidywanie zapisane przed pomiarem: bez retardowanego czlonu
// przyspieszeniowego JEDYNYM momentem radiacyjnym jest reakcja
// wlasna, a przy bramce z 278h nie zapadajacej sie do ll powinna byc
// pelna, czyli okolo -1.00 strumienia, i domkniecie wyjdzie okolo
// 1.00 zamiast 1.037.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu_=firstMass*secondMass/(firstMass+secondMass);
    const double kk=pairCoulombStrength/mu_;
    const double A=pairBohrRadius({electron,positron});
    const OsculatingElements el{-kk/(2.0*A),std::sqrt(kk*A)};
    const double period=osculatingPeriod(el.specificEnergy,kk);
    FarFieldSampling sampling;
    sampling.directionCount=194;
    sampling.controlRadius=1.0e6*bohrRadius;
    for(int channel=0;channel<2;++channel){
        const bool para=(channel==0);
        const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
        const Vec3 m2=Vec3{0,0,para?1.0:-1.0}*secondMagneticMoment;
        const State start=osculatingPeriapsisState(el,kk,m1,m2,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
        acc.reactionModel=
            ChargeRadiationReactionModel::individualLandauLifshitz;
        acc.computeOutwardFlux=true;
        acc.farFieldSampling=sampling;
        acc.useRetardedExternalForces=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start; bool ok=true;
        for(int i=0;i<200&&ok;++i) ok=engine.advance(st,period/400);
        if(!ok){ std::printf("%s przerwane\n",para?"para":"orto");
                 continue; }
        enum { COUL,DARW,CDIP,SELF,COUNT };
        const char* names[COUNT]={"mutualForces (Coulomb)","darwinForces",
            "chargeDipoleForces","reakcja wlasna (LL)"};
        Vec3 torque[COUNT]{}; double power[COUNT]{};
        Vec3 fluxAngular{}; double fluxEnergy=0.0; int n=0;
        for(int i=0;i<200&&ok;++i){
            ok=engine.advance(st,period/400);
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
            const ParticleMultipoleRadiation radiation=
                particleMultipoleRadiation(st,total,h,true,
                    ChargeRadiationReactionModel::individualLandauLifshitz,
                    false,sampling);
            const auto add=[&](int slot,const Vec3& f1,const Vec3& f2){
                torque[slot]=torque[slot]+cross(st.firstPosition,f1)
                                         +cross(st.secondPosition,f2);
                power[slot]+=dot(f1,st.firstVelocity)
                            +dot(f2,st.secondVelocity);
            };
            add(COUL,coulomb.first,coulomb.second);
            add(DARW,darwin.first,darwin.second);
            add(CDIP,chargeDipole.first,chargeDipole.second);
            add(SELF,radiation.chargeReaction.first,
                     radiation.chargeReaction.second);
            fluxAngular=fluxAngular+radiation.outwardFlux.angularMomentum;
            fluxEnergy+=radiation.outwardFlux.energy;
            ++n;
        }
        if(n<1) continue;
        const double inv=1.0/n;
        std::printf("\n=== %s, sciezka ZACHOWAWCZA, %d probek ===\n",
            para?"PARA":"ORTO",n);
        std::printf("%26s %13s %13s %10s\n","czlon","moment z","moc [W]",
            "/strumien");
        Vec3 torqueSum{}; double powerSum=0.0;
        const double scale=std::abs(fluxAngular.z)>0.0
            ?std::abs(fluxAngular.z*inv):1.0;
        for(int k=0;k<COUNT;++k){
            std::printf("%26s %13.4e %13.4e %10.4f\n",names[k],
                torque[k].z*inv,power[k]*inv,torque[k].z*inv/scale);
            torqueSum=torqueSum+torque[k]; powerSum+=power[k];
        }
        std::printf("%26s %13.4e %13.4e %10.4f\n","SUMA",
            torqueSum.z*inv,powerSum*inv,torqueSum.z*inv/scale);
        std::printf("%26s %13.4e %13.4e %10.4f\n","strumien pola",
            fluxAngular.z*inv,fluxEnergy*inv,1.0);
        std::printf("%26s %13.4e %10s %10.4f\n","RESIDUUM",
            (torqueSum.z+fluxAngular.z)*inv,"",
            (torqueSum.z+fluxAngular.z)*inv/scale);
    }
}
