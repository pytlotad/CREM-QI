// Audit 287: skad ujemny znak -0.0063 z 284f?
//
// Przy sile ZACHOWAWCZEJ w trybie stochastycznym iloraz
// -dE_mech/E_wypromieniowana wyszedl -0.0063, czyli orbita formalnie
// ZYSKUJE 0.63 procenta energii.  284f nazwalo trzech kandydatow.
// Jeden odpada przed pomiarem: conservativeParticleEnergy ZAWIERA
// darwinInteractionEnergy (electrodynamics.hpp:5434), wiec brakujacy
// potencjal Darwina nie jest wyjasnieniem.
//
// Zostaje ablacja dwuwymiarowa, bo kazdy kandydat reaguje inaczej:
//   blad calkowania        -> maleje z tolerancja
//   sektor dipolowy        -> zapada sie przy zgaszonych momentach
//   niezgodnosc sila/energia Darwina -> zostaje w obu
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
    std::printf("%10s %8s %15s %15s %11s\n","tolerancja","momenty",
        "dE_mech [J]","E_wyprom. [J]","iloraz");
    for(double tol:{1.0e-8,1.0e-10,1.0e-12}){
        for(int ablate=0;ablate<2;++ablate){
            const double scale=ablate?1.0e-6:1.0;
            const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*scale);
            const State start=osculatingPeriapsisState(el,kk,m1,
                m1*(secondMagneticMoment/firstMagneticMoment),
                Vec3{0,0,1},Vec3{1,0,0},0.0);
            ClassicalTrajectoryEngine::Accuracy acc;
            acc.relativeTolerance=tol; acc.maximumDepth=26;
            acc.reactionModel=
                ChargeRadiationReactionModel::stochasticElectricDipole;
            acc.computeOutwardFlux=true;
            acc.farFieldSampling.directionCount=194;
            acc.farFieldSampling.controlRadius=1.0e6*bohrRadius;
            acc.useRetardedExternalForces=false;   // sila zachowawcza
            ClassicalTrajectoryEngine engine(start,acc);
            State st=start; bool ok=true;
            for(int i=0;i<200&&ok;++i) ok=engine.advance(st,period/400);
            if(!ok){ std::printf("%10.0e %8s przerwane\n",tol,
                ablate?"off":"ON"); continue; }
            const double weightFirst=retardedDipoleSectorWeight(st);
            const double mechanicalFirst=
                conservativeParticleEnergy(st,weightFirst);
            const double radiatedFirst=st.radiatedEnergy;
            for(int i=0;i<400&&ok;++i) ok=engine.advance(st,period/400);
            if(!ok){ std::printf("%10.0e %8s przerwane\n",tol,
                ablate?"off":"ON"); continue; }
            const double weightLast=retardedDipoleSectorWeight(st);
            const double deltaMechanical=
                conservativeParticleEnergy(st,weightLast)-mechanicalFirst;
            const double radiated=st.radiatedEnergy-radiatedFirst;
            std::printf("%10.0e %8s %15.6e %15.6e %11.4f\n",tol,
                ablate?"off":"ON",deltaMechanical,radiated,
                radiated!=0.0?-deltaMechanical/radiated:0.0);
        }
    }
    std::printf("# maleje z tolerancja = calkowanie; zapada z momentami"
        " = dipol; stale = Darwin\n");
}
