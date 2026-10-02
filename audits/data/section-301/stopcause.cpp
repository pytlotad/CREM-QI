// Audit 292: czy produkcyjne orbity schodza pod bariere Comptona?
//
// 291f(ii) zostawilo to jako otwarte.  Czytanie kodu daje odpowiedz
// bogatsza, niz bylo pytanie: pole result.stopCause juz to rejestruje,
// a komentarz przy CremCollapseEstimate twierdzi dwie rzeczy:
//   (A) bariera konczy tylko 24 procent trajektorii, a 76 procent
//       konczy sie na stosunku okres/przelot swiatla <= 150, ktory
//       jest NUMERYCZNYM marginesem, nie skala fizyczna;
//   (B) jeden foton przenosi perycentrum o czynnik 10-60 (zmierzone
//       21,7 r* -> 0,35 r*), wiec trajektorie wychodza Z WNETRZA
//       obszaru, ktory limity wykluczaja, a stosunek siada na 36
//       przeciw progowi 150.
// Lekcja z 291 mowi tego nie przyjmowac na slowo.  Ta sonda mierzy
// oba twierdzenia wprost, czytajac instrumenty, ktore model juz ma.
//
// Zarejestrowane PRZED uruchomieniem, zeby bylo falsyfikowalne:
//   R292a: udzial ComptonBarrier bedzie mniejszoscia (< 50 procent).
//   R292b: terminalPeriapsisOverBarrier bedzie w wiekszosci
//          PONIZEJ 1, czyli para konczy pod bariera, a nie na niej --
//          bo (B) mowi o przeskoku, nie o dojsciu.
//   R292c: terminalPeriodToLightCrossing bedzie ponizej progu 150 dla
//          trajektorii konczonych RetardationLimit (to jest definicja)
//          i rowniez ponizej 150 dla czesci konczonych ComptonBarrier
//          (to NIE jest definicja i jest tresciwym przewidywaniem (B)).
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const int seeds=(argc>1)?std::atoi(argv[1]):12;
    const double budget=(argc>2)?std::atof(argv[2]):40.0;
    const int phenomenon=(argc>3)?std::atoi(argv[3]):1;   // 1 = para
    // DOMYSLNY gRadiationReactionModel to individualLandauLifshitz, czyli
    // model CIAGLY, w ktorym nie ma fotonow i oba przelaczniki z 288 sa
    // bezczynne -- dlatego pierwszy przebieg dal kolumny bit-identyczne.
    // Testowane twierdzenia (24 procent ComptonBarrier, 21,7 r* -> 0,35 r*
    // na jeden foton) dotycza KASKADY FOTONOWEJ, wiec tryb trzeba wybrac
    // jawnie.  CREM_CONTINUOUS przywraca ciagly dla kontroli.
    if(std::getenv("CREM_CONTINUOUS")==nullptr)
        gRadiationReactionModel=
            ChargeRadiationReactionModel::stochasticElectricDipole;
    // 295c: konfiguracja REJESTRU, zrekonstruowana w 295b przez
    // wykluczenie.  README (linia 8885) kaze podac jawnie
    // --ground-state-floor --bohr-photon-energy --emission poisson,
    // ale 295b pokazalo, ze punkt 547,5 r* z wlaczona podloga jest
    // dokladnie a_pair, czyli orbita n=1 -- wiec rejestr, ktory nie
    // raportuje ani jednej przyczyny GroundStateFloor, musi byc BEZ
    // podlogi.  Stad ten zestaw: poziom i podloga jak dzis, a emisja i
    // kwant jak wtedy.  Kazda flaga osobno, zeby dalo sie je rozdzielic.
    if(std::getenv("CREM_PROBE_POISSON")) gDeterministicEmission=false;
    if(std::getenv("CREM_PROBE_BOHR")) gBohrLevelPhotonEnergy=true;
    if(std::getenv("CREM_PROBE_FLOOR")) gGroundStateEmissionFloor=true;
    if(const char* l=std::getenv("CREM_PROBE_LEVEL"))
        gInitialPrincipalLevel=std::atoi(l);
    std::printf("# poziom=%d podloga=%d bohr=%d emisja=%s\n",
        gInitialPrincipalLevel,gGroundStateEmissionFloor?1:0,
        gBohrLevelPhotonEnergy?1:0,
        gDeterministicEmission?"deterministyczna":"poisson");
    std::printf("# model reakcji: %s\n",
        gRadiationReactionModel
            ==ChargeRadiationReactionModel::stochasticElectricDipole
        ?"stochasticElectricDipole":"individualLandauLifshitz");
    std::printf("# ziaren=%d budzet=%.0f s/ziarno zjawisko=%d\n",
        seeds,budget,phenomenon);
    std::printf("# bariera r*=%.6e m, L_kontakt=sqrt(mu k r*)=%.4f hbar\n",
        comptonBarrierRadius,
        std::sqrt(firstMass*secondMass/(firstMass+secondMass)
            *pairCoulombStrength*comptonBarrierRadius)
        /hbar);
    std::printf("%5s %16s %10s %10s %9s %8s %11s %11s\n","ziarno",
        "przyczyna","r_p/r*","okres/prz","L [hbar]","1-e^2",
        "zycie [s]","t_kontakt");
    std::printf("# nan w zyciu = przebieg ocenzurowany budzetem,"
        " nie wynik fizyczny\n");
    int counts[4]={0,0,0,0};
    std::vector<double> ratios,periapses;
    int belowBarrier=0,ratioBelowThreshold=0,barrierAndBelow=0;
    int contactReached=0;
    for(int i=0;i<seeds;++i){
        const CremCollapseEstimate e=
            estimateCremCollapse(static_cast<std::uint64_t>(
                (std::getenv("CREM_PROBE_SEED")
                    ?std::atoi(std::getenv("CREM_PROBE_SEED")):40)+i),
                phenomenon,budget);
        const char* name="?";
        switch(e.stopCause){
            case CollapseStopCause::None: name="None"; break;
            case CollapseStopCause::ComptonBarrier:
                name="ComptonBarrier"; break;
            case CollapseStopCause::RetardationLimit:
                name="RetardationLimit"; break;
            case CollapseStopCause::GroundStateFloor:
                name="GroundStateFloor"; break;
        }
        counts[static_cast<int>(e.stopCause)]++;
        // 296: liczniki odmow -- R296a przewiduje, ze po spadku n ponizej
        // 1/2 kazda proba emisji jest odmawiana przez ceiling.
        std::printf("# fotony=%lld odmowy ceiling=%llu kinematyka=%llu"
            " odrzut=%llu\n",e.emittedPhotonCount,e.refusedByCeiling,
            e.refusedByKinematics,e.refusedByRecoil);
        // 301: wynik kalibracji odroznia cenzure zegarowa/limitowa
        // (ObservationLimit) od awarii i od konca fizycznego (ReachedCutoff).
        std::printf("# outcome=%d (0=ReachedCutoff 1=ObservationLimit"
            " 2=NumericalFailure) calib=%.6e s\n",
            static_cast<int>(e.calibrationOutcome),e.calibrationSeconds);
        std::printf("%5d %16s %10.4f %10.2f %9.5f %8.4f %11.4e %11.4e\n",
            40+i,name,e.terminalPeriapsisOverBarrier,
            e.terminalPeriodToLightCrossing,e.terminalAngularMomentum,
            e.terminalKeplerConsistency,e.lifetimeSeconds,
            e.contactPassageAtBarrierSeconds);
        if(std::isfinite(e.contactPassageAtBarrierSeconds)) ++contactReached;
        if(std::isfinite(e.terminalPeriapsisOverBarrier)) {
            periapses.push_back(e.terminalPeriapsisOverBarrier);
            if(e.terminalPeriapsisOverBarrier<1.0) ++belowBarrier;
        }
        if(std::isfinite(e.terminalPeriodToLightCrossing)) {
            ratios.push_back(e.terminalPeriodToLightCrossing);
            if(e.terminalPeriodToLightCrossing<150.0)
                ++ratioBelowThreshold;
            if(e.stopCause==CollapseStopCause::ComptonBarrier
               &&e.terminalPeriodToLightCrossing<150.0)
                ++barrierAndBelow;
        }
    }
    std::printf("\n# R292a udzial ComptonBarrier: %d/%d"
        "  (RetardationLimit %d, GroundStateFloor %d, None %d)\n",
        counts[1],seeds,counts[2],counts[3],counts[0]);
    std::printf("# R292b konczy PONIZEJ bariery: %d z %d skonczonych\n",
        belowBarrier,static_cast<int>(periapses.size()));
    std::printf("# R292c stosunek < 150: %d z %d; z tego konczonych"
        " BARIERA: %d\n",ratioBelowThreshold,
        static_cast<int>(ratios.size()),barrierAndBelow);
    std::printf("# R292d L doszlo do kontaktu (L <= 0.0427 hbar):"
        " %d/%d\n",contactReached,seeds);
    if(!periapses.empty()){
        double lo=1e300,hi=-1e300;
        for(double v:periapses){lo=std::min(lo,v);hi=std::max(hi,v);}
        std::printf("# r_p/r* zakres: %.4f .. %.4f\n",lo,hi);
    }
    if(!ratios.empty()){
        double lo=1e300,hi=-1e300;
        for(double v:ratios){lo=std::min(lo,v);hi=std::max(hi,v);}
        std::printf("# okres/przelot zakres: %.2f .. %.2f"
            " (prog 150)\n",lo,hi);
    }
}
