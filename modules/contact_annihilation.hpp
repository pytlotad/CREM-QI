#pragma once

// Statistical experiment 6: contact annihilation from near-1s starts with
// quantized spins (audits 314-317).
//
// WHAT IT IS.  The para/ortho lifetime ratio from the one place in the model
// where the two channels separate with the physical sign: the smeared Fermi
// contact term of the dipole field, which pulls parallel moments (para) in
// and pushes antiparallel ones (ortho) out at r <~ r* (audit 312b), combined
// with the model's own 2 gamma / 3 gamma selection rule read at the moment of
// contact.  Circular orbits never sample the contact term (they live at
// ~5-547 r*, audit 312c); near-radial orbits with periapsis just above r* do.
//
// THE TWO IMPORTS, named rather than hidden.
//   (1) Annihilation at the first entry into r <= r* (terminalSeparation).
//       The model has no annihilation rate of its own; this rule gives it an
//       event.
//   (2) The Ore-Powell suppression of the 3 gamma channel,
//       eps = 4 (pi^2 - 9) alpha / (9 pi) = 1/1113.9, a QED number.
// Everything else is the model: the orbit, the contact force, the moments at
// the instant of contact and the 2 gamma weight w = (|mu1+mu2|/2 mu)^2, whose
// complement is the net-spin (3 gamma) weight exactly (README, 2g/3g section).
//
// CONFIGURATION, fixed by the experiment and restored afterwards:
//   * quantized spins (cos = +-1), so ortho carries w = 0 EXACTLY at contact
//     -- an exact symmetry, not a fit (audit 317: 7 of 7 entries, w = 0.0000).
//     CREM_CONTACT_FREE_SPINS=1 keeps the free (conditioned) draw instead;
//   * conservative dynamics (radiation reaction disabled): a 1s state does not
//     radiate, and with reaction on the near-radial orbits collapse
//     classically within a fraction of an orbit (audit 301g), drowning the
//     contact mechanism;
//   * the near-1s start: Bohr level n (--level), L = J0 times circular
//     (--contact-j0, default 0.07: Kepler periapsis just above r*).
//
// WHAT THE TIMES MEAN.  Rule (1) annihilates on every entry, so absolute
// times are femtoseconds and carry no physical meaning; the observable is the
// RATIO tau_ortho / tau_para.  Per channel the estimate is
//   tau = 1 / ( nu * < w + (1 - w) eps > ),   nu = entries / exposure,
// i.e. the contact entry rate times the mean annihilation probability per
// entry; the trajectory stops at the first entry, so the weight at later
// entries is assumed equal to the first.
//
// MEASURED (audit 317, J0 = 0.07, 30 pairs): quantized, tau_o/tau_p ~ 2.1e4
// = 19.2 (contact geometry) x 1114 (selection rule); the real ratio is 1135.
// The geometric factor has no quantum counterpart -- |psi(0)|^2 does not
// depend on spin at leading order -- and is a classical overstatement: r* is
// exactly the radius at which the dipole energy equals the Coulomb energy.
// Free spins give ~7.5.  Read the output with both facts in mind.

#include "configuration_panel.hpp"
#include "crem_trajectory.hpp"
#include "physical_constants.hpp"
#include "sampling_utilities.hpp"
#include "simulation_interface.hpp"
#include "vector3.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

// Default near-1s start of experiment 6 (audit 316: the separation is
// present at J0 = 0.07 and absent at 0.08).
inline double gContactInitialAngularMomentum = 0.07;
// Observation window in orbits; a trajectory without contact entry inside it
// is right-censored and contributes its full window to the exposure.
inline double gContactOrbitWindow = 40.0;

inline double orePowellSuppression() {
    return 4.0*(pi*pi-9.0)*fineStructureConstant/(9.0*pi);
}

struct ContactAnnihilationEvent {
    bool reachedContact=false;   // first entry into r <= r* inside the window
    double exposureSeconds=0.0;  // entry time, or the full window if censored
    double exposureOrbits=0.0;
    double startWeightTwoPhoton=0.0;
    double weightTwoPhoton=0.0;  // at the entry (or at the window end)
    double minimumSeparationOverBarrier=0.0;
    bool valid=false;
};

struct ContactAnnihilationPair {
    ContactAnnihilationEvent para;
    ContactAnnihilationEvent ortho;
};

inline double contactTwoPhotonWeight(const Vec3& first,const Vec3& second) {
    const double scale=first.norm()+second.norm();
    if(!(scale>0.0)) return 0.0;
    const double ratio=(first+second).norm()/scale;
    return ratio*ratio;
}

inline double contactOrbitalPeriod() {
    const double reducedMass=firstMass*secondMass/(firstMass+secondMass);
    const double level=static_cast<double>(std::max(1,gInitialPrincipalLevel));
    const double a=level*level*pairBohrRadius(activePair);
    return 2.0*pi*std::sqrt(a*a*a*reducedMass/pairCoulombStrength);
}

// One trajectory, one channel.  Expects the experiment's globals to be set
// (see runContactAnnihilationExperiment).
inline ContactAnnihilationEvent estimateContactAnnihilation(
        std::uint64_t seed,int phenomenon) {
    ContactAnnihilationEvent event;
    const double period=contactOrbitalPeriod();
    SimulationOptions options;
    options.collectFrames=true;
    options.frameCount=2;
    options.radiatedEnergyBookkeeping=false;
    options.observationTime=(gContactOrbitWindow+0.05)*period;
    options.terminalSeparation=comptonBarrierRadius;   // import (1)
    const SimulationResult run=simulate(seed,phenomenon,options);
    if(run.frames.empty()||run.outcome==SimulationOutcome::NumericalFailure)
        return event;
    const Frame& start=run.frames.front();
    // At ReachedCutoff the last frame is the interpolated crossing state.
    const Frame& last=run.frames.back();
    event.reachedContact=run.outcome==SimulationOutcome::ReachedCutoff;
    event.exposureSeconds=event.reachedContact
        ?run.elapsedTime:std::min(run.elapsedTime,gContactOrbitWindow*period);
    event.exposureOrbits=event.exposureSeconds/period;
    event.startWeightTwoPhoton=
        contactTwoPhotonWeight(start.firstDipole,start.secondDipole);
    event.weightTwoPhoton=
        contactTwoPhotonWeight(last.firstDipole,last.secondDipole);
    event.minimumSeparationOverBarrier=
        run.minimumSeparation/comptonBarrierRadius;
    event.valid=std::isfinite(event.exposureSeconds)
        &&event.exposureSeconds>0.0;
    return event;
}

// Paired ensemble: trajectory i uses splitMix64(masterSeed + i) in BOTH
// channels, as runCremCollapseExperiment does.  Sets and restores the
// experiment's configuration around the run.
inline std::vector<ContactAnnihilationPair> runContactAnnihilationExperiment(
        std::uint64_t masterSeed,int runCount) {
    const bool savedQuantization=gSpinQuantization;
    const ChargeRadiationReactionModel savedReaction=gRadiationReactionModel;
    const double savedFraction=gInitialAngularMomentumFraction;
    gSpinQuantization=std::getenv("CREM_CONTACT_FREE_SPINS")==nullptr;
    gRadiationReactionModel=ChargeRadiationReactionModel::disabled;
    gInitialAngularMomentumFraction=gContactInitialAngularMomentum;

    std::vector<ContactAnnihilationPair> pairs(
        static_cast<size_t>(std::max(0,runCount)));
    std::atomic<int> nextIndex{0};
    std::atomic<int> completed{0};
    std::mutex outputMutex;
    const int workerCount=std::max(1,std::min(runCount,
        static_cast<int>(std::max(1u,std::thread::hardware_concurrency()))));
    std::cout<<"Running "<<runCount<<" paired contact-annihilation trajectories"
             <<" on "<<workerCount<<" worker"<<(workerCount==1?"":"s")
             <<" (J0 = "<<gContactInitialAngularMomentum<<", window "
             <<gContactOrbitWindow<<" orbits, spins "
             <<(gSpinQuantization?"quantized":"free")<<").\n"
             <<"  Radiation reaction is DISABLED for this experiment"
               " (conservative dynamics, whatever the header above says).\n";
    const auto worker=[&]() {
        while(true) {
            const int index=nextIndex.fetch_add(1);
            if(index>=runCount) break;
            const std::uint64_t trajectorySeed=
                splitMix64(masterSeed+static_cast<std::uint64_t>(index));
            ContactAnnihilationPair& pair=pairs[static_cast<size_t>(index)];
            pair.para=estimateContactAnnihilation(trajectorySeed,1);
            pair.ortho=estimateContactAnnihilation(trajectorySeed,2);
            const int done=completed.fetch_add(1)+1;
            if(done%5==0||done==runCount) {
                std::lock_guard<std::mutex> lock(outputMutex);
                std::cout<<"contact pairs: "<<done<<"/"<<runCount<<'\n';
            }
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(workerCount));
    for(int index=0;index<workerCount;++index) workers.emplace_back(worker);
    for(std::thread& thread:workers) thread.join();

    gSpinQuantization=savedQuantization;
    gRadiationReactionModel=savedReaction;
    gInitialAngularMomentumFraction=savedFraction;
    return pairs;
}

struct ContactChannelSummary {
    int trajectories=0;
    int entries=0;
    double exposureOrbits=0.0;
    double entryRatePerOrbit=0.0;
    double meanWeightAtEntry=0.0;
    double meanAnnihilationProbability=0.0;   // < w + (1-w) eps >
    double lifetimeOrbits=0.0;                 // 1 / (nu < p >)
};

inline ContactChannelSummary summarizeContactChannel(
        const std::vector<ContactAnnihilationPair>& pairs,bool para) {
    ContactChannelSummary summary;
    const double eps=orePowellSuppression();
    double weightSum=0.0,probabilitySum=0.0;
    for(const ContactAnnihilationPair& pair:pairs) {
        const ContactAnnihilationEvent& e=para?pair.para:pair.ortho;
        if(!e.valid) continue;
        ++summary.trajectories;
        summary.exposureOrbits+=e.exposureOrbits;
        if(e.reachedContact) {
            ++summary.entries;
            weightSum+=e.weightTwoPhoton;
            probabilitySum+=e.weightTwoPhoton+(1.0-e.weightTwoPhoton)*eps;
        }
    }
    if(summary.exposureOrbits>0.0)
        summary.entryRatePerOrbit=summary.entries/summary.exposureOrbits;
    if(summary.entries>0) {
        summary.meanWeightAtEntry=weightSum/summary.entries;
        summary.meanAnnihilationProbability=probabilitySum/summary.entries;
        summary.lifetimeOrbits=1.0/(summary.entryRatePerOrbit
            *summary.meanAnnihilationProbability);
    }
    return summary;
}

inline int reportContactAnnihilationExperiment(std::uint64_t masterSeed,
                                               int runCount) {
    const std::vector<ContactAnnihilationPair> pairs=
        runContactAnnihilationExperiment(masterSeed,runCount);
    const ContactChannelSummary para=summarizeContactChannel(pairs,true);
    const ContactChannelSummary ortho=summarizeContactChannel(pairs,false);
    const double period=contactOrbitalPeriod();
    const double eps=orePowellSuppression();
    std::cout<<std::setprecision(6)
             <<"\nContact annihilation (statistical experiment 6)\n"
             <<"  imports: annihilation at first entry into r <= r* = "
             <<comptonBarrierRadius<<" m; 3 gamma suppression eps = "
             <<eps<<" = 1/"<<1.0/eps<<" (Ore-Powell)\n"
             <<"  absolute times are NOT physical (every entry annihilates);"
                " read the ratio\n\n";
    const auto line=[&](const char* name,const ContactChannelSummary& s) {
        std::cout<<"  "<<name<<": contact entries "<<s.entries<<"/"
                 <<s.trajectories<<", entry rate "<<s.entryRatePerOrbit
                 <<" /orbit, <w_2gamma> at entry "<<s.meanWeightAtEntry
                 <<", <w+(1-w)eps> "<<s.meanAnnihilationProbability;
        if(s.lifetimeOrbits>0.0)
            std::cout<<", tau = "<<s.lifetimeOrbits<<" orbits = "
                     <<s.lifetimeOrbits*period<<" s";
        std::cout<<'\n';
    };
    line("para ",para);
    line("ortho",ortho);
    if(para.lifetimeOrbits>0.0&&ortho.entries>0) {
        const double geometric=para.entryRatePerOrbit/ortho.entryRatePerOrbit;
        const double ratio=ortho.lifetimeOrbits/para.lifetimeOrbits;
        // Poisson counting error on the two entry numbers only.
        const double relative=std::sqrt(1.0/para.entries+1.0/ortho.entries);
        std::cout<<"\n  tau_ortho/tau_para = "<<ratio<<"  (+/- "
                 <<100.0*relative<<"% from entry counts)\n"
                 <<"    contact geometry alone: "<<geometric
                 <<"  (no quantum counterpart: |psi(0)|^2 is spin-independent"
                    " at leading order)\n"
                 <<"    selection rule alone:   "
                 <<para.meanAnnihilationProbability
                    /ortho.meanAnnihilationProbability<<'\n'
                 <<"    measured:               1135 (142.04 ns / 125.14 ps)\n";
    } else if(ortho.entries==0&&ortho.exposureOrbits>0.0) {
        std::cout<<"\n  ortho never reached contact in "<<ortho.exposureOrbits
                 <<" orbits of exposure; tau_ortho/tau_para is only bounded"
                    " from below\n";
    }
    return 0;
}
