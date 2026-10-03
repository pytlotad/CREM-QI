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
//
// ANNIHILATION PHOTONS (audit 320).  Each contact entry also produces its
// final state, in the pair's rest frame, energies in units of W/2:
//   * the channel is drawn with P(2 gamma) = w / (w + (1-w) eps), the same
//     two rates the lifetime estimate adds;
//   * 2 gamma: back to back along an ISOTROPIC axis, x = 1 each -- a J = 0
//     state singles out no direction;
//   * 3 gamma: energies from the joint Ore-Powell density; the decay-plane
//     normal n is tied to the net spin S = mu1/g1 + mu2/g2 at the instant of
//     contact.  Under quantization ortho carries |S| = hbar along S-hat, the
//     model's version of the m = +1 state about S-hat, and tree-level QED
//     (audit 320a, computed from the six diagrams, the spin sum reproducing
//     Ore-Powell to 1e-15) gives for m = +-1
//         dN/dcos(theta_n) ~ 1 - cos^2(theta_n)/3,   theta_n = angle(n, S),
//     POINTWISE on the Dalitz plane (normal-to-in-plane ratio exactly 2), so
//     the plane orientation is independent of the energies.  This law is the
//     third import.  NOT carried: the orientation of the photon triangle
//     inside its plane, which QED correlates with the in-plane part of S
//     (anisotropy up to ~1); it is drawn uniformly here.  The law is even in
//     cos(theta_n), so the CPT-odd correlation S.(k1 x k2) averages to zero,
//     as in QED.

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

struct ContactAnnihilationPhotons {
    int count=0;                    // 2 or 3; 0 when there was no contact
    Vec3 direction[3];              // unit momenta
    double fraction[3]={0.0,0.0,0.0};   // E_i / (W/2); sum = 2
    // 3 gamma only: cos of the angle between the plane normal and S-hat, and
    // the J-PET CPT observable S-hat.(k1 x k2)/|k1 x k2| with |k1|>|k2|>|k3|.
    double normalSpinCosine=0.0;
    double orderedSpinCorrelation=0.0;
};

struct ContactAnnihilationEvent {
    bool reachedContact=false;   // first entry into r <= r* inside the window
    double exposureSeconds=0.0;  // entry time, or the full window if censored
    double exposureOrbits=0.0;
    double startWeightTwoPhoton=0.0;
    double weightTwoPhoton=0.0;  // at the entry (or at the window end)
    double minimumSeparationOverBarrier=0.0;
    Vec3 spinAtContactOverHbar;  // S1 + S2 = mu1/g1 + mu2/g2, in hbar
    ContactAnnihilationPhotons photons;
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

inline Vec3 contactNetSpinOverHbar(const Vec3& first,const Vec3& second) {
    const double g1=firstGyromagneticRatioOf(), g2=secondGyromagneticRatioOf();
    return ((g1!=0.0?first*(1.0/g1):Vec3{})+(g2!=0.0?second*(1.0/g2):Vec3{}))
           *(1.0/hbar);
}

inline Vec3 contactIsotropicDirection(std::uint64_t& stream) {
    const double c=2.0*drawUniformUnit(stream)-1.0;
    const double phi=2.0*pi*drawUniformUnit(stream);
    const double s=std::sqrt(std::max(0.0,1.0-c*c));
    return Vec3{s*std::cos(phi),s*std::sin(phi),c};
}

// Any unit vector perpendicular to a unit axis.
inline Vec3 contactPerpendicular(const Vec3& axis) {
    const Vec3 trial=std::abs(axis.z)<0.9?Vec3{0.0,0.0,1.0}:Vec3{1.0,0.0,0.0};
    const Vec3 p=cross(axis,trial);
    return p*(1.0/p.norm());
}

// Final state of one contact annihilation (see the header, audit 320).
inline ContactAnnihilationPhotons drawContactAnnihilationPhotons(
        double weightTwoPhoton,const Vec3& spinOverHbar,std::uint64_t& stream) {
    ContactAnnihilationPhotons out;
    const double eps=orePowellSuppression();
    const double w=std::clamp(weightTwoPhoton,0.0,1.0);
    const double twoPhoton=w/(w+(1.0-w)*eps);
    if(drawUniformUnit(stream)<twoPhoton) {
        const Vec3 axis=contactIsotropicDirection(stream);
        out.count=2;
        out.direction[0]=axis; out.direction[1]=-axis;
        out.fraction[0]=out.fraction[1]=1.0;
        return out;
    }
    out.count=3;
    drawOrePowellEnergyFractions(stream,out.fraction);
    const double spin=spinOverHbar.norm();
    Vec3 normal;
    if(spin>1e-12) {
        // n about S-hat from 1 - c^2/3 (maximum 1 at c = 0) by rejection.
        const Vec3 sHat=spinOverHbar*(1.0/spin);
        double c=0.0;
        do { c=2.0*drawUniformUnit(stream)-1.0; }
        while(drawUniformUnit(stream)>1.0-c*c/3.0);
        const double phi=2.0*pi*drawUniformUnit(stream);
        const Vec3 e1=contactPerpendicular(sHat), e2=cross(sHat,e1);
        const double s=std::sqrt(std::max(0.0,1.0-c*c));
        normal=sHat*c+(e1*std::cos(phi)+e2*std::sin(phi))*s;
    } else {
        normal=contactIsotropicDirection(stream);   // no spin, no axis
    }
    // Close the triangle in the plane: k1 at a uniform angle psi, k2 turned
    // by the opening angle theta12 about n, k3 = -(k1 + k2).
    const double x1=out.fraction[0], x2=out.fraction[1], x3=out.fraction[2];
    const double c12=std::clamp((x3*x3-x1*x1-x2*x2)/(2.0*x1*x2),-1.0,1.0);
    const double s12=std::sqrt(std::max(0.0,1.0-c12*c12));
    const double psi=2.0*pi*drawUniformUnit(stream);
    const Vec3 u=contactPerpendicular(normal), v=cross(normal,u);
    const Vec3 d1=u*std::cos(psi)+v*std::sin(psi);
    const Vec3 d2=d1*c12+cross(normal,d1)*s12;
    const Vec3 k3=-(d1*x1+d2*x2);
    out.direction[0]=d1; out.direction[1]=d2; out.direction[2]=k3*(1.0/k3.norm());
    if(spin>1e-12) {
        const Vec3 sHat=spinOverHbar*(1.0/spin);
        out.normalSpinCosine=dot(normal,sHat);
        int order[3]={0,1,2};
        std::sort(order,order+3,[&](int a,int b){
            return out.fraction[a]>out.fraction[b];});
        const Vec3 n12=cross(out.direction[order[0]],out.direction[order[1]]);
        out.orderedSpinCorrelation=dot(sHat,n12)/n12.norm();
    }
    return out;
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
    event.spinAtContactOverHbar=
        contactNetSpinOverHbar(last.firstDipole,last.secondDipole);
    if(event.reachedContact) {
        std::uint64_t stream=splitMix64(seed^0xa77a1e5u)
            +static_cast<std::uint64_t>(phenomenon);
        event.photons=drawContactAnnihilationPhotons(
            event.weightTwoPhoton,event.spinAtContactOverHbar,stream);
    }
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

struct ContactPhotonClosure { double momentum=0.0, energy=0.0; };

inline ContactPhotonClosure contactPhotonClosure(const ContactAnnihilationPhotons& p) {
    Vec3 total; double energy=0.0;
    for(int i=0;i<p.count;++i) { total+=p.direction[i]*p.fraction[i]; energy+=p.fraction[i]; }
    return {total.norm()/2.0,std::abs(energy-2.0)/2.0};   // relative to W/c, W
}

// Per-channel photon report, plus a self-test of the 3 gamma generator on a
// fixed spin axis where the expectations are exact: <P2(cos theta_n)> =
// -1/20 for 1 - c^2/3, <S.(k1 x k2)> = 0, <x> = 2/3.
inline void reportContactAnnihilationPhotons(
        const std::vector<ContactAnnihilationPair>& pairs,std::uint64_t masterSeed) {
    std::cout<<"\n  Annihilation photons (rest frame; 2 gamma axis isotropic;"
               " 3 gamma plane normal ~ 1 - cos^2/3 about S at contact, audit 320)\n";
    for(int channel=0;channel<2;++channel) {
        int two=0,three=0; double worstP=0.0,worstE=0.0,spin=0.0,p2=0.0,cpt=0.0;
        for(const ContactAnnihilationPair& pair:pairs) {
            const ContactAnnihilationEvent& e=channel==0?pair.para:pair.ortho;
            if(!e.valid||!e.reachedContact) continue;
            const ContactPhotonClosure c=contactPhotonClosure(e.photons);
            worstP=std::max(worstP,c.momentum); worstE=std::max(worstE,c.energy);
            spin+=e.spinAtContactOverHbar.norm();
            if(e.photons.count==2) ++two;
            if(e.photons.count==3) {
                ++three;
                const double x=e.photons.normalSpinCosine;
                p2+=0.5*(3.0*x*x-1.0); cpt+=e.photons.orderedSpinCorrelation;
            }
        }
        std::cout<<"  "<<(channel==0?"para ":"ortho")<<": 2 gamma "<<two
                 <<", 3 gamma "<<three<<", <|S1+S2|>/hbar at contact "
                 <<(two+three>0?spin/(two+three):0.0);
        if(three>0)
            std::cout<<", <P2(n.S)> "<<p2/three<<" (QED -0.05), <S.(k1 x k2)> "
                     <<cpt/three<<" (QED 0)";
        std::cout<<"\n         max |sum k|c/W "<<worstP<<", max |sum E - W|/W "<<worstE<<'\n';
    }
    std::uint64_t stream=splitMix64(masterSeed^0x5e1f7e57u);
    const int draws=200000;
    double p2=0.0,cpt=0.0,x=0.0,worst=0.0;
    for(int i=0;i<draws;++i) {
        const ContactAnnihilationPhotons p=
            drawContactAnnihilationPhotons(0.0,Vec3{0.0,0.0,1.0},stream);
        const double c=p.normalSpinCosine;
        p2+=0.5*(3.0*c*c-1.0); cpt+=p.orderedSpinCorrelation; x+=p.fraction[0];
        worst=std::max(worst,contactPhotonClosure(p).momentum);
    }
    std::cout<<"  generator self-test ("<<draws<<" draws, S = z): <P2> "<<p2/draws
             <<" (exact -0.05, sd 0.001), <S.(k1 x k2)> "<<cpt/draws
             <<" (0), <x1> "<<x/draws<<" (2/3), max |sum k|c/W "<<worst<<'\n';
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
    reportContactAnnihilationPhotons(pairs,masterSeed);
    return 0;
}
