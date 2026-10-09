#pragma once

// LABORATORY OBSERVABLES OF EXPERIMENTS 1/2 (audit 390, stage (a) of the
// distribution redesign; README, "Przebudowa rozkładów").
//
// Each trajectory of experiments 1/2 is turned into the event a laboratory
// records: the annihilation time on the laboratory clock and the
// annihilation photons in the laboratory frame.  The time and the channel
// come from the same function experiment 6 uses (finalStateDecayFromEstimate:
// the engine's survival along the cascade, then Exp(1/Gamma) of the final
// state; P(2 gamma) = w / (w + (1-w) eps)).  The photons are those of the
// experiment-6 generator in the pair's rest frame -- energies (W/2) x
// fraction with W the invariant mass the trajectory ended with -- boosted to
// the laboratory with the source velocity (--ps-source) composed with the
// cascade recoil.  No detector response is applied.
//
// What is the model's and what is imported, per panel:
//   time spectrum, survival   model (QR rate, w, cascade) against the
//                             measured lifetime;
//   photon energy             2 gamma: W/2 is the model's (binding + spin
//                             coupling), the width is the source's Doppler;
//                             3 gamma: the Ore-Powell distribution itself is
//                             a QED import (the generator), so the panel
//                             checks the generator and the kinematics;
//   photon angles             2 gamma: acollinearity from the source and the
//                             recoil only; 3 gamma: Ore-Powell kinematics.

#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

struct LabAnnihilationSample {
    std::vector<double> times;              // laboratory clock, plot units
    std::vector<double> photonEnergyKeV;    // every annihilation photon, lab
    std::vector<double> acollinearityMrad;  // 2 gamma: pi - angle(k1, k2)
    std::vector<double> threePhotonAngleDeg;// 3 gamma: the three pair angles
    std::vector<double> halfInvariantKeV;   // W/2 per event
    std::vector<double> sourceBeta;         // |v_source|/c per event
    std::vector<double> threePhotonRestFraction; // 3 gamma: x = E/(W/2), pair frame
    int twoPhotonEvents=0, threePhotonEvents=0;
    int duringCascade=0, invalid=0;
    double maximumTwoPhotonEnergyErrorKeV=0.0; // rest source: |E - W/2|
};

// Ore-Powell single-photon spectrum, x = E / (W/2) in [0, 1] (unnormalized;
// the formula of the catalogue entry ore_powell_three_gamma).
inline double orePowellSpectrumShape(double x) {
    if(!(x>0.0)) return 0.0;
    if(x>=1.0) return 1.0;  // limit of the expression at x -> 1
    const double l=std::log(1.0-x);
    return x*(1.0-x)/((2.0-x)*(2.0-x))-2.0*(1.0-x)*(1.0-x)*l/std::pow(2.0-x,3)
        +(2.0-x)/x+2.0*(1.0-x)*l/(x*x);
}

inline double orePowellSpectrumNorm() {
    static const double norm=[]{
        const int n=200000; double s=0.0;
        for(int i=0;i<n;++i) s+=orePowellSpectrumShape((i+0.5)/n);
        return s/n;
    }();
    return norm;
}

// Kolmogorov-Smirnov distance of photon energy fractions x = E/(W/2) against
// the Ore-Powell spectrum, and its asymptotic p value.
inline ExponentialFitCheck orePowellKolmogorovSmirnov(std::vector<double> x) {
    ExponentialFitCheck out;
    if(x.empty()) return out;
    static const std::vector<double> cdf=[]{
        const int n=20000; std::vector<double> table(n+1,0.0);
        for(int i=0;i<n;++i) table[i+1]=table[i]+orePowellSpectrumShape((i+0.5)/n);
        for(double& v: table) v/=table[n];
        return table;
    }();
    const auto F=[&](double v){
        const double t=std::clamp(v,0.0,1.0)*20000.0; const int i=std::min(19999,int(t));
        return cdf[i]+(cdf[i+1]-cdf[i])*(t-i); };
    std::sort(x.begin(),x.end());
    const double n=static_cast<double>(x.size());
    for(size_t i=0;i<x.size();++i) {
        const double f=F(x[i]);
        out.d=std::max(out.d,std::max(f-i/n,(i+1)/n-f));
    }
    const double lambda=(std::sqrt(n)+0.12+0.11/std::sqrt(n))*out.d;
    double p=0.0;
    for(int k=1;k<=100;++k)
        p+=2.0*((k%2)?1.0:-1.0)*std::exp(-2.0*k*k*lambda*lambda);
    out.p=std::clamp(p,0.0,1.0);
    return out;
}

inline LabAnnihilationSample collectLabAnnihilation(
        const std::vector<CremCollapseEstimate>& estimates,
        std::uint64_t masterSeed,int phenomenon,double timeScale) {
    LabAnnihilationSample sample;
    for(size_t index=0;index<estimates.size();++index) {
        const CremCollapseEstimate& estimate=estimates[index];
        // The same per-trajectory seed runCremCollapseExperiment used.
        const std::uint64_t seed=
            splitMix64(masterSeed+static_cast<std::uint64_t>(index));
        const FinalStateDecay decay=
            finalStateDecayFromEstimate(estimate,seed,phenomenon);
        if(!decay.valid||!std::isfinite(estimate.annihilationInvariantEnergy)) {
            ++sample.invalid;
            continue;
        }
        // Own stream for the source, so --ps-source never changes the decay.
        std::uint64_t stream=splitMix64(seed^0x50c0ce5ULL);
        const Vec3 sourceVelocity=drawPsSourceVelocity(gPsSource,stream);
        const Vec3 velocity=composePsVelocities(sourceVelocity,
            estimate.centreOfMassVelocityAtStop);
        const double sourceGamma=1.0/std::sqrt(
            1.0-sourceVelocity.squaredNorm()/(c*c));
        // Laboratory clock: the cascade on the engine's lab clock (recoil
        // included), the rest of the life at the recoil speed, all dilated
        // by the source motion.
        const double recoilGamma=1.0/std::sqrt(
            1.0-estimate.centreOfMassVelocityAtStop.squaredNorm()/(c*c));
        const double cascadeLabFactor=estimate.lifetimeSeconds>0.0
            ?estimate.lifetimeSecondsLab/estimate.lifetimeSeconds:1.0;
        double restClock=decay.duringCascade
            ?decay.decaySeconds*cascadeLabFactor
            :decay.cascadeSeconds*cascadeLabFactor
             +(decay.decaySeconds-decay.cascadeSeconds)*recoilGamma;
        sample.times.push_back(sourceGamma*restClock*timeScale);
        if(decay.duringCascade) ++sample.duringCascade;
        sample.sourceBeta.push_back(sourceVelocity.norm()/c);
        const double halfW=0.5*estimate.annihilationInvariantEnergy;
        const double keV=1.0e3*eCharge;
        sample.halfInvariantKeV.push_back(halfW/keV);
        const ContactAnnihilationPhotons& photons=decay.photons;
        std::vector<LabPhoton> lab;
        for(int i=0;i<photons.count;++i)
            lab.push_back(boostPhotonToLab(halfW*photons.fraction[i],
                photons.direction[i],velocity));
        for(const LabPhoton& photon: lab)
            sample.photonEnergyKeV.push_back(photon.energyJoules/keV);
        if(photons.count==2) {
            ++sample.twoPhotonEvents;
            const double cosine=std::clamp(
                dot(lab[0].direction,lab[1].direction),-1.0,1.0);
            sample.acollinearityMrad.push_back(1.0e3*(pi-std::acos(cosine)));
            if(!(velocity.norm()>0.0))
                for(const LabPhoton& photon: lab)
                    sample.maximumTwoPhotonEnergyErrorKeV=std::max(
                        sample.maximumTwoPhotonEnergyErrorKeV,
                        std::abs(photon.energyJoules-halfW)/keV);
        } else if(photons.count==3) {
            ++sample.threePhotonEvents;
            for(int i=0;i<3;++i)
                sample.threePhotonRestFraction.push_back(photons.fraction[i]);
            for(int i=0;i<3;++i) for(int j=i+1;j<3;++j)
                sample.threePhotonAngleDeg.push_back(180.0/pi*std::acos(
                    std::clamp(dot(lab[i].direction,lab[j].direction),
                               -1.0,1.0)));
        }
    }
    return sample;
}

// ROOT objects drawn on the pads must live until the PDFs are written.
struct LabPlotKeepAlive {
    std::vector<std::unique_ptr<TObject>> objects;
    template<class T> T* keep(T* object) {
        objects.emplace_back(object);
        return object;
    }
};

inline TPaveText* labNote(LabPlotKeepAlive& keep,double x1,double y1,
                          double x2,double y2,
                          const std::vector<std::string>& lines) {
    TPaveText* box=keep.keep(new TPaveText(x1,y1,x2,y2,"NDC"));
    box->SetFillColor(kWhite);
    box->SetBorderSize(1);
    box->SetTextAlign(12);
    box->SetTextSize(0.032);
    for(const std::string& line: lines) box->AddText(line.c_str());
    box->Draw();
    return box;
}

// Pad 1: annihilation-time spectrum on the laboratory clock, log scale, with
// the maximum-likelihood exponential and the measured lifetime.
inline void drawAnnihilationTimeSpectrum(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const LabAnnihilationSample& sample,double measuredTau,
        double measuredTauError,const char* unit,const std::string& reference) {
    pad->cd();
    pad->SetGrid();
    pad->SetLogy();
    const size_t n=sample.times.size();
    double mean=0.0;
    for(double t: sample.times) mean+=t;
    mean=n?mean/n:0.0;
    const double upper=std::max({8.0*measuredTau,8.0*mean,
        n?*std::max_element(sample.times.begin(),sample.times.end()):1.0});
    const int bins=std::clamp(static_cast<int>(std::sqrt(double(n))*2),10,80);
    TH1D* h=keep.keep(new TH1D("annihilation_time_spectrum",
        (std::string("Annihilation-time spectrum (laboratory clock);t [")+unit
         +"];events / bin").c_str(),bins,0.0,upper));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    for(double t: sample.times) h->Fill(t);
    h->SetLineColor(plot_style::crem());
    h->SetLineWidth(2);
    h->SetMinimum(0.5);
    h->Draw("HIST");
    const double width=upper/bins;
    TF1* fit=keep.keep(new TF1("annihilation_time_mle","[0]*exp(-x/[1])",0.0,upper));
    fit->SetParameters(n*width/std::max(mean,1e-300),mean);
    fit->SetLineColor(plot_style::crem());
    fit->SetLineStyle(2);
    fit->Draw("SAME");
    TF1* measured=keep.keep(new TF1("annihilation_time_measured","[0]*exp(-x/[1])",0.0,upper));
    measured->SetParameters(n*width/measuredTau,measuredTau);
    measured->SetLineColor(plot_style::experimental());
    measured->SetLineWidth(3);
    measured->SetLineStyle(3);
    measured->Draw("SAME");
    char buffer[200];
    std::vector<std::string> lines;
    std::snprintf(buffer,sizeof buffer,"model: tau (MLE) = %.4g #pm %.2g %s (N = %zu)",
                  mean,n?mean/std::sqrt(double(n)):0.0,unit,n);
    lines.push_back(buffer);
    std::snprintf(buffer,sizeof buffer,"measured: %.5g #pm %.2g %s",measuredTau,
                  measuredTauError,unit);
    lines.push_back(buffer);
    lines.push_back(reference);
    labNote(keep,0.42,0.70,0.98,0.90,lines);
}

// Pad 2: survival 1 - F(t) of the laboratory decay times against the measured
// exponential, with the Kolmogorov-Smirnov distance to it.
inline void drawAnnihilationSurvival(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const LabAnnihilationSample& sample,double measuredTau,const char* unit) {
    pad->cd();
    pad->SetGrid();
    pad->SetLogy();
    std::vector<double> t=sample.times;
    std::sort(t.begin(),t.end());
    const size_t n=t.size();
    std::vector<double> x,y;
    for(size_t i=0;i<n;++i) {
        x.push_back(t[i]); y.push_back(1.0-double(i)/n);
        x.push_back(t[i]); y.push_back(std::max(1.0-double(i+1)/n,0.5/n));
    }
    const double upper=n?std::max(t.back(),5.0*measuredTau):5.0*measuredTau;
    TGraph* g=keep.keep(new TGraph(static_cast<int>(x.size()),x.data(),y.data()));
    g->SetTitle((std::string("Survival of the pair (laboratory clock);t [")+unit
                 +"];fraction not yet annihilated").c_str());
    g->SetLineColor(plot_style::crem());
    g->SetLineWidth(2);
    g->Draw("AL");
    g->GetXaxis()->SetLimits(0.0,upper);
    g->SetMinimum(0.3/std::max<size_t>(n,1));
    g->SetMaximum(1.2);
    TF1* measured=keep.keep(new TF1("survival_measured","exp(-x/[0])",0.0,upper));
    measured->SetParameter(0,measuredTau);
    measured->SetLineColor(plot_style::experimental());
    measured->SetLineWidth(3);
    measured->SetLineStyle(3);
    measured->Draw("SAME");
    const ExponentialFitCheck ks=exponentialKolmogorovSmirnov(sample.times,measuredTau);
    char buffer[200];
    std::snprintf(buffer,sizeof buffer,"KS vs measured Exp(tau): D = %.3f, p = %.3g",
                  ks.d,ks.p);
    labNote(keep,0.42,0.80,0.98,0.90,{buffer});
}

// Pad 3: laboratory energy of the annihilation photons.
inline void drawAnnihilationPhotonEnergy(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const LabAnnihilationSample& sample,bool para) {
    pad->cd();
    pad->SetGrid();
    const double electronRestKeV=electronMass*c*c/(1.0e3*eCharge);
    double meanHalfW=0.0;
    for(double v: sample.halfInvariantKeV) meanHalfW+=v;
    meanHalfW=sample.halfInvariantKeV.empty()?electronRestKeV
        :meanHalfW/sample.halfInvariantKeV.size();
    if(para) {
        // A line: zoom on its width (Doppler) or on the binding shift.
        double spread=0.0;
        for(double e: sample.photonEnergyKeV)
            spread=std::max(spread,std::abs(e-meanHalfW));
        const double half=std::max({1.5*spread,2.0*(electronRestKeV-meanHalfW),1.0e-3});
        TH1D* h=keep.keep(new TH1D("annihilation_photon_energy",
            "Annihilation photon energy, laboratory (two-photon line);E_{#gamma} [keV];photons / bin",
            80,meanHalfW-half,meanHalfW+half));
        h->SetDirectory(nullptr);
    h->SetStats(false);
        for(double e: sample.photonEnergyKeV) h->Fill(e);
        h->SetLineColor(plot_style::crem());
        h->SetLineWidth(2);
        h->Draw("HIST");
        char buffer[200];
        std::snprintf(buffer,sizeof buffer,"W/2 = %.6f keV (model: binding + spin coupling)",meanHalfW);
        std::vector<std::string> lines{buffer};
        std::snprintf(buffer,sizeof buffer,"m_{e}c^{2} = %.6f keV",electronRestKeV);
        lines.push_back(buffer);
        lines.push_back("source: "+describePsSource(gPsSource)+"; no detector response");
        labNote(keep,0.14,0.75,0.98,0.90,lines);
        return;
    }
    TH1D* h=keep.keep(new TH1D("annihilation_photon_energy",
        "Annihilation photon energy, laboratory (three photons);E_{#gamma} [keV];photons / bin",
        60,0.0,1.02*meanHalfW));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    for(double e: sample.photonEnergyKeV) h->Fill(e);
    h->SetLineColor(plot_style::crem());
    h->SetLineWidth(2);
    h->Draw("HIST");
    const double width=1.02*meanHalfW/60.0;
    const double total=static_cast<double>(sample.photonEnergyKeV.size());
    TF1* ore=keep.keep(new TF1("ore_powell_spectrum",
        [meanHalfW](double* x,double* p){
            return p[0]*orePowellSpectrumShape(x[0]/meanHalfW);},
        0.0,meanHalfW,1));
    ore->SetParameter(0,total*width/(meanHalfW*orePowellSpectrumNorm()));
    ore->SetNpx(400);
    ore->SetLineColor(plot_style::qed());
    ore->SetLineWidth(2);
    ore->SetLineStyle(2);
    ore->Draw("SAME");
    labNote(keep,0.14,0.75,0.70,0.90,{
        "Ore-Powell (LO QED): the generator's import",
        "measured continuum consistent with QED:",
        "Chang, Tang & Li, PLB 157 (1985); detector response needed"});
}

// Pad 4: angles between the annihilation photons in the laboratory.
inline void drawAnnihilationPhotonAngles(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const LabAnnihilationSample& sample,bool para) {
    pad->cd();
    pad->SetGrid();
    if(para) {
        double maximum=0.0;
        for(double a: sample.acollinearityMrad) maximum=std::max(maximum,a);
        TH1D* h=keep.keep(new TH1D("annihilation_acollinearity",
            "Two-photon acollinearity, laboratory;#pi - #theta_{12} [mrad];events / bin",
            60,0.0,std::max(1.2*maximum,1.0e-6)));
        h->SetDirectory(nullptr);
    h->SetStats(false);
        for(double a: sample.acollinearityMrad) h->Fill(a);
        h->SetLineColor(plot_style::crem());
        h->SetLineWidth(2);
        h->Draw("HIST");
        labNote(keep,0.40,0.80,0.98,0.90,{
            "from the source motion and the cascade recoil only",
            "source: "+describePsSource(gPsSource)});
        return;
    }
    TH1D* h=keep.keep(new TH1D("annihilation_photon_angles",
        "Three-photon opening angles, laboratory;#theta_{ij} [deg];photon pairs / bin",
        60,0.0,180.0));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    for(double a: sample.threePhotonAngleDeg) h->Fill(a);
    h->SetLineColor(plot_style::crem());
    h->SetLineWidth(2);
    h->Draw("HIST");
    labNote(keep,0.14,0.80,0.75,0.90,{
        "coplanar at rest (sum = 360 deg); Ore-Powell kinematics"});
}

// Photon multiplicity of the decays: the model's 2 gamma / 3 gamma split.
inline void drawPhotonMultiplicity(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const LabAnnihilationSample& sample) {
    pad->cd();
    pad->SetGrid();
    TH1D* h=keep.keep(new TH1D("photon_multiplicity",
        "Annihilation photon multiplicity;number of photons;events",4,0.5,4.5));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    h->SetBinContent(2,sample.twoPhotonEvents);
    h->SetBinContent(3,sample.threePhotonEvents);
    h->SetFillColor(plot_style::crem());
    h->SetLineColor(plot_style::crem());
    h->Draw("HIST");
    const double total=sample.twoPhotonEvents+sample.threePhotonEvents;
    char buffer[200];
    std::snprintf(buffer,sizeof buffer,"2#gamma: %d, 3#gamma: %d (fraction 3#gamma %.4g)",
                  sample.twoPhotonEvents,sample.threePhotonEvents,
                  total>0?sample.threePhotonEvents/total:0.0);
    labNote(keep,0.30,0.78,0.98,0.90,{buffer,
        "P(2#gamma) = w/(w + (1-w)#varepsilon), #varepsilon = Ore-Powell (import)"});
}
