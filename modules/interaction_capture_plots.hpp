#pragma once

// CAPTURE OBSERVABLES OF EXPERIMENT 5 (audit 393, stage (d) of the
// distribution redesign; README, "Przebudowa rozkładów").
//
// Experiment 5 collides free e+ e- and classifies the outcome; a capture is
// the classical counterpart of radiative recombination, e+ + e- -> Ps + photon.
// No laboratory measures that process for free e+ on free e-, so the panels
// here are compared with THEORY, labelled as such:
//   * Kramers' semiclassical radiative-recombination cross section,
//         sigma_n = (32 pi / 3 sqrt 3) alpha^3 a^2 eta^4 / (n (eta^2 + n^2)),
//         eta^2 = Ry/E (Sommerfeld parameter),
//     as quoted by Kotelnikov & Milstein, Phys. Scr. (2019),
//     doi:10.1088/1402-4896/ab060a, Eq. (23) (Kramers 1923); below the exact
//     (Stobbe) result by 20 % at n = 1, 12 % at n = 2, 7.5 % at n = 4.  For
//     a pair the Bohr radius and Rydberg are those of the reduced mass
//     (a_Ps = 2 a_0, Ry_Ps = 6.8 eV): the radiating dipole is e r_rel, as in
//     the model.
//   * the statistical spin weights of the formed state, singlet : triplet =
//     1 : 3, against the model's para/ortho classification.
//
// Cross section from a Gaussian beam.  The transverse offset is
// (b_y, b_z) ~ N_2(0, sigma^2 I) truncated at R = R_match/2, i.e. b has the
// density P(b) = (b/sigma^2) exp(-b^2/(2 sigma^2)) / Z, Z = 1 - exp(-R^2/(2
// sigma^2)).  sigma_cap = int p_cap(b) 2 pi b db is estimated without
// assuming small b (Horvitz-Thompson): each event weighs
// 2 pi b / P(b) = 2 pi sigma^2 Z exp(b^2/(2 sigma^2)) and sigma_cap is the
// mean of (captured ? weight : 0); valid while captures stay inside R.

#include <TH1D.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>

inline double kramersPartialCrossSection(int n,double etaSquared,double bohr) {
    return 32.0*pi/(3.0*std::sqrt(3.0))*std::pow(fineStructureConstant,3)
        *bohr*bohr*etaSquared*etaSquared/(n*(etaSquared+double(n)*n));
}

inline double kramersTotalCrossSection(double etaSquared,double bohr) {
    double sum=0.0;
    for(int n=1;n<=200000;++n) {
        const double term=kramersPartialCrossSection(n,etaSquared,bohr);
        sum+=term;
        if(n>100&&term<1e-12*sum) break;
    }
    return sum;
}

struct CaptureSample {
    std::vector<double> energyEv, captured;      // per event: E_CM, 1 if captured
    std::vector<double> nEffective, orbitalHbar; // captures only
    std::vector<double> capturedImpact;          // captures only, m
    std::vector<double> weight;                  // per event: 2 pi b / P(b), m^2
    int para=0, ortho=0;
};

inline CaptureSample collectCaptures(const std::vector<InteractionEvent>& events,
                                    const InteractionConfiguration& configuration) {
    CaptureSample s;
    const double sigma=configuration.impactParameterSigma;
    const double half=0.5*configuration.matchingRadius;
    const double z=1.0-std::exp(-half*half/(2.0*sigma*sigma));
    const double rydberg=pairCoulombStrength*pairCoulombStrength*pairReducedMass
        /(2.0*hbar*hbar);
    for(const InteractionEvent& e: events) {
        if(e.outcome==InteractionOutcome::Unresolved
           ||e.outcome==InteractionOutcome::NumericalFailure
           ||!std::isfinite(e.kineticEnergyEv)) continue;
        const bool bound=e.outcome==InteractionOutcome::ParaPositronium
                       ||e.outcome==InteractionOutcome::OrthoPositronium;
        s.energyEv.push_back(e.kineticEnergyEv);
        s.captured.push_back(bound?1.0:0.0);
        s.weight.push_back(2.0*pi*sigma*sigma*z*std::exp(
            e.impactParameter*e.impactParameter/(2.0*sigma*sigma)));
        if(!bound) continue;
        if(e.outcome==InteractionOutcome::ParaPositronium) ++s.para; else ++s.ortho;
        if(std::isfinite(e.finalRelativeEnergyEv)&&e.finalRelativeEnergyEv<0.0)
            s.nEffective.push_back(std::sqrt(rydberg/(-e.finalRelativeEnergyEv*eCharge)));
        if(std::isfinite(e.finalOrbitalAngularMomentumHbar))
            s.orbitalHbar.push_back(e.finalOrbitalAngularMomentumHbar);
        s.capturedImpact.push_back(e.impactParameter);
    }
    return s;
}

// Pad 1: capture cross section against E_CM with Kramers' total.
inline void drawCaptureCrossSection(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const CaptureSample& s,const InteractionConfiguration& configuration) {
    pad->cd();
    pad->SetGrid();
    pad->SetLogx();
    pad->SetLogy();
    double lo=1e300,hi=0.0;
    for(double e: s.energyEv) { lo=std::min(lo,e); hi=std::max(hi,e); }
    if(!(hi>lo)) { lo=0.1; hi=10.0; }
    const int bins=8;
    const double barnUnit=1.0e-28;
    std::vector<double> x,y,exl,exh,eyl,eyh;
    for(int b=0;b<bins;++b) {
        const double e0=lo*std::pow(hi/lo,double(b)/bins), e1=lo*std::pow(hi/lo,double(b+1)/bins);
        double n=0,sum=0,sum2=0;
        for(size_t i=0;i<s.energyEv.size();++i)
            if(s.energyEv[i]>=e0&&(s.energyEv[i]<e1||(b==bins-1&&s.energyEv[i]<=e1))) {
                const double v=s.captured[i]*s.weight[i];
                n+=1; sum+=v; sum2+=v*v; }
        if(n<5) continue;
        const double centre=std::sqrt(e0*e1), mean=sum/n;
        const double sd=std::sqrt(std::max(0.0,sum2/n-mean*mean)/n);
        // An empty bin is drawn at an upper bound of one capture of the
        // bin's mean weight rather than at zero on a log axis.
        double meanWeight=0.0;
        for(size_t i=0;i<s.energyEv.size();++i)
            if(s.energyEv[i]>=e0&&s.energyEv[i]<=e1) meanWeight+=s.weight[i];
        meanWeight/=n;
        const double value=mean>0.0?mean:meanWeight/n;
        x.push_back(centre); exl.push_back(centre-e0); exh.push_back(e1-centre);
        y.push_back(value/barnUnit);
        eyl.push_back((mean>0.0?std::min(sd,0.9*mean):0.9*value)/barnUnit);
        eyh.push_back((mean>0.0?sd:0.0)/barnUnit);
    }
    TGraphAsymmErrors* model=keep.keep(new TGraphAsymmErrors(static_cast<int>(x.size()),
        x.data(),y.data(),exl.data(),exh.data(),eyl.data(),eyh.data()));
    model->SetTitle("Capture cross section e^{+}e^{-} #rightarrow Ps (classical radiative capture);"
                    "E_{CM} [eV];#sigma_{capture} [barn]");
    model->SetMarkerStyle(20);
    model->SetMarkerColor(plot_style::crem());
    model->SetLineColor(plot_style::crem());
    std::vector<double> tx,ty;
    const double rydberg=pairCoulombStrength*pairCoulombStrength*pairReducedMass/(2.0*hbar*hbar);
    const double bohr=pairBohrRadius(activePair);
    for(int i=0;i<=60;++i) {
        const double e=0.5*lo*std::pow(4.0*hi/lo,i/60.0);
        tx.push_back(e);
        ty.push_back(kramersTotalCrossSection(rydberg/(e*eCharge),bohr)/barnUnit);
    }
    TGraph* kramers=keep.keep(new TGraph(static_cast<int>(tx.size()),tx.data(),ty.data()));
    kramers->SetLineColor(plot_style::theory());
    kramers->SetLineStyle(2);
    kramers->SetLineWidth(2);
    double ymin=*std::min_element(ty.begin(),ty.end()), ymax=*std::max_element(ty.begin(),ty.end());
    for(double v: y) { ymin=std::min(ymin,v); ymax=std::max(ymax,v); }
    model->SetMinimum(0.1*ymin);
    model->SetMaximum(3000.0*ymax);   // room for the note above the points
    model->Draw("AP");
    model->GetXaxis()->SetLimits(0.5*lo,2.0*hi);
    kramers->Draw("L SAME");
    double bmax=0.0;
    for(double b: s.capturedImpact) bmax=std::max(bmax,b);
    char buffer[200];
    std::snprintf(buffer,sizeof buffer,"captures %zu of %zu; max captured b %.3g pm; #sigma_{b} %.3g pm",
                  s.capturedImpact.size(),s.energyEv.size(),bmax*1e12,
                  configuration.impactParameterSigma*1e12);
    labNote(keep,0.14,0.70,0.98,0.90,{buffer,
        "dashed: Kramers, semiclassical theory (a_{Ps}, Ry_{Ps})",
        "Kotelnikov & Milstein, Phys. Scr. 94, 055403 (2019)",
        "no measurement of free e^{+} on free e^{-}"});
}

// Pad 2: effective principal number of the captured state against Kramers'
// partial cross sections at the same energies.
inline void drawCapturedLevels(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const CaptureSample& s,const std::vector<InteractionEvent>& events) {
    pad->cd();
    pad->SetGrid();
    double top=10.0;
    for(double n: s.nEffective) top=std::max(top,n);
    const int bins=std::min(60,static_cast<int>(std::ceil(top))+1);
    TH1D* h=keep.keep(new TH1D("captured_level_distribution",
        "Captured state: n_{eff} = (Ry/|E|)^{1/2};n_{eff};captures / bin",bins,0.5,bins+0.5));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    for(double n: s.nEffective) h->Fill(n);
    h->SetLineColor(plot_style::crem());
    h->SetLineWidth(2);
    // Kramers expectation: per captured event at its energy, the share of
    // level n in the total; summed and drawn as a histogram.
    TH1D* expect=keep.keep(new TH1D("captured_level_kramers","Kramers",bins,0.5,bins+0.5));
    expect->SetDirectory(nullptr);
    const double rydberg=pairCoulombStrength*pairCoulombStrength*pairReducedMass/(2.0*hbar*hbar);
    const double bohr=pairBohrRadius(activePair);
    for(const InteractionEvent& e: events) {
        if(e.outcome!=InteractionOutcome::ParaPositronium
           &&e.outcome!=InteractionOutcome::OrthoPositronium) continue;
        const double eta2=rydberg/(e.kineticEnergyEv*eCharge);
        const double total=kramersTotalCrossSection(eta2,bohr);
        for(int n=1;n<=bins;++n)
            expect->Fill(n,kramersPartialCrossSection(n,eta2,bohr)/total);
    }
    expect->SetLineColor(plot_style::theory());
    expect->SetLineStyle(2);
    expect->SetLineWidth(2);
    h->SetMaximum(1.4*std::max(h->GetMaximum(),expect->GetMaximum()));
    h->Draw("HIST");
    expect->Draw("HIST SAME");
    labNote(keep,0.35,0.80,0.98,0.90,{"dashed: Kramers #sigma_{n} shares (theory)"});
}

// Pad 3: orbital angular momentum of the captured state, with the Langer
// grid (l + 1/2) hbar marked.
inline void drawCapturedAngularMomentum(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const CaptureSample& s) {
    pad->cd();
    pad->SetGrid();
    double top=5.0;
    for(double l: s.orbitalHbar) top=std::max(top,l);
    TH1D* h=keep.keep(new TH1D("captured_angular_momentum",
        "Captured state: orbital L;L/#hbar;captures / bin",
        std::min(80,static_cast<int>(4*std::ceil(top))),0.0,std::ceil(top)));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    for(double l: s.orbitalHbar) h->Fill(l);
    h->SetLineColor(plot_style::crem());
    h->SetLineWidth(2);
    h->SetMaximum(1.4*std::max(1.0,h->GetMaximum()));
    h->Draw("HIST");
    labNote(keep,0.35,0.82,0.98,0.90,{"classical L is continuous; Langer grid (l+1/2)#hbar for reference"});
}

// Pad 4: para/ortho split against the statistical weights 1 : 3.
inline void drawSpinChannelSplit(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const CaptureSample& s) {
    pad->cd();
    pad->SetGrid();
    TH1D* h=keep.keep(new TH1D("capture_spin_channels",
        "Captured spin channel;;captures",2,0.0,2.0));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    h->GetXaxis()->SetBinLabel(1,"para (singlet)");
    h->GetXaxis()->SetBinLabel(2,"ortho (triplet)");
    h->SetBinContent(1,s.para);
    h->SetBinContent(2,s.ortho);
    h->SetFillColor(plot_style::crem());
    TH1D* weights=keep.keep(new TH1D("capture_spin_weights","1:3",2,0.0,2.0));
    weights->SetDirectory(nullptr);
    const double total=s.para+s.ortho;
    weights->SetBinContent(1,0.25*total);
    weights->SetBinContent(2,0.75*total);
    weights->SetLineColor(plot_style::theory());
    weights->SetLineStyle(2);
    weights->SetLineWidth(3);
    h->SetMaximum(1.5*std::max(1.0,std::max(h->GetMaximum(),weights->GetMaximum())));
    h->Draw("HIST");
    weights->Draw("HIST SAME");
    char buffer[160];
    std::snprintf(buffer,sizeof buffer,"model para:ortho = %d:%d; statistical weights 1:3 (dashed)",
                  s.para,s.ortho);
    labNote(keep,0.14,0.80,0.98,0.90,{buffer});
}
