#pragma once

// LABORATORY (FIXED-TARGET) VIEW OF EXPERIMENTS 3/4 (audit 392, stage (c) of
// the distribution redesign; README, "Przebudowa rozkładów").
//
// The beam experiments are computed in the centre-of-momentum frame.  A
// laboratory scatters a positron beam on electrons at rest, so the same
// events are shown here in that frame: the positron's laboratory angle and,
// for experiment 4, the recoil-electron energy, which is what the Bhabha
// measurements recorded.  Equal masses only (the frame relations below use
// m1 = m2).  The relations are exact relativistically:
//     sqrt(s) = 2 m c^2 + K_CM,  T_lab = s/(2 m c^2) - 2 m c^2,
//     tan(theta_lab) = tan(theta*/2) / gamma_cm,  gamma_cm = (E_lab + m c^2)/sqrt(s),
//     y = T_recoil / T_lab = (1 - cos theta*)/2.
// Each angular bin of the CM histogram maps one-to-one onto a laboratory bin,
// so the cross section per bin is the same and only the solid angle (or the
// width in y) changes.
//
// The reference curves are THEORY: Rutherford (for distinguishable e+ e- the
// exact non-relativistic quantum result), tree-level Bhabha
// (qedElasticDifferentialCrossSection), and for experiment 3 Dirac's
// annihilation-in-flight cross section.  No measurement overlaps the model's
// domain: the Bhabha measurements (Ashkin, Page & Woodward, Phys. Rev. 94,
// 357 (1954): Bhabha's formula verified within 10 % at 0.6-1.0 MeV) are at
// angles where the closest approach lies far below the Compton barrier, where
// the model declares itself invalid (audit 392).  They validate the Bhabha
// curve, not the model.

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>

inline constexpr double beamLabBarn=1.0e-28;   // m^2

struct BeamLabFrame {
    double restEnergy=0.0;       // m c^2 of one particle
    double labKinetic=0.0;       // T_lab of the beam particle
    double gammaCm=1.0;
    bool valid=false;
};

inline BeamLabFrame beamLabFrame(double centreOfMassKineticEnergy,double mass) {
    BeamLabFrame f;
    f.restEnergy=mass*c*c;
    const double sqrtS=2.0*f.restEnergy+centreOfMassKineticEnergy;
    const double s=sqrtS*sqrtS;
    const double labEnergy=s/(2.0*f.restEnergy)-f.restEnergy;  // total, beam particle
    f.labKinetic=labEnergy-f.restEnergy;
    f.gammaCm=(labEnergy+f.restEnergy)/sqrtS;
    f.valid=std::isfinite(f.labKinetic)&&f.labKinetic>0.0;
    return f;
}

inline double labAngleFromCm(double thetaCm,const BeamLabFrame& f) {
    return std::atan(std::tan(0.5*thetaCm)/f.gammaCm);
}

// Dirac's two-photon annihilation cross section of a positron of Lorentz
// factor gamma on an electron at rest (QED, lowest order); r_e the
// classical electron radius.  Low-energy limit: pi r_e^2 c/v.
inline double diracAnnihilationCrossSection(double gamma) {
    const double re=coulombConstant*eCharge*eCharge/(electronMass*c*c);
    if(!(gamma>1.0)) return std::numeric_limits<double>::infinity();
    const double root=std::sqrt(gamma*gamma-1.0);
    return pi*re*re/(gamma+1.0)*((gamma*gamma+4.0*gamma+1.0)/(gamma*gamma-1.0)
        *std::log(gamma+root)-(gamma+3.0)/root);
}

// Per-bin cross sections over the CM bin edges (radians): model counts,
// Rutherford (closed form) and Bhabha (Simpson in cos theta*).
struct BeamBinCrossSections {
    std::vector<double> edges;           // theta* edges, radians
    std::vector<double> model, modelError, rutherford, qed;  // m^2 per bin
};

inline BeamBinCrossSections beamBinCrossSections(const std::vector<double>& edgesDeg,
        const std::vector<double>& thetaCm,double sampledArea,int runCount,
        double coulombLength,double centreOfMassKineticEnergy,double m1,double m2) {
    BeamBinCrossSections b;
    for(double e: edgesDeg) b.edges.push_back(e*pi/180.0);
    const size_t bins=b.edges.size()-1;
    std::vector<double> counts(bins,0.0);
    for(double t: thetaCm) {
        auto it=std::upper_bound(b.edges.begin(),b.edges.end(),t);
        if(it==b.edges.begin()||it==b.edges.end()) continue;
        counts[static_cast<size_t>(it-b.edges.begin()-1)]+=1.0;
    }
    for(size_t i=0;i<bins;++i) {
        const double lo=b.edges[i], hi=b.edges[i+1];
        const double scale=sampledArea/runCount;
        b.model.push_back(counts[i]*scale);
        b.modelError.push_back(std::sqrt(counts[i]*(1.0-counts[i]/runCount))*scale);
        b.rutherford.push_back(pi*coulombLength*coulombLength
            *(1.0/std::pow(std::tan(0.5*lo),2)-1.0/std::pow(std::tan(0.5*hi),2)));
        const double cosLo=std::cos(lo), cosHi=std::cos(hi);
        double sum=0.0;
        constexpr int n=8;
        for(int k=0;k<=n;++k) {
            const double ct=cosHi+(cosLo-cosHi)*k/n;
            const double w=(k==0||k==n)?1.0:(k%2?4.0:2.0);
            sum+=w*qedElasticDifferentialCrossSection(centreOfMassKineticEnergy,
                std::acos(std::clamp(ct,-1.0,1.0)),m1,m2);
        }
        b.qed.push_back(2.0*pi*(cosLo-cosHi)*sum/(3.0*n));
    }
    return b;
}

// Pad: d sigma / d Omega_lab of the scattered positron.
inline void drawLabDifferentialCrossSection(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const BeamBinCrossSections& b,const BeamLabFrame& f) {
    pad->cd();
    pad->SetGrid();
    pad->SetLogy();
    std::vector<double> labEdges;
    for(double t: b.edges) labEdges.push_back(labAngleFromCm(t,f)*180.0/pi);
    const int bins=static_cast<int>(labEdges.size())-1;
    TH1D* model=keep.keep(new TH1D("lab_cross_section",
        "Elastic cross section, fixed target (e^{+} beam on e^{-} at rest);"
        "#theta_{lab} of the positron [deg];d#sigma/d#Omega_{lab} [barn/sr]",bins,labEdges.data()));
    TH1D* ruth=keep.keep(new TH1D("lab_rutherford","Rutherford",bins,labEdges.data()));
    TH1D* qed=keep.keep(new TH1D("lab_bhabha","Bhabha",bins,labEdges.data()));
    for(TH1D* h: {model,ruth,qed}) { h->SetDirectory(nullptr); h->SetStats(false); }
    double minimum=1e300;
    for(int i=0;i<bins;++i) {
        const double lo=labEdges[i]*pi/180.0, hi=labEdges[i+1]*pi/180.0;
        const double omega=2.0*pi*(std::cos(lo)-std::cos(hi));
        model->SetBinContent(i+1,b.model[i]/(omega*beamLabBarn));
        model->SetBinError(i+1,b.modelError[i]/(omega*beamLabBarn));
        ruth->SetBinContent(i+1,b.rutherford[i]/(omega*beamLabBarn));
        qed->SetBinContent(i+1,b.qed[i]/(omega*beamLabBarn));
        if(b.rutherford[i]>0.0) minimum=std::min(minimum,b.rutherford[i]/(omega*beamLabBarn));
    }
    model->SetLineColor(plot_style::crem());
    model->SetMarkerColor(plot_style::crem());
    model->SetMarkerStyle(20);
    ruth->SetLineColor(plot_style::theory());
    ruth->SetLineStyle(2);
    ruth->SetLineWidth(2);
    qed->SetLineColor(plot_style::qed());
    qed->SetLineStyle(3);
    qed->SetLineWidth(3);
    model->SetMinimum(0.3*minimum);
    model->Draw("E1");
    ruth->Draw("HIST SAME");
    qed->Draw("HIST SAME");
    char buffer[200];
    std::snprintf(buffer,sizeof buffer,"T_{lab} = %.4g eV, #gamma_{cm} = %.6g",
                  f.labKinetic/eCharge,f.gammaCm);
    labNote(keep,0.30,0.70,0.98,0.90,{buffer,
        "dashed: Rutherford (exact QM), dotted: Bhabha (QED)",
        "Bhabha to 10% at 0.6-1 MeV: Ashkin et al., PR 94 (1954)"});
}

// Pad (experiment 4): recoil-electron energy spectrum d sigma / d y.
inline void drawRecoilEnergySpectrum(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const BeamBinCrossSections& b,const BeamLabFrame& f) {
    pad->cd();
    pad->SetGrid();
    pad->SetLogy();
    std::vector<double> yEdges;
    for(double t: b.edges) yEdges.push_back(0.5*(1.0-std::cos(t)));
    const int bins=static_cast<int>(yEdges.size())-1;
    TH1D* model=keep.keep(new TH1D("recoil_energy_spectrum",
        "Recoil-electron energy, fixed target;y = T_{e}/T_{lab};d#sigma/dy [barn]",
        bins,yEdges.data()));
    TH1D* ruth=keep.keep(new TH1D("recoil_rutherford","Rutherford",bins,yEdges.data()));
    TH1D* qed=keep.keep(new TH1D("recoil_bhabha","Bhabha",bins,yEdges.data()));
    for(TH1D* h: {model,ruth,qed}) { h->SetDirectory(nullptr); h->SetStats(false); }
    double minimum=1e300;
    for(int i=0;i<bins;++i) {
        const double width=yEdges[i+1]-yEdges[i];
        model->SetBinContent(i+1,b.model[i]/(width*beamLabBarn));
        model->SetBinError(i+1,b.modelError[i]/(width*beamLabBarn));
        ruth->SetBinContent(i+1,b.rutherford[i]/(width*beamLabBarn));
        qed->SetBinContent(i+1,b.qed[i]/(width*beamLabBarn));
        if(b.rutherford[i]>0.0) minimum=std::min(minimum,b.rutherford[i]/(width*beamLabBarn));
    }
    model->SetLineColor(plot_style::crem());
    model->SetMarkerColor(plot_style::crem());
    model->SetMarkerStyle(20);
    ruth->SetLineColor(plot_style::theory());
    ruth->SetLineStyle(2);
    ruth->SetLineWidth(2);
    qed->SetLineColor(plot_style::qed());
    qed->SetLineStyle(3);
    qed->SetLineWidth(3);
    model->SetMinimum(0.3*minimum);
    model->Draw("E1");
    ruth->Draw("HIST SAME");
    qed->Draw("HIST SAME");
    char buffer[160];
    std::snprintf(buffer,sizeof buffer,"T_{lab} = %.4g eV; y = (1 - cos#theta*)/2",
                  f.labKinetic/eCharge);
    labNote(keep,0.30,0.78,0.98,0.90,{buffer,
        "the quantity Bhabha experiments record (recoil spectrum)"});
}

// Pad (experiment 3): Dirac's annihilation-in-flight cross section against
// the model's cross section for reaching the Compton barrier (its proxy for
// annihilation; the classical model has no in-flight annihilation itself).
inline void drawAnnihilationInFlight(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const BeamLabFrame& f,double modelCrossSection,double modelError) {
    pad->cd();
    pad->SetGrid();
    pad->SetLogx();
    pad->SetLogy();
    std::vector<double> x,y;
    for(int i=0;i<=120;++i) {
        const double tEv=std::pow(10.0,-1.0+i*0.08);  // 0.1 eV .. 1e8.6 eV
        const double gamma=1.0+tEv*eCharge/f.restEnergy;
        x.push_back(tEv);
        y.push_back(diracAnnihilationCrossSection(gamma)/beamLabBarn);
    }
    TGraph* dirac=keep.keep(new TGraph(static_cast<int>(x.size()),x.data(),y.data()));
    dirac->SetTitle("Annihilation in flight (Dirac, QED) vs the model's barrier cross section;"
                    "T_{lab} of the positron [eV];#sigma [barn]");
    dirac->SetLineColor(plot_style::qed());
    dirac->SetLineWidth(2);
    dirac->Draw("AL");
    dirac->SetMinimum(0.3*y.back());
    dirac->SetMaximum(30.0*std::max(y.front(),modelCrossSection/beamLabBarn));
    const double tLabEv=f.labKinetic/eCharge;
    const double xm[1]={tLabEv}, ym[1]={modelCrossSection/beamLabBarn},
                 ex[1]={0.0}, ey[1]={modelError/beamLabBarn};
    TGraphErrors* model=keep.keep(new TGraphErrors(1,xm,ym,ex,ey));
    model->SetMarkerColor(plot_style::crem());
    model->SetMarkerStyle(20);
    model->SetMarkerSize(1.4);
    model->Draw("P SAME");
    char buffer[200];
    std::snprintf(buffer,sizeof buffer,"model at T_{lab} = %.4g eV: %.4g barn; Dirac: %.4g barn",
                  tLabEv,modelCrossSection/beamLabBarn,
                  diracAnnihilationCrossSection(1.0+f.labKinetic/f.restEnergy)/beamLabBarn);
    labNote(keep,0.14,0.15,0.90,0.30,{buffer,
        "theory only: measured at >= 50 MeV, outside the model"});
}
