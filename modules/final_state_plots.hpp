#pragma once

// PLOTS OF EXPERIMENT 6 (audit 391, stage (b) of the distribution redesign;
// README, "Przebudowa rozkładów").  Until now the experiment printed text
// only.  Every quantity is in the pair's rest frame, where the laboratory
// analyses quote them (J-PET reconstructs the decay plane and the spin per
// event).  What is the model's and what is imported:
//   decay times          model (QR rate, w, cascade) against the measured
//                        lifetimes;
//   Dalitz, spectrum,    the Ore-Powell distribution and the 1 - cos^2/3
//   plane normal         plane-normal law are QED imports of the generator,
//                        so these panels check the generator and its
//                        kinematics, not the dynamics;
//   CPT correlation      O_CPT = S.(k1 x k2)/|k1 x k2|, |k1| > |k2| > |k3|:
//                        the generator is CPT-symmetric, so its mean must be
//                        zero within statistics, as J-PET measured
//                        (Moskal et al., Nat. Commun. 12, 5658 (2021):
//                        <O_CPT> = 0.00025 +/- 0.00036, C_CPT = <O_CPT>/P =
//                        0.00067 +/- 0.00095).

#include <TCanvas.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>

inline void drawFinalStateDecayTimes(TVirtualPad* pad,LabPlotKeepAlive& keep,
        const std::vector<double>& times,double measuredTau,double measuredError,
        const char* unit,const char* name,const std::string& title,
        const std::string& reference) {
    LabAnnihilationSample sample;
    sample.times=times;
    drawAnnihilationTimeSpectrum(pad,keep,sample,measuredTau,measuredError,unit,
                                 reference);
    if(TH1D* h=dynamic_cast<TH1D*>(pad->GetPrimitive("annihilation_time_spectrum"))) {
        h->SetName(name);
        h->SetTitle((title+";t ["+unit+"];events / bin").c_str());
    }
}

inline int plotFinalStateAnnihilation(const std::vector<FinalStatePair>& pairs) {
    std::vector<double> paraTimes,orthoTimes,fractions,normalCosine,cpt;
    std::vector<double> dalitzMax,dalitzMid;
    for(const FinalStatePair& pair: pairs) {
        if(pair.para.valid) paraTimes.push_back(pair.para.decaySeconds*1.0e12);
        if(pair.ortho.valid) orthoTimes.push_back(pair.ortho.decaySeconds*1.0e9);
        for(const FinalStateDecay* decay: {&pair.para,&pair.ortho}) {
            if(!decay->valid||decay->photons.count!=3) continue;
            double x[3]={decay->photons.fraction[0],decay->photons.fraction[1],
                         decay->photons.fraction[2]};
            for(double v: x) fractions.push_back(v);
            std::sort(x,x+3);
            dalitzMax.push_back(x[2]);
            dalitzMid.push_back(x[1]);
            if(decay->spinOverHbar.norm()>1e-12) {
                normalCosine.push_back(decay->photons.normalSpinCosine);
                cpt.push_back(decay->photons.orderedSpinCorrelation);
            }
        }
    }
    gROOT->SetBatch(kTRUE);
    root_export::preparePdfExporter();
    gStyle->SetOptStat(0);
    LabPlotKeepAlive keep;
    TCanvas canvas("final_state_statistics","Final-state annihilation",1280,900);
    TPad decayPage("final_state_decay_page","Decay and three-photon kinematics",
                   0.0,0.0,1.0,1.0);
    TPad spinPage("final_state_spin_page","Spin correlations",0.0,0.0,1.0,1.0);
    for(TPad* page: {&decayPage,&spinPage}) {
        page->SetFillColor(kWhite);
        canvas.cd();
        page->Draw();
    }
    decayPage.cd();
    decayPage.Divide(2,2,0.006,0.006);
    spinPage.cd();
    spinPage.Divide(2,1,0.006,0.006);

    const auto& para=statistics_archive::scientificValue("para_lifetime_from_rate");
    const auto& ortho=statistics_archive::scientificValue("ortho_lifetime_from_rate");
    drawFinalStateDecayTimes(decayPage.GetPad(1),keep,paraTimes,para.value,
        para.totalUncertainty,"ps","decay_time_para",
        "p-Ps decay time (pair frame)","Al-Ramadhan & Gidley, PRL 72, 1632 (1994)");
    drawFinalStateDecayTimes(decayPage.GetPad(2),keep,orthoTimes,ortho.value,
        ortho.totalUncertainty,"ns","decay_time_ortho",
        "o-Ps decay time (pair frame)","Vallery, Zitzewitz & Gidley, PRL 90, 203402 (2003)");

    decayPage.cd(3);
    gPad->SetGrid();
    TH2D* dalitz=keep.keep(new TH2D("three_photon_dalitz",
        "Three-photon Dalitz plot (Ore-Powell generator);x_{max} = E_{max}/(W/2);"
        "x_{mid} = E_{mid}/(W/2)",40,2.0/3.0,1.0,40,0.5,1.0));
    dalitz->SetDirectory(nullptr);
    for(size_t i=0;i<dalitzMax.size();++i) dalitz->Fill(dalitzMax[i],dalitzMid[i]);
    dalitz->Draw("COLZ");

    decayPage.cd(4);
    gPad->SetGrid();
    TH1D* spectrum=keep.keep(new TH1D("three_photon_energy",
        "Three-photon energy fractions (pair frame);x = E_{#gamma}/(W/2);photons / bin",
        50,0.0,1.0));
    spectrum->SetDirectory(nullptr);
    for(double v: fractions) spectrum->Fill(v);
    spectrum->SetLineColor(plot_style::crem());
    spectrum->SetLineWidth(2);
    spectrum->Draw("HIST");
    TF1* ore=keep.keep(new TF1("three_photon_ore_powell",
        [](double* x,double* p){ return p[0]*orePowellSpectrumShape(x[0]); },0.0,1.0,1));
    ore->SetParameter(0,fractions.size()*(1.0/50.0)/orePowellSpectrumNorm());
    ore->SetNpx(400);
    ore->SetLineColor(plot_style::qed());
    ore->SetLineStyle(2);
    ore->Draw("SAME");
    {
        const ExponentialFitCheck ks=orePowellKolmogorovSmirnov(fractions);
        char buffer[160];
        std::snprintf(buffer,sizeof buffer,"KS vs Ore-Powell: D = %.3g, p = %.3g (%zu photons)",
                      ks.d,ks.p,fractions.size());
        labNote(keep,0.14,0.78,0.80,0.90,{buffer,
            "measured: Chang, Tang & Li, PLB 157 (1985) -- consistent with QED"});
    }

    spinPage.cd(1);
    gPad->SetGrid();
    TH1D* normal=keep.keep(new TH1D("plane_normal_vs_spin",
        "Decay-plane normal vs o-Ps spin;cos#theta_{n,S};decays / bin",20,-1.0,1.0));
    normal->SetDirectory(nullptr);
    for(double v: normalCosine) normal->Fill(v);
    normal->SetLineColor(plot_style::crem());
    normal->SetLineWidth(2);
    normal->SetMinimum(0.0);
    normal->SetMaximum(1.6*std::max(1.0,normal->GetMaximum()));
    normal->Draw("HIST");
    TF1* law=keep.keep(new TF1("plane_normal_law","[0]*(1-x*x/3)",-1.0,1.0));
    law->SetParameter(0,normalCosine.size()*0.1/(16.0/9.0));
    law->SetLineColor(plot_style::qed());
    law->SetLineStyle(2);
    law->Draw("SAME");
    labNote(keep,0.14,0.80,0.80,0.89,{"1 - cos^{2}#theta/3: QED law of the generator (import)"});

    spinPage.cd(2);
    gPad->SetGrid();
    TH1D* correlation=keep.keep(new TH1D("cpt_correlation",
        "CPT-odd correlation O_{CPT} = #hat{S}#upoint(k_{1}#times k_{2})/|k_{1}#times k_{2}|;"
        "O_{CPT};decays / bin",20,-1.0,1.0));
    correlation->SetDirectory(nullptr);
    double mean=0.0,var=0.0;
    for(double v: cpt) { correlation->Fill(v); mean+=v; }
    mean=cpt.empty()?0.0:mean/cpt.size();
    for(double v: cpt) var+=(v-mean)*(v-mean);
    const double error=cpt.size()>1?std::sqrt(var/(cpt.size()-1)/cpt.size()):0.0;
    correlation->SetLineColor(plot_style::crem());
    correlation->SetLineWidth(2);
    correlation->SetMinimum(0.0);
    correlation->SetMaximum(1.7*std::max(1.0,correlation->GetMaximum()));
    correlation->Draw("HIST");
    {
        char buffer[200];
        std::snprintf(buffer,sizeof buffer,"model: <O_{CPT}> = %.4f #pm %.4f (N = %zu, spin known: P = 1)",
                      mean,error,cpt.size());
        labNote(keep,0.12,0.72,0.98,0.89,{buffer,
            "J-PET: <O_{CPT}> = 0.00025 #pm 0.00036, C_{CPT} = 0.00067 #pm 0.00095",
            "Moskal et al., Nat. Commun. 12, 5658 (2021)"});
    }
    std::cout<<std::setprecision(6)<<"Experiment 6 plots (audit 391): <O_CPT> = "<<mean
             <<" +/- "<<error<<" over "<<cpt.size()<<" o-Ps 3 gamma decays; "
             <<fractions.size()<<" three-photon energy fractions\n";
    canvas.Modified();
    canvas.Update();
    reportExports(root_export::saveStatisticalPlots(6,{
        {decayPage.GetPad(1),1,'b',1,"decay_time_para"},
        {decayPage.GetPad(2),1,'b',2,"decay_time_ortho"},
        {decayPage.GetPad(3),1,'b',3,"three_photon_dalitz"},
        {decayPage.GetPad(4),1,'b',4,"three_photon_energy"},
        {spinPage.GetPad(1),2,'b',1,"plane_normal_vs_spin"},
        {spinPage.GetPad(2),2,'b',2,"cpt_correlation"}}));
    return reportArchiveOperation(statistics_archive::writeScientificReferencesText(),
                                  "scientific-reference catalogue")?0:3;
}
