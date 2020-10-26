#define STAND_ALONE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/HistMaker.C"
#include "../../utils/PlotUtils.C"
#include "../../utils/root_files.h"

const int type = FLAG_MUONS;
const int year = 2018;
const bool write_out = false;
char *plot_dir = "Paper_plots/";


void make_cosmics_plots(TTree *t1, TH1F *h_dphi, TH1F *h_eta, 
        bool is_data=false, int flag1 = FLAG_MUONS,
        int year = 2016, Double_t m_low = 150., Double_t m_high = 9999999., bool ss = false){
    //read event data
        TempMaker tm(t1, is_data, year);
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.do_RC = true;


        tm.setup();
        int nEvents=0;

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);
            bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;
            //bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.has_no_bjets && tm.not_cosmic;
            //bool good_eta = (fabs(tm.el1_eta) > 1.5 && fabs(tm.el2_eta) < 1.4) || (fabs(tm.el2_eta) > 1.5 && fabs(tm.el1_eta) < 1.4);
            //pass = pass && good_eta;


            if(pass){
                nEvents++;
                tm.doCorrections();
                tm.getEvtWeight();
                h_dphi->Fill(std::fabs(tm.lep_p->Phi() - tm.lep_m->Phi()), tm.evt_weight);
                h_eta->Fill(std::fabs(tm.lep_p->Phi() + tm.lep_m->Phi()), tm.evt_weight);
            }


        
        }
        printf("Selected %i events \n", nEvents);


    tm.finish();
}


void draw_cosmics_plot(){

    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    setup_all_SFs(year);


    int n_dphi_bins = 40;
    float dphi_max = TMath::Pi() + 0.01;
    TH1F *data_dphi = new TH1F("data_dphi", "Data Dimuon Mass Distribution", n_dphi_bins, 0., dphi_max);
    TH1F *mc_dphi = new TH1F("mc_dphi", "MC Signal (qqbar, qglu, qbarglu)", n_dphi_bins, 0., dphi_max);
    TH1F *mc_tautau_dphi = new TH1F("mc_tautau_dphi", "MC no signal (qq, gluglu qbarqbar)", n_dphi_bins, 0., dphi_max);
    TH1F *ttbar_dphi = new TH1F("ttbar_dphi", "TTBar Background", n_dphi_bins, 0., dphi_max);
    TH1F *diboson_dphi = new TH1F("diboson_dphi", "DiBoson (WW, WZ, ZZ)", n_dphi_bins, 0., dphi_max);
    TH1F *gg_dphi = new TH1F("gg_dphi", "QCD", n_dphi_bins, 0., dphi_max);
    TH1F *WJets_dphi = new TH1F("WJets_dphi", "WJets", n_dphi_bins, 0., dphi_max);
    TH1F *wt_dphi = new TH1F("wt_dphi", "tw + #bar{t}w", n_dphi_bins, 0., dphi_max);


    int n_eta_bins = 40;
    float eta_max = 2.;
    TH1F *data_eta = new TH1F("data_eta", "Data Dimuon Mass Distribution", n_eta_bins, 0., eta_max);
    TH1F *mc_eta = new TH1F("mc_eta", "MC Signal (qqbar, qglu, qbarglu)", n_eta_bins, 0., eta_max);
    TH1F *mc_tautau_eta = new TH1F("mc_tautau_eta", "MC no signal (qq, gluglu qbarqbar)", n_eta_bins, 0., eta_max);
    TH1F *ttbar_eta = new TH1F("ttbar_eta", "TTBar Background", n_eta_bins, 0., eta_max);
    TH1F *diboson_eta = new TH1F("diboson_eta", "DiBoson (WW, WZ, ZZ)", n_eta_bins, 0., eta_max);
    TH1F *gg_eta = new TH1F("gg_eta", "QCD", n_eta_bins, 0., eta_max);
    TH1F *WJets_eta = new TH1F("WJets_eta", "WJets", n_eta_bins, 0., eta_max);
    TH1F *wt_eta = new TH1F("wt_eta", "tw + #bar{t}w", n_eta_bins, 0., eta_max);



    mc_tautau_dphi->SetFillColor(kMagenta);
    mc_tautau_eta->SetFillColor(kMagenta);

    mc_dphi->SetFillColor(kRed+1);
    mc_eta->SetFillColor(kRed+1);

    ttbar_dphi->SetFillColor(kBlue);
    ttbar_eta->SetFillColor(kBlue);


    wt_dphi->SetFillColor(kOrange+7);
    wt_eta->SetFillColor(kOrange+7);

    diboson_dphi->SetFillColor(kGreen+3);
    diboson_eta->SetFillColor(kGreen+3);

    gg_dphi->SetFillColor(kOrange);
    gg_eta->SetFillColor(kOrange);

    float m_low = 150.;
    float m_high = 1000.;


    make_cosmics_plots(t_mumu_data, data_dphi, data_eta, true, type, year, m_low, m_high);
    make_cosmics_plots(t_mumu_mc, mc_dphi, mc_eta, false, type, year, m_low, m_high);
    make_cosmics_plots(t_mumu_ttbar, ttbar_dphi, ttbar_eta, false, type, year, m_low, m_high);
    make_cosmics_plots(t_mumu_diboson, diboson_dphi, diboson_eta, false, type, year, m_low, m_high);
    make_cosmics_plots(t_mumu_wt, wt_dphi, wt_eta, false, type, year, m_low, m_high);
    make_cosmics_plots(t_mumu_gamgam, gg_dphi, gg_eta, false, type, year, m_low, m_high);


    THStack *dphi_stack = new THStack("dphi_stack", "MuMu Phi Distribution: Data vs MC ; |#Delta #phi|");
    dphi_stack->Add(diboson_dphi);
    dphi_stack->Add(wt_dphi);
    dphi_stack->Add(ttbar_dphi);
    dphi_stack->Add(gg_dphi);
    dphi_stack->Add(mc_tautau_dphi);
    dphi_stack->Add(mc_dphi);


    THStack *eta_stack = new THStack("eta_stack", "#eta Distribution: Data vs MC; |#eta_{1} + #eta_{2}|");
    eta_stack->Add(diboson_eta);
    eta_stack->Add(wt_eta);
    eta_stack->Add(ttbar_eta);
    eta_stack->Add(gg_eta);
    eta_stack->Add(mc_tautau_eta);
    eta_stack->Add(mc_eta);

    gStyle->SetLegendBorderSize(0);
    float x_size = 0.3;
    float y_size = 0.3;
    TLegend *leg1 = new TLegend(x_size, y_size);
    leg1->AddEntry(data_dphi, "data", "p");
    leg1->AddEntry(mc_dphi, "DY Signal", "f");
    leg1->AddEntry(mc_tautau_dphi, "DY to #tau#tau", "f");
    leg1->AddEntry(gg_dphi, "#gamma#gamma to #mu#mu", "f");
    leg1->AddEntry(ttbar_dphi, "t#bar{t}", "f");
    leg1->AddEntry(wt_dphi, "tW + #bar{t}W", "f");
    leg1->AddEntry(diboson_dphi, "WW + WZ + ZZ", "f");
    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");


    TCanvas *c_dphi, *c_eta;
    TPad *p_dphi, *p_eta;
    int iPeriod = 4; 
    writeExtraText = false;
    char plt_file[100];

    std::tie(c_dphi, p_dphi) = make_stack_ratio_plot(data_dphi, dphi_stack, leg1, "dphi", "|#Delta #phi|", -1., false, false);
    CMS_lumi(p_dphi, year, 33);
    sprintf(plt_file, "%sMuMu%i_cosmics_dphi.pdf", plot_dir, year % 2000);
    if(write_out) c_dphi->Print(plt_file);

    std::tie(c_eta, p_eta) = make_stack_ratio_plot(data_eta, eta_stack, leg1, "eta", "|#eta_{1} + #eta_{2}|", -1., false, false);
    CMS_lumi(p_eta, year, 33);
    sprintf(plt_file, "%sMuMu%i_cosmics_eta.pdf", plot_dir, year % 2000);
    if(write_out) c_eta->Print(plt_file);

}
