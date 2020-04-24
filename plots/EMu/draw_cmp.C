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

int year = 2018;
const bool write_out = true;
char *plot_dir = "Paper_plots/";





void draw_cmp(){

    printf("Year is %i \n", year);
    setTDRStyle();
    init_emu(year);
    init_emu_indv_bkgs(year);

                                

    TH1F *data_m = new TH1F("data_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *ttbar_m = new TH1F("ttbar_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *diboson_m = new TH1F("diboson_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *wt_m = new TH1F("wt_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *dy_m = new TH1F("dy_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *qcd_m = new TH1F("qcd_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);

    TH1F *data_pt = new TH1F("data_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *wt_pt = new TH1F("wt_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *dy_pt = new TH1F("dy_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *qcd_pt = new TH1F("qcd_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);

    int n_cost_bins = 10;
    TH1F *data_cost = new TH1F("data_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *diboson_cost = new TH1F("diboson_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *wt_cost = new TH1F("wt_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *dy_cost = new TH1F("dy_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *qcd_cost = new TH1F("qcd_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, 0, 100.);

    Double_t m_low = 150;
    Double_t m_high = 10000;

    bool ss = false;

    make_emu_m_cost_pt_xf_hist(t_emu_data, data_m, data_pt, data_cost, h_dummy, true,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_emu_ttbar, ttbar_m, ttbar_pt, ttbar_cost, h_dummy, false,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_emu_diboson, diboson_m, diboson_pt, diboson_cost, h_dummy, false,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_emu_wt, wt_m, wt_pt, wt_cost, h_dummy, false,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_emu_dy, dy_m, dy_pt, dy_cost, h_dummy, false,  year, m_low, m_high, ss);

    Fakerate_est_emu(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam, qcd_m,  qcd_pt, qcd_cost, FLAG_MUONS, year, m_low, m_high);

    float qcd_err = 0.5;
    setHistError(qcd_m, qcd_err);
    setHistError(qcd_cost, qcd_err);
    setHistError(qcd_pt, qcd_err);

    Double_t data_count = data_m->Integral();
    Double_t fake_count = qcd_m->Integral();
    Double_t mc_count = ttbar_m->Integral() + diboson_m->Integral() + wt_m->Integral() + dy_m->Integral();

    Double_t fake_unc = 0.40 * fake_count;
    double bkg_scale_unc = pow(0.05*0.05 + 0.025*0.025, 0.5);
    Double_t mc_unc = bkg_scale_unc* mc_count;

    printf("Data count %.0f \n", data_count);
    printf("MC count %.0f \n", mc_count);
    printf("Fake count %.0f \n", fake_count);
    Double_t ratio = data_count / (mc_count + fake_count);
    Double_t unc = sqrt( (data_count/(mc_count + fake_count)/(mc_count + fake_count)) +  
                          pow(data_count/(mc_count + fake_count)/(mc_count + fake_count), 2) * (fake_unc*fake_unc + mc_unc*mc_unc));
    printf("Ratio is %1.3f +/- %1.3f \n", ratio, unc);



    dy_m->SetFillColor(kRed+1);
    ttbar_m->SetFillColor(kBlue);
    wt_m->SetFillColor(kOrange+7); 
    diboson_m->SetFillColor(kGreen+3);
    qcd_m->SetFillColor(kRed -7);

    dy_pt->SetFillColor(kRed+1);
    ttbar_pt->SetFillColor(kBlue);
    wt_pt->SetFillColor(kOrange+7); 
    diboson_pt->SetFillColor(kGreen+3);
    qcd_pt->SetFillColor(kRed -7);

    dy_cost->SetFillColor(kRed+1);
    ttbar_cost->SetFillColor(kBlue);
    wt_cost->SetFillColor(kOrange+7); 
    diboson_cost->SetFillColor(kGreen+3);
    qcd_cost->SetFillColor(kRed -7);

    THStack *m_stack = new THStack("m_stack", "EMu Mass Distribution: Data vs MC ; m_{e#mu} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(qcd_m);
    m_stack->Add(dy_m);
    m_stack->Add(wt_m);
    m_stack->Add(ttbar_m);

    THStack *pt_stack = new THStack("pt_stack", "EMu Mass Distribution: Data vs MC ; pt_{e#mu} (GeV)");
    pt_stack->Add(diboson_pt);
    pt_stack->Add(qcd_pt);
    pt_stack->Add(dy_pt);
    pt_stack->Add(wt_pt);
    pt_stack->Add(ttbar_pt);

    THStack *cost_stack = new THStack("cost_stack", "EMu Cos(theta) Distribution: Data vs MC ; cos(#theta)");
    cost_stack->Add(diboson_cost);
    cost_stack->Add(qcd_cost);
    cost_stack->Add(dy_cost);
    cost_stack->Add(wt_cost);
    cost_stack->Add(ttbar_cost);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(dy_m, "DY #rightarrow #tau#tau", "f");
    leg1->AddEntry(qcd_m, "QCD and W+Jets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");

    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");

    TCanvas *c_m, *c_cost, *c_pt, *c_xf, *c_phi, *c_rap;
    TPad *p_m, *p_cost, *p_pt, *p_xf, *p_phi, *p_rap;
    int iPeriod = 4; 
    writeExtraText = false;
    char plt_file[100];
    
    bool logy = true;
    bool logx = false;

    sprintf(plt_file, "%sEMu%i_m_cmp.pdf", plot_dir, year % 2000);

    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{e#mu} (GeV)", -1., true);
    CMS_lumi(p_m, year, 33 );
    if(write_out) c_m->Print(plt_file);

    logy = false;
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "Cos(#theta)", -1., logy,logx);
    CMS_lumi(p_cost, year, 33);
    sprintf(plt_file, "%sEMu%i_cost_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_cost->Print(plt_file);

    logy = true;
    std::tie(c_pt, p_pt) = make_stack_ratio_plot(data_pt, pt_stack, leg3, "pt", "dimuon pt (GeV)", -1., logy, logx);
    CMS_lumi(p_pt, year, 33);
    sprintf(plt_file, "%sEMu%i_pt_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_pt->Print(plt_file);


}









