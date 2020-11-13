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
const bool write_out = false;
char *plot_dir = "Paper_plots/";





void draw_cmp(){

    printf("Year is %i \n", year);
    setTDRStyle();

    init_emu(year);
    init_emu_indv_bkgs(year);
    setup_all_SFs(year);

                                

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

    int n_cost_bins = 8;
    TH1F *data_cost = new TH1F("data_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *diboson_cost = new TH1F("diboson_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *wt_cost = new TH1F("wt_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *dy_cost = new TH1F("dy_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *qcd_cost = new TH1F("qcd_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);


    int n_rap_bins = 20;
    TH1F *data_rap = new TH1F("data_rap", "Data", n_rap_bins, -2.5,2.5);
    TH1F *ttbar_rap = new TH1F("ttbar_rap", "TTbar Background", n_rap_bins, -2.5,2.5);
    TH1F *diboson_rap = new TH1F("diboson_rap", "DiBoson (WW, WZ,ZZ)", n_rap_bins, -2.5,2.5);
    TH1F *wt_rap = new TH1F("wt_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *dy_rap = new TH1F("dy_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *qcd_rap = new TH1F("QCD_rap", "QCD", n_rap_bins, -2.5,2.5);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, 0, 100.);

    Double_t m_low = 150;
    Double_t m_high = 10000;

    bool ss = false;

    bool do_emu_cost_rw = false;
    make_emu_m_cost_pt_rap_hist(t_emu_data, data_m, data_cost, data_pt, data_rap, true,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_rap_hist(t_emu_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_rap, false,  year, m_low, m_high, ss, do_emu_cost_rw);
    make_emu_m_cost_pt_rap_hist(t_emu_diboson, diboson_m, diboson_cost, diboson_pt, diboson_rap, false,  year, m_low, m_high, ss, do_emu_cost_rw);
    make_emu_m_cost_pt_rap_hist(t_emu_wt, wt_m, wt_cost, wt_pt, wt_rap, false,  year, m_low, m_high, ss, do_emu_cost_rw);
    make_emu_m_cost_pt_rap_hist(t_emu_dy, dy_m, dy_cost, dy_pt, dy_rap, false,  year, m_low, m_high, ss);

    bool fakes_reweight = true;
    Fakerate_est_emu(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam, qcd_m,  qcd_cost, qcd_pt, qcd_rap, FLAG_MUONS, year, m_low, m_high, fakes_reweight);


    Double_t data_count = data_m->Integral();
    Double_t fake_count = qcd_m->Integral();
    Double_t mc_count = ttbar_m->Integral() + diboson_m->Integral() + wt_m->Integral() + dy_m->Integral();

    //includes lumi unc
    float fake_err = 0.5;
    double top_unc = pow(0.05*0.05 + 0.025*0.025, 0.5);
    double diboson_unc = pow(0.09*0.09 + 0.025*0.025, 0.5);
    double dy_unc = pow(0.03*0.03 + 0.025*0.025, 0.5);
    Double_t fake_unc = fake_err * fake_count;
    Double_t mc_unc = (ttbar_m->Integral() + wt_m->Integral()) * top_unc + diboson_m->Integral() * diboson_unc + dy_m->Integral()*dy_unc;

    printf("Data count %.0f +/- %.0f \n", data_count, sqrt(data_count));
    printf("MC count %.0f +/- %0.f \n", mc_count, mc_unc);
    printf("Fake count %.0f +/- %.0f \n", fake_count, fake_unc);
    Double_t ratio = data_count / (mc_count + fake_count);
    Double_t unc = sqrt( (data_count/(mc_count + fake_count)/(mc_count + fake_count)) +  
                          pow(data_count/(mc_count + fake_count)/(mc_count + fake_count), 2) * (fake_unc*fake_unc + mc_unc*mc_unc));
    printf("Ratio is %1.3f +/- %1.3f \n", ratio, unc);

    setHistError(qcd_m, fake_err);
    setHistError(qcd_cost, fake_err);
    setHistError(qcd_pt, fake_err);

    setHistError(diboson_m, diboson_unc);
    setHistError(diboson_cost, diboson_unc);
    setHistError(diboson_pt, diboson_unc);

    setHistError(dy_m, dy_unc);
    setHistError(dy_cost, dy_unc);
    setHistError(dy_pt, dy_unc);

    setHistError(ttbar_m, top_unc);
    setHistError(ttbar_cost, top_unc);
    setHistError(ttbar_pt, top_unc);

    setHistError(wt_m, top_unc);
    setHistError(wt_cost, top_unc);
    setHistError(wt_pt, top_unc);


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

    dy_rap->SetFillColor(kRed+1);
    ttbar_rap->SetFillColor(kBlue);
    wt_rap->SetFillColor(kOrange+7); 
    diboson_rap->SetFillColor(kGreen+3);
    qcd_rap->SetFillColor(kRed -7);

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

    symmetrize1d(diboson_cost);
    symmetrize1d(qcd_cost);
    symmetrize1d(wt_cost);
    symmetrize1d(dy_cost);
    symmetrize1d(ttbar_cost);

    THStack *cost_stack = new THStack("cost_stack", "EMu Cos(theta) Distribution: Data vs MC ; cos(#theta)");
    cost_stack->Add(diboson_cost);
    cost_stack->Add(qcd_cost);
    cost_stack->Add(dy_cost);
    cost_stack->Add(wt_cost);
    cost_stack->Add(ttbar_cost);

    THStack *rap_stack = new THStack("rap_stack", "EMu Cos(theta) Distribution: Data vs MC ; cos(#theta)");
    rap_stack->Add(diboson_rap);
    rap_stack->Add(qcd_rap);
    rap_stack->Add(dy_rap);
    rap_stack->Add(wt_rap);
    rap_stack->Add(ttbar_rap);

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
    std::tie(c_pt, p_pt) = make_stack_ratio_plot(data_pt, pt_stack, leg3, "pt", "dilepton pt (GeV)", -1., logy, logx);
    CMS_lumi(p_pt, year, 33);
    sprintf(plt_file, "%sEMu%i_pt_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_pt->Print(plt_file);

    logy = true;
    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg3, "rap", "dilepton rapidity", -1., logy, logx);
    CMS_lumi(p_rap, year, 33);
    sprintf(plt_file, "%sEMu%i_rap_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_rap->Print(plt_file);


}









