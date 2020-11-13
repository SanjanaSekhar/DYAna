
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
#include "../../utils/root_files.h"
#include "../../utils/PlotUtils.C"

const int type = FLAG_ELECTRONS;
const int year = 2018;
const bool write_out = false;
char *plot_dir = "Paper_plots/";


void draw_cmp(){
    setTDRStyle();
    init(year);
    init_indv_bkgs(year);

    int n_pt_bins1 = 7;
    Float_t pt_bins1[] = {0., 10., 20., 30., 50., 70., 100., 300., 700. };

    TH1F *mc_pt = new TH1F("mc_pt", "MC signal", n_pt_bins1, pt_bins1);
    TH1F *mc_tautau_pt = new TH1F("mc_tautau_pt", "MC signal", n_pt_bins1, pt_bins1);
    TH1F *data_pt = new TH1F("data_pt", "MC signal", n_pt_bins1, pt_bins1);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "MC signal", n_pt_bins1, pt_bins1);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC signal", n_pt_bins1, pt_bins1);
    TH1F *wt_pt = new TH1F("wt_pt", "MC signal", n_pt_bins1, pt_bins1);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC signal", n_pt_bins1, pt_bins1);
    TH1F *gg_pt = new TH1F("gg_pt", "MC signal", n_pt_bins1, pt_bins1);

    int n_xf_bins1 = 4;
    float xf_bins1[] = {0.,0.04, 0.07, 0.1, 0.5};
    TH1F *mc_xf = new TH1F("mc_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *mc_tautau_xf = new TH1F("mc_tautau_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *data_xf = new TH1F("data_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *diboson_xf = new TH1F("diboson_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *wt_xf = new TH1F("wt_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *QCD_xf = new TH1F("QCD_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *gg_xf = new TH1F("gg_xf", "MC signal", n_xf_bins1,  xf_bins1);

    /*
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_m_bins, m_bins);
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, m_bins);
    TH1F *mc_tautau_m = new TH1F("mc_tautau_m", "MC no signal (qq, gluglu qbarqbar)", n_m_bins, m_bins);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_m_bins, m_bins);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_m_bins, m_bins);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_m_bins, m_bins);
    TH1F *gg_m = new TH1F("gg_m", "QCD", n_m_bins, m_bins);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_m_bins, m_bins);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_m_bins, m_bins);
    */
    int n_m_bins = 30;
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_m_bins, 150, 2000);
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, 150, 2000);
    TH1F *mc_tautau_m = new TH1F("mc_tautau_m", "MC no signal (qq, gluglu qbarqbar)", n_m_bins, 150, 2000);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_m_bins, 150, 2000);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_m_bins, 150, 2000);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_m_bins, 150, 2000);
    TH1F *gg_m = new TH1F("gg_m", "QCD", n_m_bins, 150, 2000);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_m_bins, 150, 2000);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_m_bins, 150, 2000);

    int num_cost_bins = 10;
    TH1F *data_cost = new TH1F("data_cost", "Data", n_cost_bins, -1.,1.);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1,1);
    TH1F *mc_tautau_cost = new TH1F("mc_tautau_cost", "MC no signal (qq, gluglu qbarqbar)", n_cost_bins, -1.,1.);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", n_cost_bins, -1.,1.);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", n_cost_bins, -1,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", n_cost_bins, -1,1);
    TH1F *gg_cost = new TH1F("gg_cost", "QCD", n_cost_bins, -1,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1,1);
    TH1F *wt_cost = new TH1F("wt_cost", "tw + #bar{t}w", n_cost_bins, -1,1);

    int n_phi_bins = 20;
    TH1F *data_phi = new TH1F("data_phi", "Data", n_phi_bins, -4.,4.);
    TH1F *mc_phi = new TH1F("mc_phi", "MC Signal (qqbar, qglu, qbarglu)", n_phi_bins, -4,4);
    TH1F *mc_tautau_phi = new TH1F("mc_tautau_phi", "MC no signal (qq, gluglu qbarqbar)", n_phi_bins, -4.,4.);
    TH1F *ttbar_phi = new TH1F("ttbar_phi", "TTbar Background", n_phi_bins, -4.,4.);
    TH1F *diboson_phi = new TH1F("diboson_phi", "DiBoson (WW, WZ,ZZ)", n_phi_bins, -4,4);
    TH1F *QCD_phi = new TH1F("QCD_phi", "QCD", n_phi_bins, -4,4);
    TH1F *gg_phi = new TH1F("gg_phi", "QCD", n_phi_bins, -4,4);
    TH1F *WJets_phi = new TH1F("WJets_phi", "WJets", n_phi_bins, -4,4);
    TH1F *wt_phi = new TH1F("wt_phi", "tw + #bar{t}w", n_phi_bins, -4,4);

    int n_rap_bins = 20;
    TH1F *data_rap = new TH1F("data_rap", "Data", n_rap_bins, -2.5,2.5);
    TH1F *mc_rap = new TH1F("mc_rap", "MC Signal (qqbar, qglu, qbarglu)", n_rap_bins, -2.5,2.5);
    TH1F *mc_tautau_rap = new TH1F("mc_tautau_rap", "MC no signal (qq, gluglu qbarqbar)", n_rap_bins, -2.5,2.5);
    TH1F *ttbar_rap = new TH1F("ttbar_rap", "TTbar Background", n_rap_bins, -2.5,2.5);
    TH1F *diboson_rap = new TH1F("diboson_rap", "DiBoson (WW, WZ,ZZ)", n_rap_bins, -2.5,2.5);
    TH1F *QCD_rap = new TH1F("QCD_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *gg_rap = new TH1F("gg_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *WJets_rap = new TH1F("WJets_rap", "WJets", n_rap_bins, -2.5,2.5);
    TH1F *wt_rap = new TH1F("wt_rap", "tw + #bar{t}w", n_rap_bins, -2.5,2.5);

    mc_tautau_cost->SetFillColor(kMagenta);
    mc_tautau_m->SetFillColor(kMagenta);
    mc_tautau_pt->SetFillColor(kMagenta);
    mc_tautau_xf->SetFillColor(kMagenta);
    mc_tautau_phi->SetFillColor(kMagenta);
    mc_tautau_rap->SetFillColor(kMagenta);

    mc_cost->SetFillColor(kRed+1);
    mc_m->SetFillColor(kRed+1);
    mc_pt->SetFillColor(kRed+1);
    mc_xf->SetFillColor(kRed+1);
    mc_phi->SetFillColor(kRed+1);
    mc_rap->SetFillColor(kRed+1);

    ttbar_cost->SetFillColor(kBlue);
    ttbar_m->SetFillColor(kBlue);
    ttbar_pt->SetFillColor(kBlue);
    ttbar_xf->SetFillColor(kBlue);
    ttbar_phi->SetFillColor(kBlue);
    ttbar_rap->SetFillColor(kBlue);

    wt_cost->SetFillColor(kOrange+7);
    wt_m->SetFillColor(kOrange+7);
    wt_pt->SetFillColor(kOrange+7);
    wt_xf->SetFillColor(kOrange+7);
    wt_phi->SetFillColor(kOrange+7);
    wt_rap->SetFillColor(kOrange+7);

    diboson_cost->SetFillColor(kGreen+3);
    diboson_m->SetFillColor(kGreen+3);
    diboson_pt->SetFillColor(kGreen+3);
    diboson_xf->SetFillColor(kGreen+3);
    diboson_phi->SetFillColor(kGreen+3);
    diboson_rap->SetFillColor(kGreen+3);

    QCD_cost->SetFillColor(kRed-7);
    QCD_m->SetFillColor(kRed-7);
    QCD_pt->SetFillColor(kRed-7);
    QCD_xf->SetFillColor(kRed-7);
    QCD_phi->SetFillColor(kRed-7);
    QCD_rap->SetFillColor(kRed-7);

    gg_cost->SetFillColor(kOrange);
    gg_m->SetFillColor(kOrange);
    gg_pt->SetFillColor(kOrange);
    gg_xf->SetFillColor(kOrange);
    gg_phi->SetFillColor(kOrange);
    gg_rap->SetFillColor(kOrange);

    float m_low = 150.;
    float m_high = 10000.;

    make_m_cost_pt_xf_hist(t_elel_data, data_m, data_cost, data_pt, data_xf, data_phi, data_rap, true, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_mc, mc_m, mc_cost, mc_pt, mc_xf, mc_phi, mc_rap, false, type,   year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_tautau, mc_tautau_m, mc_tautau_cost, mc_tautau_pt, mc_tautau_xf, mc_tautau_phi, mc_tautau_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_xf, ttbar_phi, ttbar_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_wt, wt_m, wt_cost, wt_pt, wt_xf, wt_phi, wt_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_gamgam, gg_m, gg_cost, gg_pt, gg_xf, gg_phi, gg_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, diboson_phi, diboson_rap, false, type,   year, m_low, m_high);


    symmetrize1d(gg_cost);

    //gg_cost->Scale(0.);
    //gg_m->Scale(0.);

    bool ss_qcd = false;
    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, QCD_phi, QCD_rap, 
            type, year, m_low, m_high, ss_qcd);


    /*
    QCD_m->Scale(0.6);
    QCD_cost->Scale(0.6);
    QCD_pt->Scale(0.6);
    QCD_xf->Scale(0.6);
    */


    printf("Data integral is %.2f \n", data_cost->Integral());
    printf("DY integral is %.2f \n", mc_cost->Integral());
    printf("ttbar integral is %.2f \n", ttbar_cost->Integral());






    float qcd_err = 0.4;
    setHistError(QCD_m, qcd_err);
    setHistError(QCD_cost, qcd_err);
    setHistError(QCD_pt, qcd_err);
    setHistError(QCD_xf, qcd_err);
    setHistError(QCD_rap, qcd_err);
    setHistError(QCD_phi, qcd_err);

    float gam_err = 0.4;
    setHistError(gg_m, gam_err);
    setHistError(gg_cost, gam_err);
    setHistError(gg_pt, gam_err);
    setHistError(gg_xf, gam_err);
    setHistError(gg_rap, gam_err);
    setHistError(gg_phi, gam_err);



    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(QCD_m);
    m_stack->Add(wt_m);
    m_stack->Add(ttbar_m);
    m_stack->Add(gg_m);
    m_stack->Add(mc_tautau_m);
    m_stack->Add(mc_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; Cos(#theta)_{r}");
    cost_stack->Add(diboson_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(wt_cost);
    cost_stack->Add(ttbar_cost);
    cost_stack->Add(gg_cost);
    cost_stack->Add(mc_tautau_cost);
    cost_stack->Add(mc_cost);

    THStack *pt_stack = new THStack("pt_stack", "ElEl Pt Distribution: Data vs MC; DiElectron Pt (GeV); Events / GeV");
    pt_stack->Add(diboson_pt);
    pt_stack->Add(QCD_pt);
    pt_stack->Add(wt_pt);
    pt_stack->Add(ttbar_pt);
    pt_stack->Add(gg_pt);
    pt_stack->Add(mc_tautau_pt);
    pt_stack->Add(mc_pt);

    binwidth_normalize(pt_stack);
    binwidth_normalize(data_pt);

    THStack *xf_stack = new THStack("xf_stack", "Di-electron x_F Distribution: Data vs MC; x_F");
    xf_stack->Add(diboson_xf);
    xf_stack->Add(QCD_xf);
    xf_stack->Add(wt_xf);
    xf_stack->Add(ttbar_xf);
    xf_stack->Add(gg_xf);
    xf_stack->Add(mc_tautau_xf);
    xf_stack->Add(mc_xf);

    THStack *phi_stack = new THStack("phi_stack", "DiElectron Phi Distribution: Data vs MC; #phi");
    phi_stack->Add(diboson_phi);
    phi_stack->Add(QCD_phi);
    phi_stack->Add(wt_phi);
    phi_stack->Add(ttbar_phi);
    phi_stack->Add(gg_phi);
    phi_stack->Add(mc_tautau_phi);
    phi_stack->Add(mc_phi);

    THStack *rap_stack = new THStack("rap_stack", "DiElectron Rapidity Distribution: Data vs MC; y");
    rap_stack->Add(diboson_rap);
    rap_stack->Add(QCD_rap);
    rap_stack->Add(wt_rap);
    rap_stack->Add(ttbar_rap);
    rap_stack->Add(gg_rap);
    rap_stack->Add(mc_tautau_rap);
    rap_stack->Add(mc_rap);


    gStyle->SetLegendBorderSize(0);
    float x_size = 0.3;
    float y_size = 0.3;
    TLegend *leg1 = new TLegend(x_size, y_size);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(mc_m, "DY Signal", "f");
    leg1->AddEntry(mc_tautau_m, "DY to #tau#tau", "f");
    leg1->AddEntry(gg_m, "#gamma#gamma to #mu#mu", "f");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");
    TLegend *leg4 = (TLegend *) leg1->Clone("leg4");
    TLegend *leg5 = (TLegend *) leg1->Clone("leg5");
    TLegend *leg6 = (TLegend *) leg1->Clone("leg6");

 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    TCanvas *c_m, *c_cost, *c_pt, *c_xf, *c_phi, *c_rap;
    TPad *p_m, *p_cost, *p_pt, *p_xf, *p_phi, *p_rap;
    int iPeriod = 4; 
    writeExtraText = false;
    char plt_file[100];


    bool logy = true;
    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{ee} (GeV)", -1., logy);
    CMS_lumi(p_m, year, 33 );
    sprintf(plt_file, "%sElEl%i_m_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_m->Print(plt_file);

    

    float x_start = 0.45;
    float y_start = 0.2;
    leg2->SetX1(x_start);
    leg2->SetX2(x_start+x_size);
    leg2->SetY1(y_start);
    leg2->SetY2(y_start+y_size);

    logy = false;
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "Cos(#theta)", -1., logy);
    CMS_lumi(p_cost, year, 33);
    sprintf(plt_file, "%sElEl%i_cost_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_cost->Print(plt_file);

    logy = true;
    std::tie(c_pt, p_pt) = make_stack_ratio_plot(data_pt, pt_stack, leg3, "pt", "dielectron pt (GeV)", -1., logy);
    CMS_lumi(p_pt, year, 33);
    sprintf(plt_file, "%sElEl%i_pt_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_pt->Print(plt_file);

    std::tie(c_xf, p_xf) = make_stack_ratio_plot(data_xf, xf_stack, leg4, "xf", "x_F (GeV)", -1., logy);
    CMS_lumi(p_xf, year, 33);
    sprintf(plt_file, "%sElEl%i_xf_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_xf->Print(plt_file);

    std::tie(c_phi, p_phi) = make_stack_ratio_plot(data_phi, phi_stack, leg5, "phi", "dielectron #phi", -1., logy);
    CMS_lumi(p_phi, year, 33);
    sprintf(plt_file, "%sElEl%i_phi_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_phi->Print(plt_file);

    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg6, "rap", "dielectron Y", -1., logy);
    CMS_lumi(p_rap, year, 33);
    sprintf(plt_file, "%sElEl%i_rap_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_rap->Print(plt_file);


}

    
    
