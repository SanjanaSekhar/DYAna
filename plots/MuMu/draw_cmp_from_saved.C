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
#include "../../utils/Colors.h"

const int type = FLAG_MUONS;
const int year = 2016;
const bool write_out = true;
char *plot_dir = "Paper_plots/prefit_kinematics/";
char *fin_name = "MuMu/saved_hists.root";
char *plot_label = "";


void draw_cmp_from_saved(){

    char year_str[80];
    sprintf(year_str, "y%i", year);

    TFile *fin = new TFile(fin_name, "READ");

    fin->cd(year_str);


    TH1F *data_m = (TH1F *) gDirectory->Get("data_m");
    TH1F *diboson_m = (TH1F *) gDirectory->Get("diboson_m");
    TH1F *QCD_m = (TH1F *) gDirectory->Get("QCD_m");
    TH1F *ttbar_m = (TH1F *) gDirectory->Get("ttbar_m");
    TH1F *gg_m = (TH1F *) gDirectory->Get("gg_m");
    TH1F *dy_m = (TH1F *) gDirectory->Get("dy_m");
    TH1F *wt_m = (TH1F *) gDirectory->Get("wt_m");
    TH1F *dy_tautau_m = (TH1F *) gDirectory->Get("dy_tautau_m");

    TH1F *data_cost = (TH1F *) gDirectory->Get("data_cost");
    TH1F *diboson_cost = (TH1F *) gDirectory->Get("diboson_cost");
    TH1F *QCD_cost = (TH1F *) gDirectory->Get("QCD_cost");
    TH1F *ttbar_cost = (TH1F *) gDirectory->Get("ttbar_cost");
    TH1F *gg_cost = (TH1F *) gDirectory->Get("gg_cost");
    TH1F *dy_cost = (TH1F *) gDirectory->Get("dy_cost");
    TH1F *wt_cost = (TH1F *) gDirectory->Get("wt_cost");
    TH1F *dy_tautau_cost = (TH1F *) gDirectory->Get("dy_tautau_cost");


    TH1F *data_rap = (TH1F *) gDirectory->Get("data_rap");
    TH1F *diboson_rap = (TH1F *) gDirectory->Get("diboson_rap");
    TH1F *QCD_rap = (TH1F *) gDirectory->Get("QCD_rap");
    TH1F *ttbar_rap = (TH1F *) gDirectory->Get("ttbar_rap");
    TH1F *gg_rap = (TH1F *) gDirectory->Get("gg_rap");
    TH1F *dy_rap = (TH1F *) gDirectory->Get("dy_rap");
    TH1F *wt_rap = (TH1F *) gDirectory->Get("wt_rap");
    TH1F *dy_tautau_rap = (TH1F *) gDirectory->Get("dy_tautau_rap");

    setTDRStyle();
    gStyle->SetLegendBorderSize(0);
    gStyle->SetErrorX(0);


    dy_tautau_cost->SetFillColor(tautau_c);
    dy_tautau_m->SetFillColor(tautau_c);
    dy_tautau_rap->SetFillColor(tautau_c);

    dy_cost->SetFillColor(DY_c);
    dy_m->SetFillColor(DY_c);
    dy_rap->SetFillColor(DY_c);

    ttbar_cost->SetFillColor(ttbar_c);
    ttbar_m->SetFillColor(ttbar_c);
    ttbar_rap->SetFillColor(ttbar_c);


    wt_cost->SetFillColor(wt_c);
    wt_m->SetFillColor(wt_c);
    wt_rap->SetFillColor(wt_c);

    diboson_cost->SetFillColor(diboson_c);
    diboson_m->SetFillColor(diboson_c);
    diboson_rap->SetFillColor(diboson_c);

    QCD_cost->SetFillColor(qcd_c);
    QCD_m->SetFillColor(qcd_c);
    QCD_rap->SetFillColor(qcd_c);

    gg_cost->SetFillColor(gamgam_c);
    gg_m->SetFillColor(gamgam_c);
    gg_rap->SetFillColor(gamgam_c);



    dy_tautau_cost->SetLineColor(tautau_c);
    dy_tautau_m->SetLineColor(tautau_c);
    dy_tautau_rap->SetLineColor(tautau_c);

    dy_cost->SetLineColor(DY_c);
    dy_m->SetLineColor(DY_c);
    dy_rap->SetLineColor(DY_c);

    ttbar_cost->SetLineColor(ttbar_c);
    ttbar_m->SetLineColor(ttbar_c);
    ttbar_rap->SetLineColor(ttbar_c);


    wt_cost->SetLineColor(wt_c);
    wt_m->SetLineColor(wt_c);
    wt_rap->SetLineColor(wt_c);

    diboson_cost->SetLineColor(diboson_c);
    diboson_m->SetLineColor(diboson_c);
    diboson_rap->SetLineColor(diboson_c);

    QCD_cost->SetLineColor(qcd_c);
    QCD_m->SetLineColor(qcd_c);
    QCD_rap->SetLineColor(qcd_c);

    gg_cost->SetLineColor(gamgam_c);
    gg_m->SetLineColor(gamgam_c);
    gg_rap->SetLineColor(gamgam_c);




    dy_tautau_cost->SetFillStyle(tautau_style);
    dy_tautau_m->SetFillStyle(tautau_style);
    dy_tautau_rap->SetFillStyle(tautau_style);

    dy_cost->SetFillStyle(DY_style);
    dy_m->SetFillStyle(DY_style);
    dy_rap->SetFillStyle(DY_style);

    ttbar_cost->SetFillStyle(ttbar_style);
    ttbar_m->SetFillStyle(ttbar_style);
    ttbar_rap->SetFillStyle(ttbar_style);


    wt_cost->SetFillStyle(wt_style);
    wt_m->SetFillStyle(wt_style);
    wt_rap->SetFillStyle(wt_style);

    diboson_cost->SetFillStyle(diboson_style);
    diboson_m->SetFillStyle(diboson_style);
    diboson_rap->SetFillStyle(diboson_style);

    QCD_cost->SetFillStyle(qcd_style);
    QCD_m->SetFillStyle(qcd_style);
    QCD_rap->SetFillStyle(qcd_style);

    gg_cost->SetFillStyle(gamgam_style);
    gg_m->SetFillStyle(gamgam_style);
    gg_rap->SetFillStyle(gamgam_style);



    TH1F *top_comb_cost = (TH1F *) ttbar_cost->Clone("top_comb_cost");
    TH1F *top_comb_m = (TH1F *) ttbar_m->Clone("top_comb_m");
    TH1F *top_comb_rap = (TH1F *) ttbar_rap->Clone("top_comb_rap");

    top_comb_cost->Add(wt_cost);
    top_comb_m->Add(wt_m);
    top_comb_rap->Add(wt_rap);






    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(gg_m);
    m_stack->Add(diboson_m);
    m_stack->Add(QCD_m);
    //m_stack->Add(wt_m);
    //m_stack->Add(ttbar_m);
    m_stack->Add(top_comb_m);
    //m_stack->Add(dy_tautau_m);
    m_stack->Add(dy_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; Cos(#theta)_{r}");
    cost_stack->Add(gg_cost);
    cost_stack->Add(diboson_cost);
    cost_stack->Add(QCD_cost);
    //cost_stack->Add(wt_cost);
    //cost_stack->Add(ttbar_cost);
    cost_stack->Add(top_comb_cost);
    //cost_stack->Add(dy_tautau_cost);
    cost_stack->Add(dy_cost);


    THStack *rap_stack = new THStack("rap_stack", "DiElectron Rapidity Distribution: Data vs MC; y");
    rap_stack->Add(gg_rap);
    rap_stack->Add(diboson_rap);
    rap_stack->Add(QCD_rap);
    //rap_stack->Add(wt_rap);
    //rap_stack->Add(ttbar_rap);
    rap_stack->Add(top_comb_rap);
    //rap_stack->Add(dy_tautau_rap);
    rap_stack->Add(dy_rap);






    float x_size = 0.4;
    float y_size = 0.3;


    //TLegend *leg1 = new TLegend(x_center - x_size/2, y_center - y_size/2, x_center + x_size/2, y_center + y_size/2);
    TLegend *leg1 = new TLegend(x_size, y_size);
    leg1->SetNColumns(2);
    leg1->SetHeader("Dimuon Signal Region");
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(dy_m, "DY Signal", "f");
    //leg1->AddEntry(dy_tautau_m, "DY #rightarrow #tau#tau", "f");
    //leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    //leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(top_comb_m, "t#bar{t} + tW + #bar{t}W", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ  ", "f");
    leg1->AddEntry(gg_m, "#gamma#gamma #rightarrow #mu#mu", "f");
    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");

    leg1->SetX1NDC(0.7);
    leg1->SetX2NDC(0.7);

 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    TCanvas *c_m, *c_cost, *c_rap;
    TPad *p_m, *p_cost, *p_rap;
    int iPeriod = 4; 
    writeExtraText = false;
    char plt_file[100], y_ax_label[100];



    bool logy = true;

    bool logx = false;
    bool draw_sys_uncs = true;
    float ratio_range = 0.5;



    float x_start_m = 0.55;
    float y_start_m = 0.6;
    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);


    float mbin_base = 10.;
    float hmax = 100000;
    if(year == 2016)
        hmax *= 0.625;
    if(year == 2017)
        hmax *= 0.75;

    sprintf(y_ax_label, "Events/%.0f GeV", mbin_base);
    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{#mu#mu} (GeV)",y_ax_label, plot_label, hmax, logy, logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_m, year, 11 );
    sprintf(plt_file, "%sMuMu%i_m_cmp.png", plot_dir, year % 2000);
    if(write_out) c_m->Print(plt_file);

    
    //float x_start_c = 0.57 - x_size/2;
    //float y_start_c = 0.14;
    float x_start_c = 0.55;
    float y_start_c = 0.6;
    leg2->SetX1(x_start_c);
    leg2->SetX2(x_start_c+x_size);
    leg2->SetY1(y_start_c);
    leg2->SetY2(y_start_c+y_size);

    
    logy = false;
    int n_cost_bins = 10;
    float cost_bin_size = 2./n_cost_bins;
    hmax = 40000;

    if(year == 2016)
        hmax *= 0.625;
    if(year == 2017)
        hmax *= 0.75;

    sprintf(y_ax_label, "Events/%.1f", cost_bin_size);
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "cos#theta_{r}",y_ax_label, plot_label,  hmax, logy,logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_cost, year, 11);
    sprintf(plt_file, "%sMuMu%i_cost_cmp.png", plot_dir, year % 2000);
    if(write_out) c_cost->Print(plt_file);


    leg3->SetX1(x_start_c);
    leg3->SetX2(x_start_c+x_size);
    leg3->SetY1(y_start_c);
    leg3->SetY2(y_start_c+y_size);

    int n_rap_bins = 20;
    float rap_bin_size = 5. / n_rap_bins;

    hmax = 28000;
    if(year == 2016)
        hmax *= 0.625;
    if(year == 2017)
        hmax *= 0.75;

    sprintf(y_ax_label, "Events/%.2f", rap_bin_size);
    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg3, "rap", "dimuon rapidity",y_ax_label, plot_label, hmax, logy, logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_rap, year, 11);
    sprintf(plt_file, "%sMuMu%i_rap_cmp.png", plot_dir, year % 2000);
    if(write_out) c_rap->Print(plt_file);
}
