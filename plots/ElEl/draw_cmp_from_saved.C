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

const int type = FLAG_ELECTRONS;
const bool write_out = true;
bool prelim = false;
char *plot_dir = "AN_plots/";
//char *plot_dir = "PAS_plots/prefit_kinematics/";
char *fin_name = "ElEl/LQ_saved_hists.root";
char *plot_label = "";

//gStyle->SetOptStat(0);
//gROOT->SetBatch(1);
void draw_cmp_from_saved(){

    TH1F *data_m ,*diboson_m, *QCD_m, *top_m, *dy_m, *wt_m, *gg_m, *LQu_m, *LQu_vec_m;
    TH1F *data_cost, *diboson_cost, *QCD_cost, *top_cost, *dy_cost, *wt_cost, *gg_cost, *LQu_cost, *LQu_vec_cost;
    TH1F *data_rap, *diboson_rap, *QCD_rap, *top_rap, *dy_rap, *wt_rap, *gg_rap, *LQu_rap, *LQu_vec_rap;

    TFile *fin = new TFile(fin_name, "READ");
    gStyle->SetOptStat(0);
    gROOT->SetBatch(1);
    for(int year = 2016; year <= 2017; year++){

        char year_str[80];
        sprintf(year_str, "y%i", year);


        fin->cd(year_str);


        TH1F *data_m_ = (TH1F *) gDirectory->Get("data_m");
        TH1F *diboson_m_ = (TH1F *) gDirectory->Get("diboson_m");
        TH1F *QCD_m_ = (TH1F *) gDirectory->Get("QCD_m");
        TH1F *ttbar_m_ = (TH1F *) gDirectory->Get("ttbar_m");
        TH1F *dy_m_ = (TH1F *) gDirectory->Get("dy_m");
        TH1F *wt_m_ = (TH1F *) gDirectory->Get("wt_m");
        TH1F *gg_m_ = (TH1F *) gDirectory->Get("gg_m");
        TH1F *LQu_m_ = (TH1F *) gDirectory->Get("LQu_m");
        TH1F *LQu_vec_m_ = (TH1F *) gDirectory->Get("LQu_vec_m");
        

        TH1F *data_cost_ = (TH1F *) gDirectory->Get("data_cost");
        TH1F *diboson_cost_ = (TH1F *) gDirectory->Get("diboson_cost");
        TH1F *QCD_cost_ = (TH1F *) gDirectory->Get("QCD_cost");
        TH1F *ttbar_cost_ = (TH1F *) gDirectory->Get("ttbar_cost");
        TH1F *dy_cost_ = (TH1F *) gDirectory->Get("dy_cost");
        TH1F *wt_cost_ = (TH1F *) gDirectory->Get("wt_cost");
        TH1F *gg_cost_ = (TH1F *) gDirectory->Get("gg_cost");
        TH1F *LQu_cost_ = (TH1F *) gDirectory->Get("LQu_cost");
        TH1F *LQu_vec_cost_ = (TH1F *) gDirectory->Get("LQu_vec_cost");
        

        TH1F *data_rap_ = (TH1F *) gDirectory->Get("data_rap");
        TH1F *diboson_rap_ = (TH1F *) gDirectory->Get("diboson_rap");
        TH1F *QCD_rap_ = (TH1F *) gDirectory->Get("QCD_rap");
        TH1F *ttbar_rap_ = (TH1F *) gDirectory->Get("ttbar_rap");
        TH1F *dy_rap_ = (TH1F *) gDirectory->Get("dy_rap");
        TH1F *wt_rap_ = (TH1F *) gDirectory->Get("wt_rap");
        TH1F *gg_rap_ = (TH1F *) gDirectory->Get("gg_rap");
        TH1F *LQu_rap_ = (TH1F *) gDirectory->Get("LQu_rap");
        TH1F *LQu_vec_rap_ = (TH1F *) gDirectory->Get("LQu_vec_rap");



        TH1F *top_cost_ = (TH1F *) ttbar_cost_->Clone("top_cost");
        TH1F *top_m_ = (TH1F *) ttbar_m_->Clone("top_m");
        TH1F *top_rap_ = (TH1F *) ttbar_rap_->Clone("top_rap");

        top_cost_->Add(wt_cost_);
        top_m_->Add(wt_m_);
        top_rap_->Add(wt_rap_);

        if(year == 2016){
            printf("%i \n", year);
            data_m = data_m_; diboson_m = diboson_m_; QCD_m = QCD_m_; top_m = top_m_; dy_m = dy_m_; wt_m = wt_m_; gg_m = gg_m_; LQu_m = LQu_m_; LQu_vec_m = LQu_vec_m_;
            data_cost = data_cost_; diboson_cost = diboson_cost_; QCD_cost = QCD_cost_; top_cost = top_cost_; dy_cost = dy_cost_; wt_cost = wt_cost_; gg_cost = gg_cost_; LQu_cost = LQu_cost_; LQu_vec_cost = LQu_vec_cost_;
            data_rap = data_rap_; diboson_rap = diboson_rap_; QCD_rap = QCD_rap_; top_rap = top_rap_; dy_rap = dy_rap_; wt_rap = wt_rap_; gg_rap = gg_rap_; LQu_rap = LQu_rap_; LQu_vec_rap = LQu_vec_rap_;

        }
        else{
            printf("%i \n", year);
            data_m->Add(data_m_);
            diboson_m->Add(diboson_m_);
            QCD_m->Add(QCD_m_);
            top_m->Add(top_m_);
            dy_m->Add(dy_m_);
            wt_m->Add(wt_m_);
            gg_m->Add(gg_m_);
            LQu_m->Add(LQu_m_);
            LQu_vec_m->Add(LQu_vec_m_);


            data_cost->Add(data_cost_);
            diboson_cost->Add(diboson_cost_);
            QCD_cost->Add(QCD_cost_);
            top_cost->Add(top_cost_);
            dy_cost->Add(dy_cost_);
            wt_cost->Add(wt_cost_);
            gg_cost->Add(gg_cost_);
            LQu_cost->Add(LQu_cost_);
            LQu_vec_cost->Add(LQu_vec_cost_);           

            data_rap->Add(data_rap_);
            diboson_rap->Add(diboson_rap_);
            QCD_rap->Add(QCD_rap_);
            top_rap->Add(top_rap_);
            dy_rap->Add(dy_rap_);
            wt_rap->Add(wt_rap_);
            gg_rap->Add(gg_rap_);
            LQu_rap->Add(LQu_rap_);
            LQu_vec_rap->Add(LQu_vec_rap_);
        }
    }


    setTDRStyle();
   gStyle->SetLegendBorderSize(0);


    setHistError(QCD_cost, qcd_sys_unc);
    setHistError(QCD_rap, qcd_sys_unc);

    setHistError(dy_cost, dy_sys_unc);
    setHistError(dy_rap, dy_sys_unc);

    setHistError(diboson_cost, diboson_sys_unc);
    setHistError(diboson_rap, diboson_sys_unc);

    setHistError(top_cost, top_sys_unc);
    setHistError(top_rap, top_sys_unc);

    setHistError(wt_cost, top_sys_unc);
    setHistError(wt_rap, top_sys_unc);


    setHistError(gg_cost, gam_sys_unc);
    setHistError(gg_rap, gam_sys_unc);

    setHistError(LQu_cost, LQu_sys_unc);
    setHistError(LQu_rap, LQu_sys_unc);

    setHistError(LQu_vec_cost, LQu_vec_sys_unc);
    setHistError(LQu_vec_rap, LQu_vec_sys_unc);




    setHistError(QCD_m, qcd_sys_unc);

    setHistMassDepError(dy_m);
    setHistMassDepError(diboson_m);
    setHistMassDepError(top_m);
    setHistMassDepError(wt_m);
    setHistMassDepError(wt_m);
    setHistMassDepError(gg_m);
    setHistMassDepError(LQu_m);
    setHistMassDepError(LQu_vec_m);
   

    float mbin_base = 10.;
    binwidth_normalize(data_m, mbin_base);
    binwidth_normalize(diboson_m, mbin_base);
    binwidth_normalize(QCD_m, mbin_base);
    binwidth_normalize(wt_m, mbin_base);
    binwidth_normalize(top_m, mbin_base);
    binwidth_normalize(gg_m, mbin_base);
    binwidth_normalize(dy_m, mbin_base);
    binwidth_normalize(LQu_m, mbin_base);
    binwidth_normalize(LQu_vec_m, mbin_base);



    dy_cost->SetFillColor(DY_c);
    dy_m->SetFillColor(DY_c);
    dy_rap->SetFillColor(DY_c);

    top_cost->SetFillColor(ttbar_c);
    top_m->SetFillColor(ttbar_c);
    top_rap->SetFillColor(ttbar_c);


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




    dy_cost->SetLineColor(DY_c);
    dy_m->SetLineColor(DY_c);
    dy_rap->SetLineColor(DY_c);

    top_cost->SetLineColor(ttbar_c);
    top_m->SetLineColor(ttbar_c);
    top_rap->SetLineColor(ttbar_c);


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


    LQu_m->SetLineColor(kGreen);
    LQu_cost->SetLineColor(kGreen);
    LQu_rap->SetLineColor(kGreen);

    LQu_vec_m->SetLineColor(kRed);
    LQu_vec_cost->SetLineColor(kRed);
    LQu_vec_rap->SetLineColor(kRed);

    LQu_m->SetLineWidth(4);
    LQu_cost->SetLineWidth(4);
    LQu_rap->SetLineWidth(4);

    LQu_vec_m->SetLineWidth(4);
    LQu_vec_cost->SetLineWidth(4);
    LQu_vec_rap->SetLineWidth(4);

    dy_cost->SetFillStyle(DY_style);
    dy_m->SetFillStyle(DY_style);
    dy_rap->SetFillStyle(DY_style);

    top_cost->SetFillStyle(ttbar_style);
    top_m->SetFillStyle(ttbar_style);
    top_rap->SetFillStyle(ttbar_style);


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









    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(gg_m);
    m_stack->Add(diboson_m);
    m_stack->Add(QCD_m);
    m_stack->Add(top_m);
    m_stack->Add(dy_m);
     //m_stack->Add(LQu_m);
     //m_stack->Add(LQu_vec_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; Cos(#theta)_{r}");
    cost_stack->Add(gg_cost);
    cost_stack->Add(diboson_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(top_cost);
    cost_stack->Add(dy_cost);
    //cost_stack->Add(LQu_cost);
    //cost_stack->Add(LQu_vec_cost);


    THStack *rap_stack = new THStack("rap_stack", "DiElectron Rapidity Distribution: Data vs MC; y");
    rap_stack->Add(gg_rap);
    rap_stack->Add(diboson_rap);
    rap_stack->Add(QCD_rap);
    rap_stack->Add(top_rap);
    rap_stack->Add(dy_rap);
    //rap_stack->Add(LQu_rap);
    //rap_stack->Add(LQu_vec_rap);








    float x_size = 0.7;
    float y_size = 0.3;


    //TLegend *leg1 = new TLegend(x_center - x_size/2, y_center - y_size/2, x_center + x_size/2, y_center + y_size/2);
    TLegend *leg1 = new TLegend(x_size, y_size);
    leg1->SetNColumns(2);
    leg1->SetHeader("Dielectron signal region");

    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");

     data_m->SetLineWidth(2);
     leg1->AddEntry(data_m, "Data", "lpe");
     data_m->SetLineWidth(2);
     leg2->AddEntry(data_m, "Data", "pe");
     leg3->AddEntry(data_m, "Data", "pe");

    leg1->AddEntry(dy_m, "DY ", "f");
    leg1->AddEntry(top_m, "t#bar{t} + single t", "f");
    leg1->AddEntry(QCD_m, "QCD and W+jets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ  ", "f");
    //leg1->AddEntry(gg_m, "#gamma#gamma #rightarrow #bf{ee}", "f");
    leg1->AddEntry(gg_m, "#gamma#gamma #rightarrow ee", "f");
    leg1->AddEntry(LQu_m, "2.5 TeV S_{eu} (y_{eu}=2.0)");
    leg1->AddEntry(LQu_vec_m, "2.5 TeV V_{eu} (g_{eu}=1.0)");
    leg1->SetTextSize(0.05);


    leg2->AddEntry(dy_m, "DY ", "f");
    leg2->AddEntry(top_m, "t#bar{t} + single t", "f");
    leg2->AddEntry(QCD_m, "QCD and W+jets", "f");
    leg2->AddEntry(diboson_m, "WW + WZ + ZZ  ", "f");
    //leg2->AddEntry(gg_m, "#gamma#gamma #rightarrow #bf{ee}", "f");
    leg2->AddEntry(gg_m, "#gamma#gamma #rightarrow ee", "f");
    leg2->AddEntry(LQu_m, "2.5 TeV S_{eu} (y_{eu}=2.0)");
    leg2->AddEntry(LQu_vec_m, "2.5 TeV V_{eu} (g_{eu}=1.0)");
    leg2->SetTextSize(0.05);


    leg3->AddEntry(dy_m, "DY ", "f");
    leg3->AddEntry(top_m, "t#bar{t} + single t", "f");
    leg3->AddEntry(QCD_m, "QCD and W+jets", "f");
    leg3->AddEntry(diboson_m, "WW + WZ + ZZ  ", "f");
    //leg3->AddEntry(gg_m, "#gamma#gamma #rightarrow #bf{ee}", "f");
    leg3->AddEntry(gg_m, "#gamma#gamma #rightarrow ee", "f");
    leg3->AddEntry(LQu_m, "2.5 TeV S_{eu} (y_{eu}=2.0)");
    leg3->AddEntry(LQu_vec_m, "2.5 TeV V_{eu} (g_{eu}=1.0)");
    leg3->SetTextSize(0.05);

    leg1->SetX1NDC(0.7);
    leg1->SetX2NDC(0.7);

 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    TCanvas *c_m, *c_cost, *c_rap;
    TPad *p_m, *p_cost, *p_rap;
    int iPeriod = 4; 
    if(prelim) writeExtraText = true;
    else writeExtraText = false;
    char plt_file[100], y_ax_label[100];



    bool logy = true;

    bool logx = false;
    bool draw_sys_uncs = true;
    float ratio_range = 0.5;



    float x_start_m = 0.2;
    float y_start_m = 0.55;
    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);


    int year = -1;
    float hmax = 700;
    //if(year == 2016)
    //    hmax *= 0.625;
    //if(year == 2017)
    //    hmax *= 0.75;

    float hmin = 0.001;
    sprintf(y_ax_label, "Events / %.0f GeV", mbin_base);
    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "m_{ee} (GeV)",y_ax_label, plot_label, hmax, logy, logx, draw_sys_uncs, ratio_range, false, hmin);
    //CMS_lumi(p_m, year, 11 );
    CMS_lumi(c_m, year, 11 );
    p_m->cd();
    LQu_m->Draw("hist  same");
    LQu_m->Print("range");
    LQu_vec_m->Draw("hist  same");
    sprintf(plt_file, "%sElElComb_m_cmp.png", plot_dir);
    if(write_out) c_m->Print(plt_file);
    sprintf(plt_file, "%sElElComb_m_cmp.pdf", plot_dir);
    if(write_out) c_m->Print(plt_file);

    
    //float x_start_c = 0.57 - x_size/2;
    //float y_start_c = 0.14;
    float x_start_c = 0.2;
    float y_start_c = 0.55;
    leg2->SetX1(x_start_c);
    leg2->SetX2(x_start_c+x_size);
    leg2->SetY1(y_start_c);
    leg2->SetY2(y_start_c+y_size);

    
    logy = false;
    int n_cost_bins = 10;
    float cost_bin_size = 2./n_cost_bins;
    hmax = 1300;

    //if(year == 2016)
    //    hmax *= 0.625;
    //if(year == 2017)
    //    hmax *= 0.75;




    sprintf(y_ax_label, "Events / %.1f", cost_bin_size);
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "cos #theta_{R}",y_ax_label, plot_label,  hmax, logy,logx, draw_sys_uncs, ratio_range);
    //CMS_lumi(p_cost, year, 11);
    CMS_lumi(c_cost, year, 11 );
    p_cost->cd();
    LQu_cost->Draw("hist  same");
    LQu_vec_cost->Draw("hist  same");
    sprintf(plt_file, "%sElElComb_cost_cmp.png", plot_dir);
    if(write_out) c_cost->Print(plt_file);
    sprintf(plt_file, "%sElElComb_cost_cmp.pdf", plot_dir);
    if(write_out) c_cost->Print(plt_file);


    leg3->SetX1(x_start_c);
    leg3->SetX2(x_start_c+x_size);
    leg3->SetY1(y_start_c);
    leg3->SetY2(y_start_c+y_size);

    int n_rap_bins = 20;
    float rap_bin_size = 5. / n_rap_bins;

    hmax = 1100;

    //if(year == 2016)
    //    hmax *= 0.625;
    //if(year == 2017)
    //    hmax *= 0.75;



    sprintf(y_ax_label, "Events / %.2f", rap_bin_size);
    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg3, "rap", "Dielectron rapidity",y_ax_label, plot_label, hmax, logy, logx, draw_sys_uncs, ratio_range);
    //CMS_lumi(p_rap, year, 11);
    CMS_lumi(c_rap, year, 11 );
    p_rap->cd();
    LQu_rap->Draw("hist  same");
    LQu_vec_rap->Draw("hist  same");
    sprintf(plt_file, "%sElElComb_rap_cmp.png", plot_dir);
    if(write_out) c_rap->Print(plt_file);
    sprintf(plt_file, "%sElElComb_rap_cmp.pdf", plot_dir);
    if(write_out) c_rap->Print(plt_file);
}
