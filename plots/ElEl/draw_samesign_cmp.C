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
#include "../../utils/Colors.h"

const int type = FLAG_ELECTRONS;
int year = 2016;
bool write_out = true;
//char *plot_dir = "Misc_plots/samesign_cmp_scaled/";
char *plot_dir = "Paper_plots/prefit_kinematics";


void draw_samesign_cmp(){
    init(year);




    setTDRStyle();

    TH1F *data_m = new TH1F("data_m", "Data Dielectron Mass Distribution", 30, 150, 2000);


    int pt_bins = 20.;
    TH1F *data_pt = new TH1F("data_pt", "MC signal", pt_bins, 0, 800);
    TH1F *back_pt = new TH1F("back_pt", "MC signal", pt_bins, 0, 800);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC signal", pt_bins, 0, 800);
    TH1F *DY_pt = new TH1F("DY_pt", "MC signal", pt_bins, 0, 800);

    int xf_nbins = 16;
    TH1F *data_xf = new TH1F("data_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *QCD_xf = new TH1F("QCD_xf", "MC signal", xf_nbins, 0, 0.8);


    int n_rap_bins = 20;
    TH1F *data_rap = new TH1F("data_rap", "Data", n_rap_bins, -2.5,2.5);
    TH1F *DY_rap = new TH1F("mc_rap", "MC Signal (qqbar, qglu, qbarglu)", n_rap_bins, -2.5,2.5);
    TH1F *back_rap = new TH1F("diboson_rap", "DiBoson (WW, WZ,ZZ)", n_rap_bins, -2.5,2.5);
    TH1F *QCD_rap = new TH1F("QCD_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *wt_rap = new TH1F("wt_rap", "tw + #bar{t}w", n_rap_bins, -2.5,2.5);


    TH1F *back_m = new TH1F("back_m", "back (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *back_xf = new TH1F("back_xf", "MC signal", xf_nbins, 0, 0.8);

    TH1F *DY_m = new TH1F("DY_m", "DY (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *DY_xf = new TH1F("DY_xf", "MC signal", xf_nbins, 0, 0.8);

    TH1F *QCD_m = new TH1F("QCD_m", "QCD", 30, 150, 2000);

    TH1F *WJets_m = new TH1F("WJets_m", "WJets", 30, 150, 2000);

    int n_cost_bins = 8;
    TH1F *data_cost = new TH1F("data_cost", "Data", n_cost_bins, -1.,1.);
    TH1F *back_cost = new TH1F("back_cost", "back (WW, WZ,ZZ)", n_cost_bins, -1.,1);
    TH1F *DY_cost = new TH1F("DY_cost", "DY (WW, WZ,ZZ)", n_cost_bins, -1.,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", n_cost_bins, -1.,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1.,1);


    TH1F *DY_phi = new TH1F("DY_phi", "", 20, -4.0, 4.0);
    TH1F *data_phi = new TH1F("data_phi", "", 20, -4.0, 4.0);
    TH1F *QCD_phi = new TH1F("QCD_phi", "", 20, -4.0, 4.0);
    TH1F *back_phi = new TH1F("back_phi", "", 20, -4.0, 4.0);

    TH1F * dummy = new TH1F("dummy", "", 100, 0., 100.);

    back_m->SetFillColor(diboson_c);
    back_cost->SetFillColor(diboson_c);
    back_xf->SetFillColor(diboson_c);
    back_pt->SetFillColor(diboson_c);
    back_phi->SetFillColor(diboson_c);
    back_rap->SetFillColor(diboson_c);

    QCD_xf->SetFillColor(qcd_c);
    QCD_m->SetFillColor(qcd_c);
    QCD_cost->SetFillColor(qcd_c);
    QCD_pt->SetFillColor(qcd_c);
    QCD_phi->SetFillColor(qcd_c);
    QCD_rap->SetFillColor(qcd_c);



    DY_xf->SetFillColor(DY_c);
    DY_pt->SetFillColor(DY_c);
    DY_m->SetFillColor(DY_c);
    DY_cost->SetFillColor(DY_c);
    DY_phi->SetFillColor(DY_c);
    DY_rap->SetFillColor(DY_c);

    float m_low = 150.;
    float m_high = 10000.;
    bool ss = true;

    make_m_cost_pt_xf_hist(t_elel_ss_data, data_m, data_cost, data_pt, data_xf, data_phi, data_rap,  true, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_ttbar, back_m, back_cost, back_pt, back_xf, back_phi, back_rap, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_wt, back_m, back_cost, back_pt, back_xf, back_phi, back_rap, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_diboson, back_m, back_cost, back_pt, back_xf, back_phi, back_rap, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_dy, DY_m, DY_cost, DY_pt, DY_xf, DY_phi, DY_rap, false, type,   year, m_low, m_high, ss);



    bool cost_reweight = false;
    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, QCD_phi, QCD_rap, type,  year, m_low, m_high, ss, cost_reweight);
    
    bool corr_ss = true;
    if(corr_ss){
        //apply DY samesign corrections
        char fn_ss_dy[100];
        sprintf(fn_ss_dy, "../analyze/SFs/%i/dy_ss_rw.root", year);

        TFile *f = new TFile(fn_ss_dy, "READ");
        f->cd();
        TH1F *h_dy_ss_corrs = (TH1F *)gDirectory->Get("h_ss_ratio");
        TH1F *h_dy_ss_mlow_cost = (TH1F *)gDirectory->Get("DY_ss_mlow_cost");
        h_dy_ss_mlow_cost->Scale(1./h_dy_ss_mlow_cost->Integral());
        float pre_scale = DY_cost->Integral();
        for(int i=1; i<= n_cost_bins; i++){
            float corr = h_dy_ss_corrs->GetBinContent(i);
            float corr_stat_err = h_dy_ss_corrs->GetBinError(i);
            float cont = DY_cost->GetBinContent(i);
            float new_cont = cont * corr;
            float err = DY_cost->GetBinError(i);
            float scaled_cont = cont / pre_scale;
            float scaled_mlow = h_dy_ss_mlow_cost->GetBinContent(i);
            //fraction difference in distributions is a sys error
            float corr_sys_err = min(0.5, abs(scaled_cont - scaled_mlow)/(0.5 * (scaled_mlow + scaled_cont)));
            //corr_sys_err = 0.;

            printf("scaled_cont, scaled_mlow: %.3f, %.3f \n", scaled_cont, scaled_mlow);
            printf("err, corr_stat_err, corr_sys_err: %.3f %.3f %.3f \n", err, corr_stat_err, corr_sys_err);

            float new_err = sqrt(err*err + (corr_stat_err*corr_stat_err + corr_sys_err*corr_sys_err) *new_cont*new_cont);
            DY_cost->SetBinContent(i, new_cont);
            DY_cost->SetBinError(i, new_err);

        }
        float post_scale = DY_cost->Integral();
        float dy_scale = post_scale/pre_scale;
        DY_m->Scale(dy_scale);
        DY_pt->Scale(dy_scale);
        DY_xf->Scale(dy_scale);
        DY_phi->Scale(dy_scale);
        DY_rap->Scale(dy_scale);
    }

    DY_cost->Print("range");
    



    printf("Integrals of data, QCD, back, DY are %.2f %.2f %.2f %.2f \n", data_m->Integral(), QCD_m->Integral(), back_m->Integral(), DY_m->Integral());
    printf("Phi Data %.2f %.2f %.2f \n", data_phi->Integral(), QCD_phi->Integral(), back_phi->Integral());



    bool normalize = false;
    bool from_fit = false;



    
    if(normalize){
        Double_t n_data = data_cost->Integral();
        Double_t n_mc = back_cost->Integral() +  DY_cost->Integral();
        Double_t n_QCD = QCD_cost->Integral();
        Double_t qcd_ratio = (n_data - n_mc) / n_QCD;
        printf("Ratio of obs to expected QCD is %.2f \n", qcd_ratio);


        QCD_cost->Scale(qcd_ratio);

        /*
        QCD_m->Scale(qcd_ratio);
        QCD_pt->Scale(qcd_ratio);
        QCD_xf->Scale(qcd_ratio);
        QCD_phi->Scale(qcd_ratio);
        */
    }


    bool scale_error=false;
    float qcd_err = 0.5;
    float back_err = 0.05;
    bool add_err = true;
    if(scale_error){
        setHistError(QCD_m, qcd_err, add_err);
        //setHistError(QCD_cost, qcd_err, add_err);
        setHistError(QCD_xf, qcd_err, add_err);
        setHistError(QCD_phi, qcd_err, add_err);

        setHistError(back_m, back_err, add_err);
        setHistError(back_cost, back_err, add_err);
        setHistError(back_xf, back_err, add_err);
        setHistError(back_phi, back_err, add_err);

    }







    

    THStack *m_stack = new THStack("m_stack", "ElEl Mass Distribution: Data vs MC ; m_{e^{+}e^{-}} (GeV)");
    m_stack->Add(back_m);
    m_stack->Add(QCD_m);
    m_stack->Add(DY_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; ee cos(#theta)");
    cost_stack->Add(back_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(DY_cost);

    THStack *pt_stack = new THStack("pt_stack", "Dielectron Pt Distribution: Data vs MC; Dielectron Pt (GeV)");
    pt_stack->Add(back_pt);
    pt_stack->Add(QCD_pt);
    pt_stack->Add(DY_pt);

    THStack *xf_stack = new THStack("xf_stack", "Dielectron x_F Distribution: Data vs MC; x_F");
    xf_stack->Add(back_xf);
    xf_stack->Add(QCD_xf);
    xf_stack->Add(DY_xf);

    THStack *phi_stack = new THStack("phi_stack", "Dielectron phi Distribution: Data vs MC; #phi");
    phi_stack->Add(back_phi);
    phi_stack->Add(QCD_phi);
    phi_stack->Add(DY_phi);

    THStack *rap_stack = new THStack("rap_stack", "Dimuon rap Distribution: Data vs MC; y");
    rap_stack->Add(back_rap);
    rap_stack->Add(QCD_rap);
    rap_stack->Add(DY_rap);


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.25, 0.25);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(DY_m, "DY (miss-sign)", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(back_m, "t#bar{t} + wt + WW + WZ + ZZ", "f");


    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");
    TLegend *leg4 = (TLegend *) leg1->Clone("leg4");
    TLegend *leg5 = (TLegend *) leg1->Clone("leg5");


    TCanvas *c_m, *c_cost, *c_pt, *c_xf, *c_phi, *c_rap;
    TPad *p_m, *p_cost, *p_pt, *p_xf, *p_phi, *p_rap;
    int iPeriod = 4; 
    writeExtraText = false;
    char plt_file[100];

    bool logx = false;
    bool draw_sys_uncs = false;


    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{ee} (GeV)","", -1., true, logx, draw_sys_uncs);
    CMS_lumi(p_m, year, 33 );
    sprintf(plt_file, "%sElEl%i_ss_m_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_m->Print(plt_file);

    
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "cos(#theta)","", -1., false, logx, draw_sys_uncs);
    CMS_lumi(p_cost, year, 33);
    sprintf(plt_file, "%sElEl%i_ss_cost_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_cost->Print(plt_file);

    std::tie(c_pt, p_pt) = make_stack_ratio_plot(data_pt, pt_stack, leg3, "pt", "dielectron pt (GeV)","", -1., true, logx, draw_sys_uncs);
    CMS_lumi(p_pt, year, 33);

    std::tie(c_xf, p_xf) = make_stack_ratio_plot(data_xf, xf_stack, leg4, "xf", "x_F (GeV)","", -1., true, logx, draw_sys_uncs);
    CMS_lumi(p_xf, year, 33);

    std::tie(c_phi, p_phi) = make_stack_ratio_plot(data_phi, phi_stack, leg5, "phi", "dielectron #phi","", -1., true, logx, draw_sys_uncs);
    CMS_lumi(p_phi, year, 33);

    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg5, "rap", "dimuon Y","", -1., true, logx, draw_sys_uncs);
    CMS_lumi(p_rap, year, 33);
}

    
    
