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
int year = 2018;
char *out_file = "../analyze/SFs/2018/dy_ss_rw.root";
bool write_out = false;
char *plot_dir = "Misc_plots/samesign_cmp_mlow/";
//char *plot_dir = "Paper_plots/";


void draw_samesign_mlow_cmp(){
    init(year);
    init_mlow(year);

    gROOT->SetBatch(1);

    TFile *f_out;
    if(write_out) f_out = TFile::Open(out_file, "RECREATE");

    setTDRStyle();

    int n_m_bins = 30;
    float m_low = 70.;
    float m_high = 110.;
    TH1F *data_m = new TH1F("data_m", "Data Dielectron Mass Distribution", n_m_bins, m_low, m_high);
    TH1F *back_m = new TH1F("back_m", "back (WW, WZ, ZZ)", n_m_bins, m_low, m_high);
    TH1F *DY_m = new TH1F("DY_m", "DY (WW, WZ, ZZ)", n_m_bins, m_low, m_high);


    int pt_bins = 20.;
    TH1F *data_pt = new TH1F("data_pt", "MC signal", pt_bins, 0, 800);
    TH1F *back_pt = new TH1F("back_pt", "MC signal", pt_bins, 0, 800);
    TH1F *DY_pt = new TH1F("DY_pt", "MC signal", pt_bins, 0, 800);

    int xf_nbins = 16;
    TH1F *data_xf = new TH1F("data_xf", "MC signal", xf_nbins, 0, 0.8);


    int n_rap_bins = 20;
    TH1F *data_rap = new TH1F("data_rap", "Data", n_rap_bins, -2.5,2.5);
    TH1F *DY_rap = new TH1F("mc_rap", "MC Signal (qqbar, qglu, qbarglu)", n_rap_bins, -2.5,2.5);
    TH1F *back_rap = new TH1F("diboson_rap", "DiBoson (WW, WZ,ZZ)", n_rap_bins, -2.5,2.5);
    TH1F *wt_rap = new TH1F("wt_rap", "tw + #bar{t}w", n_rap_bins, -2.5,2.5);


    TH1F *back_xf = new TH1F("back_xf", "MC signal", xf_nbins, 0, 0.8);

    TH1F *DY_xf = new TH1F("DY_xf", "MC signal", xf_nbins, 0, 0.8);



    int n_cost_bins = 8;
    TH1F *data_cost = new TH1F("data_cost", "Data", n_cost_bins, -1.,1.);
    TH1F *back_cost = new TH1F("back_cost", "back (WW, WZ,ZZ)", n_cost_bins, -1.,1);
    TH1F *DY_cost = new TH1F("DY_ss_mlow_cost", "DY (WW, WZ,ZZ)", n_cost_bins, -1.,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1.,1);


    TH1F *DY_phi = new TH1F("DY_phi", "", 20, -4.0, 4.0);
    TH1F *data_phi = new TH1F("data_phi", "", 20, -4.0, 4.0);
    TH1F *back_phi = new TH1F("back_phi", "", 20, -4.0, 4.0);

    TH1F * dummy = new TH1F("dummy", "", 100, 0., 100.);

    back_m->SetFillColor(diboson_c);
    back_cost->SetFillColor(diboson_c);
    back_xf->SetFillColor(diboson_c);
    back_pt->SetFillColor(diboson_c);
    back_phi->SetFillColor(diboson_c);
    back_rap->SetFillColor(diboson_c);


    DY_xf->SetFillColor(DY_c);
    DY_pt->SetFillColor(DY_c);
    DY_m->SetFillColor(DY_c);
    DY_cost->SetFillColor(DY_c);
    DY_phi->SetFillColor(DY_c);
    DY_rap->SetFillColor(DY_c);

    bool ss = true;

    make_m_cost_pt_xf_hist(t_elel_ss_data_mlow, data_m, data_cost, data_pt, data_xf, data_phi, data_rap,  true, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_ttbar, back_m, back_cost, back_pt, back_xf, back_phi, back_rap, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_wt, back_m, back_cost, back_pt, back_xf, back_phi, back_rap, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_diboson, back_m, back_cost, back_pt, back_xf, back_phi, back_rap, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_dy_mlow, DY_m, DY_cost, DY_pt, DY_xf, DY_phi, DY_rap, false, type,   year, m_low, m_high, ss);



    bool cost_reweight = false;


    printf("Integrals of data, back, DY are %.2f %.2f %.2f \n", data_m->Integral(), back_m->Integral(), DY_m->Integral());



    bool normalize = true;
    bool from_fit = false;
    





    

    THStack *m_stack = new THStack("m_stack", "ElEl Mass Distribution: Data vs MC ; m_{e^{+}e^{-}} (GeV)");
    m_stack->Add(back_m);
    m_stack->Add(DY_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; ee cos(#theta)");
    cost_stack->Add(back_cost);
    cost_stack->Add(DY_cost);

    THStack *pt_stack = new THStack("pt_stack", "Dielectron Pt Distribution: Data vs MC; Dielectron Pt (GeV)");
    pt_stack->Add(back_pt);
    pt_stack->Add(DY_pt);

    THStack *xf_stack = new THStack("xf_stack", "Dielectron x_F Distribution: Data vs MC; x_F");
    xf_stack->Add(back_xf);
    xf_stack->Add(DY_xf);

    THStack *phi_stack = new THStack("phi_stack", "Dielectron phi Distribution: Data vs MC; #phi");
    phi_stack->Add(back_phi);
    phi_stack->Add(DY_phi);

    THStack *rap_stack = new THStack("rap_stack", "Dimuon rap Distribution: Data vs MC; y");
    rap_stack->Add(back_rap);
    rap_stack->Add(DY_rap);


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.25, 0.25);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(DY_m, "DY (miss-sign)", "f");
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


    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{ee} (GeV)","", -1., true);
    CMS_lumi(p_m, year, 33 );
    sprintf(plt_file, "%sElEl%i_ss_m_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_m->Print(plt_file);

    
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "cos(#theta)","", -1., false);
    CMS_lumi(p_cost, year, 33);
    sprintf(plt_file, "%sElEl%i_ss_cost_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_cost->Print(plt_file);

    std::tie(c_pt, p_pt) = make_stack_ratio_plot(data_pt, pt_stack, leg3, "pt", "dielectron pt (GeV)","", -1., true);
    CMS_lumi(p_pt, year, 33);

    std::tie(c_xf, p_xf) = make_stack_ratio_plot(data_xf, xf_stack, leg4, "xf", "x_F (GeV)","", -1., true);
    CMS_lumi(p_xf, year, 33);

    std::tie(c_phi, p_phi) = make_stack_ratio_plot(data_phi, phi_stack, leg5, "phi", "dielectron #phi","", -1., true);
    CMS_lumi(p_phi, year, 33);

    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg5, "rap", "dimuon Y", "",  -1., true);
    CMS_lumi(p_rap, year, 33);


    if(write_out){
        f_out->cd();

        TList *stackHists = cost_stack->GetHists();
        TH1* sum = (TH1*)stackHists->At(0)->Clone();
        sum->Reset();

        for (int i=0;i<stackHists->GetSize();++i) {
          sum->Add((TH1*)stackHists->At(i));
        }
        auto h_ratio = (TH1F *) data_cost->Clone("h_ss_ratio");
        h_ratio->Sumw2();
        h_ratio->Divide(sum);
        h_ratio->Print("range");
        h_ratio->Write();

        DY_cost->Write();

    }

    

}

    
    
