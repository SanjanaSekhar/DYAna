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

int year = 2016;
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

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, 0, 100.);

    Double_t m_low = 150;
    Double_t m_high = 10000;

    make_emu_m_cost_pt_xf_hist(t_emu_data, data_m, h_dummy, h_dummy, h_dummy, true,  year, m_low, m_high);
    make_emu_m_cost_pt_xf_hist(t_emu_ttbar, ttbar_m, h_dummy, h_dummy, h_dummy, false,  year, m_low, m_high);
    make_emu_m_cost_pt_xf_hist(t_emu_diboson, diboson_m, h_dummy, h_dummy, h_dummy, false,  year, m_low, m_high);
    make_emu_m_cost_pt_xf_hist(t_emu_wt, wt_m, h_dummy, h_dummy, h_dummy, false,  year, m_low, m_high);
    make_emu_m_cost_pt_xf_hist(t_emu_dy, dy_m, h_dummy, h_dummy, h_dummy, false,  year, m_low, m_high);

    Fakerate_est_emu(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam, qcd_m,  FLAG_MUONS, year, m_low, m_high);

    //correct for wrong ttbar xsec
    //ttbar_m->Scale(831.76/730.6);

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

    THStack *m_stack = new THStack("m_stack", "EMu Mass Distribution: Data vs MC ; m_{e#mu} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(qcd_m);
    m_stack->Add(dy_m);
    m_stack->Add(wt_m);
    m_stack->Add(ttbar_m);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(dy_m, "DY #rightarrow #tau#tau", "f");
    leg1->AddEntry(qcd_m, "QCD and W+Jets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");

    int iPeriod = 4; 
    writeExtraText = false;
    
    TCanvas *c_m;
    TPad *p_m;
    char plt_file[100];
    sprintf(plt_file, "%sEMu%i_m_cmp.pdf", plot_dir, year % 2000);

    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{e#mu} (GeV)", -1., true);
    CMS_lumi(p_m, year, 33 );
    if(write_out) c_m->Print(plt_file);



}









