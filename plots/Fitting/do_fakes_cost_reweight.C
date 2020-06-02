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



const int year = 2016;
char *out_file = "../analyze/SFs/2016/fakes_cost_rw.root";
const bool write_out = true;
char *plot_dir = "Misc_plots/fakes_cost_reweights/";




void do_fakes_cost_reweight(){



    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    gROOT->SetBatch(0);



    bool logy = false;
    char plot_file[100];
    char h_name[100];

    TFile *f_out;
    if(write_out) f_out = TFile::Open(out_file, "UPDATE");

    int n_bins = 8;

    TH1F *mumu_data_cost = new TH1F("mumu_data_cost", "Data", n_bins, -1.,1.);
    TH1F *mumu_other_cost = new TH1F("mumu_diboson_cost", "DiBoson (WW, WZ,ZZ)", n_bins, -1.,1);
    TH1F *mumu_QCD_cost = new TH1F("mumu_QCD_cost", "QCD", n_bins, -1.,1);

    TH1F *elel_data_cost = new TH1F("elel_data_cost", "Data", n_bins, -1.,1.);
    TH1F *elel_other_cost = new TH1F("elel_diboson_cost", "DiBoson (WW, WZ,ZZ)", n_bins, -1.,1);
    TH1F *elel_QCD_cost = new TH1F("elel_QCD_cost", "QCD", n_bins, -1.,1);

    TH1F *dummy = new TH1F("h_dummy", "", 100, 0, 100);

    int m_low = 150.;
    int m_high = 10000.;
    bool ss = true;
    bool in_os_region = false;

    make_m_cost_pt_xf_hist(t_mumu_ss_data, dummy, mumu_data_cost, dummy, dummy, dummy,  dummy, true, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_diboson, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_wt, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_ttbar, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_dy, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);

    make_fakerate_est(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, dummy, mumu_QCD_cost, dummy, dummy, dummy, dummy, FLAG_MUONS, year, m_low, m_high, ss, in_os_region);

    sprintf(h_name, "mumu%i_ss_cost_data_sub", year % 2000);
    TH1F *h_mumu_data_sub = (TH1F *) mumu_data_cost->Clone(h_name);
    h_mumu_data_sub->Add(mumu_other_cost, -1);


    TCanvas *c_mumu_plot = make_ratio_plot(string("mumu_ss_cost_comparison"), h_mumu_data_sub, "Data - Other Backgrounds", mumu_QCD_cost, "Fakes Estimate", "ratio", "samesign mumu cos(#theta)", logy, false);
    sprintf(plot_file, "%sy%i_mumu_ss_cost_rw.png", plot_dir, year - 2000,  m_low);
    c_mumu_plot->Print(plot_file);

    sprintf(h_name, "mumu%i_ss_cost_ratio", year % 2000);
    TH1F *h_mumu_ratio = (TH1F *) h_mumu_data_sub->Clone(h_name);
    h_mumu_ratio->Divide(mumu_QCD_cost);
    h_mumu_ratio->Print("range");


    dummy->Reset();

    make_m_cost_pt_xf_hist(t_elel_ss_data, dummy, elel_data_cost, dummy, dummy, dummy,  dummy, true, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_diboson, dummy, elel_other_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_wt, dummy, elel_other_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_ttbar, dummy, elel_other_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_dy, dummy, elel_other_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);

    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, dummy, elel_QCD_cost, dummy, dummy, dummy, dummy, FLAG_ELECTRONS, year, m_low, m_high, ss, in_os_region);

    sprintf(h_name, "elel%i_ss_cost_data_sub", year % 2000);
    TH1F *h_elel_data_sub = (TH1F *) elel_data_cost->Clone(h_name);
    h_elel_data_sub->Add(elel_other_cost, -1);


    TCanvas *c_elel_plot = make_ratio_plot(string("elel_ss_cost_comparison"), h_elel_data_sub, "Data - Other Backgrounds", elel_QCD_cost, "Fakes Estimate", "ratio", "samesign ee cos(#theta)", logy, false);
    sprintf(plot_file, "%sy%i_elel_ss_cost_rw.png", plot_dir, year - 2000,  m_low);
    c_elel_plot->Print(plot_file);

    sprintf(h_name, "elel%i_ss_cost_ratio", year % 2000);
    TH1F *h_elel_ratio = (TH1F *) h_elel_data_sub->Clone(h_name);
    h_elel_ratio->Divide(elel_QCD_cost);
    h_elel_ratio->Print("range");

    if(write_out){
        f_out->cd();
        printf("Writing out \n");
        f_out->Print();
        h_mumu_ratio->Write();
        h_elel_ratio->Write();
        h_elel_data_sub->Write();
        h_mumu_data_sub->Write();
    }
}
