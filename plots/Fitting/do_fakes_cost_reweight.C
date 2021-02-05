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



const int year = 2018;
char *out_file = "../analyze/SFs/2016/fakes_cost_rw.root";
const bool write_out = false;
bool corr_ss = true;
char *plot_dir = "Misc_plots/fakes_cost_reweights_test/";




void do_fakes_cost_reweight(){



    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    gROOT->SetBatch(1);



    bool logy = false;
    char plot_file[100];
    char h_name[100];
    char plt_title[100];

    TFile *f_out;
    if(write_out) f_out = TFile::Open(out_file, "RECREATE");

    int n_bins = 8;

    TH1F *mumu_data_cost = new TH1F("mumu_data_cost", "Data", n_bins, -1.,1.);
    TH1F *mumu_other_cost = new TH1F("mumu_diboson_cost", "DiBoson (WW, WZ,ZZ)", n_bins, -1.,1);
    TH1F *mumu_QCD_ss_cost = new TH1F("mumu_QCD_ss_cost", "QCD", n_bins, -1.,1);
    TH1F *mumu_QCD_os_cost = new TH1F("mumu_QCD_os_cost", "QCD", n_bins, -1.,1);

    TH1F *elel_data_cost = new TH1F("elel_data_cost", "Data", n_bins, -1.,1.);
    TH1F *elel_other_cost = new TH1F("elel_diboson_cost", "DiBoson (WW, WZ,ZZ)", n_bins, -1.,1);
    TH1F *elel_dy_ss_cost = new TH1F("elel_dy_ss_cost", "DiBoson (WW, WZ,ZZ)", n_bins, -1.,1);
    TH1F *elel_QCD_ss_cost = new TH1F("elel_QCD_ss_cost", "QCD", n_bins, -1.,1);
    TH1F *elel_QCD_os_cost = new TH1F("elel_QCD_os_cost", "QCD", n_bins, -1.,1);

    TH1F *dummy = new TH1F("h_dummy", "", 100, 0, 100);

    int m_low = 170.;
    int m_high = 10000.;
    bool ss = true;

    make_m_cost_pt_xf_hist(t_mumu_ss_data, dummy, mumu_data_cost, dummy, dummy, dummy,  dummy, true, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_diboson, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_wt, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_ttbar, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_dy, dummy, mumu_other_cost, dummy, dummy, dummy, dummy, false, FLAG_MUONS,   year, m_low, m_high, ss);

    bool reweight = false;
    bool sys_errors = false;
    make_fakerate_est(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, dummy, mumu_QCD_ss_cost, dummy, dummy, dummy, dummy, FLAG_MUONS, 
            year, m_low, m_high, ss, reweight, sys_errors);
    ss = false;
    make_fakerate_est(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, dummy, mumu_QCD_os_cost, dummy, dummy, dummy, dummy, FLAG_MUONS, 
            year, m_low, m_high, ss, reweight, sys_errors);

    symmetrize1d(mumu_QCD_os_cost);

    sprintf(h_name, "mumu%i_ss_cost_data_sub", year % 2000);
    TH1F *h_mumu_data_sub = (TH1F *) mumu_data_cost->Clone(h_name);
    h_mumu_data_sub->Add(mumu_other_cost, -1);
    symmetrize1d(h_mumu_data_sub);


    sprintf(plt_title, "Muons %i", year);
    TCanvas *c_mumu_plot = make_ratio_plot(string("mumu_ss_cost_comparison"), h_mumu_data_sub, "Data - Other Backgrounds", mumu_QCD_ss_cost, "Fakes Estimate", 
            "ratio", "samesign mumu cos(#theta)", logy, false, 0.0, 2.0, plt_title);
    sprintf(plot_file, "%sy%i_mumu_ss_cost_rw.png", plot_dir, year - 2000);
    c_mumu_plot->Print(plot_file);

    sprintf(h_name, "mumu%i_ss_cost_ratio", year % 2000);
    TH1F *h_mumu_ratio = (TH1F *) h_mumu_data_sub->Clone(h_name);
    h_mumu_ratio->Divide(mumu_QCD_ss_cost);

    // Do systematic on shape from os to ss fakes est ratio
    mumu_QCD_ss_cost->Scale(1./mumu_QCD_ss_cost->Integral());
    mumu_QCD_os_cost->Scale(1./mumu_QCD_os_cost->Integral());

    TCanvas *c_os_ss_mumu_plot = make_ratio_plot(string("mumu_ss_os_comparison"), mumu_QCD_os_cost, "OS Fakes Estimate", mumu_QCD_ss_cost, "SS Fakes Estimate", 
            "OS/SS", "mumu cos(#theta)", logy, false, 0.0,2.0, plt_title);
    sprintf(plot_file, "%sy%i_mumu_os_ss_cost_ratio.png", plot_dir, year - 2000);
    c_os_ss_mumu_plot->Print(plot_file);


    h_mumu_ratio->Print("range");
    printf("Adding sys errors from os/ss ratio \n");
    for(int bin = 1; bin<=n_bins; bin++){
        
        float corr = h_mumu_ratio->GetBinContent(bin);
        float ss_cont  = mumu_QCD_ss_cost->GetBinContent(bin);
        float os_cont  = mumu_QCD_os_cost->GetBinContent(bin);
        float sys = (fabs(os_cont - ss_cont) / os_cont) * fabs(1. - corr);
        float stat = h_mumu_ratio->GetBinError(bin);
        float tot_err = pow(sys*sys + stat*stat, 0.5);
        tot_err = std::fmin(tot_err, 0.7*corr);
        h_mumu_ratio->SetBinError(bin, tot_err);
    }
    h_mumu_ratio->Print("range");


    dummy->Reset();

    make_m_cost_pt_xf_hist(t_elel_ss_data, dummy, elel_data_cost, dummy, dummy, dummy,  dummy, true, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_diboson, dummy, elel_other_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_wt, dummy, elel_other_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_ttbar, dummy, elel_other_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_dy, dummy, elel_dy_ss_cost, dummy, dummy, dummy, dummy, false, FLAG_ELECTRONS,   year, m_low, m_high, ss);

    ss = true;
    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, dummy, elel_QCD_ss_cost, dummy, dummy, dummy, dummy, FLAG_ELECTRONS, 
            year, m_low, m_high, ss, reweight, sys_errors);
    ss = false;
    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, dummy, elel_QCD_os_cost, dummy, dummy, dummy, dummy, FLAG_ELECTRONS, 
            year, m_low, m_high, ss, reweight, sys_errors);
    symmetrize1d(elel_QCD_os_cost);



    if(corr_ss){
        //apply DY samesign corrections
        char fn_ss_dy[100];
        sprintf(fn_ss_dy, "../analyze/SFs/%i/dy_ss_rw.root", year);

        TFile *f = new TFile(fn_ss_dy, "READ");
        f->cd();
        TH1F *h_dy_ss_corrs = (TH1F *)gDirectory->Get("h_ss_ratio");
        for(int i=1; i<= n_cost_bins; i++){
            float corr = h_dy_ss_corrs->GetBinContent(i);
            float corr_err = h_dy_ss_corrs->GetBinError(i);
            float cont = elel_dy_ss_cost->GetBinContent(i);
            float err = elel_dy_ss_cost->GetBinError(i);
            float new_err = sqrt(err*err + corr_err*cont*corr_err*cont);
            elel_dy_ss_cost->SetBinContent(i, cont*corr);
            elel_dy_ss_cost->SetBinError(i, new_err);

        }
    }



    sprintf(h_name, "elel%i_ss_cost_data_sub", year % 2000);
    
    TH1F *h_elel_data_sub = (TH1F *) elel_data_cost->Clone(h_name);
    h_elel_data_sub->Add(elel_other_cost, -1);
    h_elel_data_sub->Add(elel_dy_ss_cost, -1);
    symmetrize1d(h_elel_data_sub);

    sprintf(plt_title, "Electrons %i", year);
    TCanvas *c_elel_plot = make_ratio_plot(string("elel_ss_cost_comparison"), h_elel_data_sub, "Data - Other Backgrounds", elel_QCD_ss_cost, "Fakes Estimate", 
            "ratio", "samesign ee cos(#theta)", logy, false, 0.0, 2.0, plt_title);
    sprintf(plot_file, "%sy%i_elel_ss_cost_rw.png", plot_dir, year - 2000);
    c_elel_plot->Print(plot_file);

    sprintf(h_name, "elel%i_ss_cost_ratio", year % 2000);
    TH1F *h_elel_ratio = (TH1F *) h_elel_data_sub->Clone(h_name);
    h_elel_ratio->Divide(elel_QCD_ss_cost);

    // Do systematic on shape from os to ss fakes est ratio
    elel_QCD_ss_cost->Scale(1./elel_QCD_ss_cost->Integral());
    elel_QCD_os_cost->Scale(1./elel_QCD_os_cost->Integral());

    TCanvas *c_os_ss_elel_plot = make_ratio_plot(string("elel_ss_os_comparison"), elel_QCD_os_cost, "OS Fakes Estimate", elel_QCD_ss_cost, "SS Fakes Estimate", 
            "OS/SS", "elel cos(#theta)", logy, false, 0.0, 2.0, plt_title);
    sprintf(plot_file, "%sy%i_elel_os_ss_cost_ratio.png", plot_dir, year - 2000);
    c_os_ss_elel_plot->Print(plot_file);


    h_elel_ratio->Print("range");
    printf("Adding sys errors from os/ss ratio \n");
    for(int bin = 1; bin<=n_bins; bin++){
        float corr = h_elel_ratio->GetBinContent(bin);
        float ss_cont  = elel_QCD_ss_cost->GetBinContent(bin);
        float os_cont  = elel_QCD_os_cost->GetBinContent(bin);
        float sys = (fabs(os_cont - ss_cont) / os_cont) * fabs(1. - corr);
        float stat = h_elel_ratio->GetBinError(bin);
        float tot = pow(sys*sys + stat*stat, 0.5);
        h_elel_ratio->SetBinError(bin, tot);
    }
    h_elel_ratio->Print("range");


    if(write_out){
        f_out->cd();
        printf("Writing out \n");
        h_mumu_ratio->Write();
        h_elel_ratio->Write();
        h_elel_data_sub->Write();
        h_mumu_data_sub->Write();
    }
}
