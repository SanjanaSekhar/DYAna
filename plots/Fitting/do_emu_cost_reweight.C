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
char *out_file = "../analyze/SFs/2016/emu_cost_rw.root";
const bool write_out = true;
char *plot_dir = "Misc_plots/emu_cost_reweights/";




void do_emu_cost_reweight(){



    setTDRStyle();
    init_emu(year);
    init_emu_indv_bkgs(year);
    gROOT->SetBatch(0);



    bool logy = false;
    char plot_file[100];
    char h_name[100];

    TFile *f_out;
    if(write_out) f_out = TFile::Open(out_file, "UPDATE");

    int n_bins = 8;

    TH1F *emu_data_cost = new TH1F("emu_data_cost", "Data", n_bins, -1.,1.);
    TH1F *emu_bkg_cost = new TH1F("emu_bkg_cost", "DiBoson ttbar and wt", n_bins, -1.,1);
    TH1F *emu_dy_cost = new TH1F("emu_dy_cost", "DY", n_bins, -1.,1);
    TH1F *emu_QCD_cost = new TH1F("emu_QCD_cost", "QCD", n_bins, -1.,1);

    TH1F *dummy = new TH1F("h_dummy", "", 100, 0, 100);


    int m_low = 150.;
    int m_high = 10000.;

    bool ss = false;



    make_emu_m_cost_pt_rap_hist(t_emu_data, dummy, emu_data_cost, dummy,   dummy, true,   year, m_low, m_high, ss);
    make_emu_m_cost_pt_rap_hist(t_emu_ttbar, dummy, emu_bkg_cost, dummy,  dummy, false,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_rap_hist(t_emu_diboson, dummy, emu_bkg_cost, dummy, dummy, false,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_rap_hist(t_emu_wt, dummy, emu_bkg_cost, dummy,dummy, false,   year, m_low, m_high, ss);
    make_emu_m_cost_pt_rap_hist(t_emu_dy, dummy, emu_dy_cost, dummy, dummy, false,   year, m_low, m_high, ss);


    bool fakes_reweight = true;
    dummy->Reset();
    Fakerate_est_emu(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam, dummy, emu_QCD_cost, dummy, dummy, FLAG_MUONS, year, m_low, m_high, fakes_reweight);

    symmetrize1d(emu_bkg_cost);
    symmetrize1d(emu_QCD_cost);


    sprintf(h_name, "emu%i_cost_data_sub", year % 2000);
    TH1F *h_emu_data_sub = (TH1F *) emu_data_cost->Clone(h_name);
    h_emu_data_sub->Add(emu_dy_cost, -1);
    h_emu_data_sub->Add(emu_QCD_cost, -1);

    symmetrize1d(h_emu_data_sub);

    


    TCanvas *c_emu_plot = make_ratio_plot(string("emu_cost_comparison"), h_emu_data_sub, "Data - Other Backgrounds", emu_bkg_cost, "ttbar + tW + diboson MC Estimate", "ratio", "e#mu cos(#theta)", logy, false);
    sprintf(plot_file, "%sy%i_emu_cost_rw.png", plot_dir, year - 2000);
    c_emu_plot->Print(plot_file);

    sprintf(h_name, "emu%i_cost_ratio", year % 2000);
    TH1F *h_emu_ratio = (TH1F *) h_emu_data_sub->Clone(h_name);
    h_emu_ratio->Divide(emu_bkg_cost);



    if(write_out){
        f_out->cd();
        printf("Writing out \n");
        f_out->Print();
        h_emu_ratio->Write();
        h_emu_data_sub->Write();
    }
}
