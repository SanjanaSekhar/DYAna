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
//char *plot_dir = "Misc_plots/emu_cost_reweights_test/";
char *plot_dir = "Paper_plots/emu_cost_reweights/";




void do_emu_cost_reweight(){



    setTDRStyle();
    init_emu(year);
    init_emu_indv_bkgs(year);
    gROOT->SetBatch(1);



    bool logy = false;
    char plot_file[100];
    char h_name[100];

    TFile *f_out;
    if(write_out) f_out = TFile::Open(out_file, "RECREATE");

    int n_bins = 8;

    TH1F *emu_data_cost = new TH1F("emu_data_cost", "Data", n_bins, -1.,1.);
    TH1F *emu_bkg_cost = new TH1F("emu_bkg_cost", "DiBoson ttbar and wt", n_bins, -1.,1);
    TH1F *emu_dy_cost = new TH1F("emu_dy_cost", "DY", n_bins, -1.,1);
    TH1F *emu_QCD_cost = new TH1F("emu_QCD_cost", "QCD", n_bins, -1.,1);

    TH1F *dummy = new TH1F("h_dummy", "", 100, 0, 100);


    bool ss = false;

    for(int mbin = 0; mbin < n_emu_rw_m_bins; mbin++){
        emu_data_cost->Reset();
        emu_bkg_cost->Reset();
        emu_dy_cost->Reset();
        emu_QCD_cost->Reset();

        float m_low = emu_rw_m_bins[mbin];
        float m_high = emu_rw_m_bins[mbin+1];



        make_emu_m_cost_pt_rap_hist(t_emu_data, dummy, emu_data_cost, dummy,   dummy, true,   year, m_low, m_high, ss);
        make_emu_m_cost_pt_rap_hist(t_emu_ttbar, dummy, emu_bkg_cost, dummy,  dummy, false,  year, m_low, m_high, ss);
        make_emu_m_cost_pt_rap_hist(t_emu_diboson, dummy, emu_bkg_cost, dummy, dummy, false,  year, m_low, m_high, ss);
        make_emu_m_cost_pt_rap_hist(t_emu_wt, dummy, emu_bkg_cost, dummy,dummy, false,   year, m_low, m_high, ss);
        make_emu_m_cost_pt_rap_hist(t_emu_dy, dummy, emu_dy_cost, dummy, dummy, false,   year, m_low, m_high, ss);


        bool fakes_reweight = false;
        bool sys_errors = true;
        dummy->Reset();
        Fakerate_est_emu(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam, t_emu_QCD_contam, dummy, emu_QCD_cost, dummy, dummy, FLAG_MUONS, year, m_low, m_high, fakes_reweight, sys_errors);

        float qcd_err = 0.5;
        setHistError(emu_QCD_cost, qcd_err);


        symmetrize1d(emu_bkg_cost);
        symmetrize1d(emu_QCD_cost);



        char plt_title[200];
        sprintf(plt_title, "%i: M %.0f-%0.f GeV", year, m_low, m_high);


        sprintf(h_name, "emu%i_mbin%i_cost_data_sub", year % 2000, mbin);
        TH1F *h_emu_data_sub = (TH1F *) emu_data_cost->Clone(h_name);
        h_emu_data_sub->Add(emu_dy_cost, -1);
        h_emu_data_sub->Add(emu_QCD_cost, -1);

        symmetrize1d(h_emu_data_sub);




        TCanvas *c_emu_plot = make_ratio_plot(string("emu_cost_comparison"), h_emu_data_sub, "Data - Other Backgrounds", emu_bkg_cost, "ttbar + tW + diboson MC Estimate", "ratio", "e#mu cos(#theta)", logy, false, 0.5, 1.5, plt_title);
        sprintf(plot_file, "%sy%i_m%i_emu_cost_rw.png", plot_dir, year - 2000, mbin);
        c_emu_plot->Print(plot_file);

        sprintf(h_name, "emu%i_mbin%i_cost_ratio", year % 2000, mbin);
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
}
