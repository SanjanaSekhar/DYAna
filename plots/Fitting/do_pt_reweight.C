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
char *out_file = "../analyze/SFs/2018/pt_rw.root";
const bool write_out = true;
char *plot_dir = "Misc_plots/pt_reweights/";




void do_pt_reweight(){



    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    gROOT->SetBatch(1);



    bool logy = true;
    char plot_file[100];
    char h_name[100];

    float qcd_err = 0.5;
    bool ss_qcd = true;
    bool in_os_region = true;

    TFile *f_out;
    if(write_out) f_out = TFile::Open(out_file, "RECREATE");

    for(int i=0; i< n_m_bins-1; i++){

        TH1F *data_mumu_pt = new TH1F("data_mumu_pt", "pt reweight", n_pt_bins, pt_bins);
        TH1F *mc_mumu_pt = new TH1F("mc_mumu_pt", "pt reweight", n_pt_bins, pt_bins);
        TH1F *bkg_mumu_pt = new TH1F("bkg_mumu_pt", "pt reweight", n_pt_bins, pt_bins);
        TH1F *QCD_mumu_pt = new TH1F("QCD_mumu_pt", "pt reweight", n_pt_bins, pt_bins);

        TH1F *data_elel_pt = new TH1F("data_elel_pt", "pt reweight", n_pt_bins, pt_bins);
        TH1F *mc_elel_pt = new TH1F("mc_elel_pt", "pt reweight", n_pt_bins, pt_bins);
        TH1F *bkg_elel_pt = new TH1F("bkg_elel_pt", "pt reweight", n_pt_bins, pt_bins);
        TH1F *QCD_elel_pt = new TH1F("QCD_elel_pt", "pt reweight", n_pt_bins, pt_bins);

        TH1F *h_dummy = new TH1F("h_dummy", "", 100, 0, 100);


        //setup_all_SFs(year);

        float m_low = m_bins[i];
        float m_high = m_bins[i+1];
        //merge last 2 mass bins b/c low stats
        if(i == n_m_bins-2) m_high = m_bins[i+2];
        //m_low = 150.;
        //m_high = 13000.;

        make_m_cost_pt_xf_hist(t_mumu_data, h_dummy, h_dummy, data_mumu_pt, h_dummy, h_dummy, h_dummy, true, FLAG_MUONS,  year, m_low, m_high);

        //combine all samples coming from DY MC
        make_m_cost_pt_xf_hist(t_mumu_mc, h_dummy, h_dummy, mc_mumu_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_MUONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_nosig, h_dummy, h_dummy, mc_mumu_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_MUONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_tautau, h_dummy, h_dummy, mc_mumu_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_MUONS,   year, m_low, m_high);

        //all backgrounds
        make_m_cost_pt_xf_hist(t_mumu_ttbar, h_dummy, h_dummy, bkg_mumu_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_MUONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_wt, h_dummy, h_dummy, bkg_mumu_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_MUONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_gamgam, h_dummy, h_dummy, bkg_mumu_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_MUONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_diboson, h_dummy, h_dummy, bkg_mumu_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_MUONS,   year, m_low, m_high);

        make_fakerate_est(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_dummy, h_dummy, QCD_mumu_pt, h_dummy, h_dummy, h_dummy, FLAG_MUONS, year, m_low, m_high, ss_qcd, in_os_region);
        setHistError(QCD_mumu_pt, qcd_err);
        bkg_mumu_pt->Add(QCD_mumu_pt);


        //QCD_cost->Scale(0.);

        sprintf(h_name, "mumu%i_m%i_pt_data_sub", year % 2000, i);
        TH1F *h_mumu_data_sub = (TH1F *) data_mumu_pt->Clone(h_name);
        h_mumu_data_sub->Add(bkg_mumu_pt, -1);


        printf("Data integral is %.2f \n", data_mumu_pt->Integral());
        printf("Data subtracted integral is %.2f \n", h_mumu_data_sub->Integral());
        printf("DY integral is %.2f \n", mc_mumu_pt->Integral());

        h_mumu_data_sub->Scale(1./h_mumu_data_sub->Integral());
        mc_mumu_pt->Scale(1./mc_mumu_pt->Integral());


        TCanvas *c_mumu_plot = make_ratio_plot(string("pt_comparison"), h_mumu_data_sub, "Data - Backgrounds", mc_mumu_pt, "DY MC", "ratio", "dimuon pT", logy, false);
        sprintf(plot_file, "%sy%i_m%.0f_mumu_pt_rw.png", plot_dir, year - 2000,  m_low);
        c_mumu_plot->Print(plot_file);

        sprintf(h_name, "mumu%i_m%i_pt_ratio", year % 2000, i);
        TH1F *h_mumu_ratio = (TH1F *) h_mumu_data_sub->Clone(h_name);
        h_mumu_ratio->Divide(mc_mumu_pt);
        h_mumu_ratio->Print("range");

       

        h_dummy->Reset();
        make_m_cost_pt_xf_hist(t_elel_data, h_dummy, h_dummy, data_elel_pt, h_dummy, h_dummy, h_dummy, true, FLAG_ELECTRONS,  year, m_low, m_high);

        //combine all samples coming from DY MC
        make_m_cost_pt_xf_hist(t_elel_mc, h_dummy, h_dummy, mc_elel_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_ELECTRONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_elel_nosig, h_dummy, h_dummy, mc_elel_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_ELECTRONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_elel_tautau, h_dummy, h_dummy, mc_elel_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_ELECTRONS,   year, m_low, m_high);

        //all backgrounds
        make_m_cost_pt_xf_hist(t_elel_ttbar, h_dummy, h_dummy, bkg_elel_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_ELECTRONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_elel_wt, h_dummy, h_dummy, bkg_elel_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_ELECTRONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_elel_gamgam, h_dummy, h_dummy, bkg_elel_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_ELECTRONS,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_elel_diboson, h_dummy, h_dummy, bkg_elel_pt, h_dummy, h_dummy, h_dummy,  false, FLAG_ELECTRONS,   year, m_low, m_high);

        make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_dummy, h_dummy, QCD_elel_pt, h_dummy, h_dummy, h_dummy, FLAG_ELECTRONS, year, m_low, m_high, ss_qcd, in_os_region);
        setHistError(QCD_elel_pt, qcd_err);
        bkg_elel_pt->Add(QCD_elel_pt);


        //QCD_cost->Scale(0.);

        sprintf(h_name, "elel%i_m%i_pt_data_sub", year % 2000, i);
        TH1F *h_elel_data_sub = (TH1F *) data_elel_pt->Clone(h_name);
        h_elel_data_sub->Add(bkg_elel_pt, -1);


        printf("Data integral is %.2f \n", data_elel_pt->Integral());
        printf("Data subtracted integral is %.2f \n", h_elel_data_sub->Integral());
        printf("DY integral is %.2f \n", mc_elel_pt->Integral());

        h_elel_data_sub->Scale(1./h_elel_data_sub->Integral());
        mc_elel_pt->Scale(1./mc_elel_pt->Integral());


        TCanvas *c_elel_plot = make_ratio_plot(string("pt_comparison"), h_elel_data_sub, "Data - Backgrounds", mc_elel_pt, "DY MC", "ratio", "electron pT", logy, false);
        sprintf(plot_file, "%sy%i_m%.0f_elel_pt_rw.png", plot_dir, year - 2000,  m_low);
        c_elel_plot->Print(plot_file);

        sprintf(h_name, "elel%i_m%i_pt_ratio", year % 2000, i);
        TH1F *h_elel_ratio = (TH1F *) h_elel_data_sub->Clone(h_name);
        h_elel_ratio->Divide(mc_elel_pt);
        h_elel_ratio->Print("range");

        if(write_out){
            f_out->cd();
            h_mumu_ratio->Write();
            h_elel_ratio->Write();
            h_elel_data_sub->Write();
            h_mumu_data_sub->Write();
        }
    }


}



