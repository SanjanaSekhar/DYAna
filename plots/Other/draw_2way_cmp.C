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







void draw_2way_cmp(){


    const int type = FLAG_MUONS;
    const int year = 2016;
    const bool write_out = true;
    char *plot_dir = "Paper_plots/photon_sample_cmp/";
    char *file_label = "MuMu16_Inel_only";
    char *plt_label = "Muons: 2016";

    bool do_emu = false;
    bool normalize = false;



    //char *f1n= "../analyze/output_files/2016/MuMu16_phot_ind_mar25.root";
    char *f1n= "../analyze/output_files/test.root";
    char *f2n= "../analyze/output_files/MuMu16_phot_ind_LPair_Inel_mar29.root";

    char *name1 = "Pythia";
    char *name2 = "LPair (Inelastic Only)" ;



    
    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    setup_all_SFs(year);





    TFile * f1 = TFile::Open(f1n, "READ");
    TTree *t1 = (TTree *) f1->Get("T_sig");

    TFile * f2 = TFile::Open(f2n, "READ");
    TTree *t2 = (TTree *) f2->Get("T_sig");



    int n_m_bins = 61;
    float m_bin_size = 30.;
	float m_bin_low = 170.;
	float m_bin_high = m_bin_low + n_m_bins*m_bin_size;
    TH1F *m1 = new TH1F("m1", "Data Dimuon Mass Distribution", n_m_bins, m_bin_low, m_bin_high);
    TH1F *m2 = new TH1F("m2", "Data Dimuon Mass Distribution", n_m_bins, m_bin_low, m_bin_high);

    int n_cost_bins = 16;
    float cost_bin_size = 2./n_cost_bins;
    TH1F *cost1 = new TH1F("cost1", "Data", n_cost_bins, -1.,1.);
    TH1F *cost2 = new TH1F("cost2", "Data", n_cost_bins, -1.,1.);

    int n_phi_bins = 20;
    TH1F *phi1 = new TH1F("phi1", "Data", n_phi_bins, -4.,4.);
    TH1F *phi2 = new TH1F("phi2", "Data", n_phi_bins, -4.,4.);

    int n_rap_bins = 20;
    float rap_bin_size = 5. / n_rap_bins;
    TH1F *rap1 = new TH1F("rap1", "Data", n_rap_bins, -2.5,2.5);
    TH1F *rap2 = new TH1F("rap2", "Data", n_rap_bins, -2.5,2.5);


    int n_xf_bins1 = 5;
    float xf_bins1[] = {0.,0.04, 0.07, 0.1, 0.2, 0.5};
    TH1F *xf1 = new TH1F("xf1", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *xf2 = new TH1F("xf2", "dy signal", n_xf_bins1,  xf_bins1);


    int n_pt_bins1 = 7;
    Float_t pt_bins1[] = {0., 10., 20., 30., 50., 70., 100., 300., 700. };
    //TH1F *pt1 = new TH1F("pt1", "dy signal", n_pt_bins1, pt_bins1);
    //TH1F *pt2 = new TH1F("pt2", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *pt1 = new TH1F("pt1", "", 10, 0, 400);
    TH1F *pt2 = new TH1F("pt2", "", 10, 0, 400);



    float m_low = 170.;
    float m_high = 13000;
    bool ss = false;
    bool logy;
    float ratio_min = 0.5;
    float ratio_max = 1.5;


    char plt_file1[100];

    if(!do_emu){
        make_m_cost_pt_xf_hist(t1, m1, cost1, pt1, xf1, phi1, rap1, false, type, year, m_low, m_high);
        make_m_cost_pt_xf_hist(t2, m2, cost2, pt2, xf2, phi2, rap2, false, type, year, m_low, m_high);
    }
    else{
        bool do_emu_cost_rw = false;

        make_emu_m_cost_pt_rap_hist(t1, m1, cost1, pt1, rap1, false,  year, m_low, m_high, ss, do_emu_cost_rw);
        make_emu_m_cost_pt_rap_hist(t2, m2, cost2, pt2, rap2, false,  year, m_low, m_high, ss, do_emu_cost_rw);
    }

    if(normalize){
        m1->Scale(1./m1->Integral());
        m2->Scale(1./m2->Integral());

        cost1->Scale(1./cost1->Integral());
        cost2->Scale(1./cost2->Integral());

        pt1->Scale(1./pt1->Integral());
        pt2->Scale(1./pt2->Integral());

        rap1->Scale(1./rap1->Integral());
        rap2->Scale(1./rap2->Integral());
    }


    logy = true;
    sprintf(plt_file1, "%s%s_m_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, m1, name1, m2, name2, "Ratio", "M (GeV)", logy, write_out, ratio_min, ratio_max, plt_label);

    logy = false;
    sprintf(plt_file1, "%s%s_cost_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, cost1, name1, cost2, name2, "Ratio", "cos(#theta)", logy, write_out, ratio_min, ratio_max, plt_label);


    logy = true;
    sprintf(plt_file1, "%s%s_pt_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, pt1, name1, pt2, name2, "Ratio", "p_{T} (GeV)", logy, write_out, ratio_min, ratio_max, plt_label);

    logy = true;
    sprintf(plt_file1, "%s%s_rap_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, rap1, name1, rap2, name2, "Ratio", "Rapidity", logy, write_out, ratio_min, ratio_max, plt_label);

}

    
    
