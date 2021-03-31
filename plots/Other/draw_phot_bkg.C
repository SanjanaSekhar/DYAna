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







void draw_phot_bkg(){


    const int type = FLAG_MUONS;
    const int year = 2018;
    const bool write_out = true;
    char *plot_dir = "Misc_plots/photon_bkg/";
    char *file_label = "MuMu18";
    char *plt_label = "Photon Induced Bkg, Muons, 2018";

    bool do_emu = false;



    char *f1n= "../analyze/output_files/2018/MuMu18_phot_ind_nov11.root";



    
    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    setup_all_SFs(year);





    TFile * f1 = TFile::Open(f1n, "READ");
    TTree *t1 = (TTree *) f1->Get("T_sig");



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
    TH1F *pt1 = new TH1F("pt1", "dy signal", 10, 0, 200);



    float m_low = 170.;
    float m_high = 13000;
    bool ss = false;
    bool logy;
    float ratio_min = 0.5;
    float ratio_max = 1.5;


    char plt_file1[100];

    make_m_cost_pt_xf_hist(t1, m1, cost1, pt1, xf1, phi1, rap1, false, type, year, m_low, m_high);


    logy = true;
    TCanvas *c1 = new TCanvas("c1", "", 800,800);
    c1->SetLogy();
    m1->SetLineColor(kBlue);
    m1->Draw("hist");
    m1->GetXaxis()->SetTitle("M (GeV)");
    sprintf(plt_file1, "%s%s_m_cmp.pdf", plot_dir, file_label);
    c1->Print(plt_file1);

    TCanvas *c2 = new TCanvas("c2", "", 800,800);
    m1->SetLineColor(kBlue);
    cost1->Draw("hist");
    cost1->GetXaxis()->SetTitle("cos(#theta)");
    sprintf(plt_file1, "%s%s_cost_cmp.pdf", plot_dir, file_label);
    c2->Print(plt_file1);



    TCanvas *c3 = new TCanvas("c3", "", 800,800);
    c3->SetLogy();
    pt1->SetLineColor(kBlue);
    pt1->Draw("hist");
    pt1->GetXaxis()->SetTitle("dilepton p_T (GeV)");
    sprintf(plt_file1, "%s%s_pt_cmp.pdf", plot_dir, file_label);
    c3->Print(plt_file1);


    TCanvas *c4 = new TCanvas("c4", "", 800,800);
    c4->SetLogy();
    rap1->SetLineColor(kBlue);
    rap1->Draw("hist");
    rap1->GetXaxis()->SetTitle("dilepton Y");
    sprintf(plt_file1, "%s%s_rap_cmp.pdf", plot_dir, file_label);
    c4->Print(plt_file1);


}

    
    
