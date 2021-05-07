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







void draw_phot_3sep_parts(){


    const int type = FLAG_MUONS;
    const int year = 2016;
    const bool write_out = true;
    char *plot_dir = "Misc_plots/photon_bkg/";
    char *file_label = "LPair_3sep";
    char *plt_label = "";

    bool do_emu = false;



    char *f1n= "../analyze/output_files/MuMu16_phot_ind_LPair_Elastic_mar29.root";
    char *f2n= "../analyze/output_files/MuMu16_phot_ind_LPair_SemiInel_mar29.root";
    char *f3n= "../analyze/output_files/MuMu16_phot_ind_LPair_Inel_mar29.root";



    
    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    setup_all_SFs(year);





    TFile * f1 = TFile::Open(f1n, "READ");
    TTree *t1 = (TTree *) f1->Get("T_sig");

    TFile * f2 = TFile::Open(f2n, "READ");
    TTree *t2 = (TTree *) f2->Get("T_sig");

    TFile * f3 = TFile::Open(f3n, "READ");
    TTree *t3 = (TTree *) f3->Get("T_sig");

    int n_m_bins = 61;
    float m_bin_size = 30.;
	float m_bin_low = 170.;
	float m_bin_high = m_bin_low + n_m_bins*m_bin_size;
    TH1F *m1 = new TH1F("m1", "LPair Elastic", n_m_bins, m_bin_low, m_bin_high);
    TH1F *m2 = new TH1F("m2", "LPair Semi-Inelastic", n_m_bins, m_bin_low, m_bin_high);
    TH1F *m3 = new TH1F("m3", "LPair Inelastic", n_m_bins, m_bin_low, m_bin_high);

    int n_cost_bins = 16;
    float cost_bin_size = 2./n_cost_bins;
    TH1F *cost1 = new TH1F("cost1", "LPair Elastic", n_cost_bins, -1.,1.);
    TH1F *cost2 = new TH1F("cost2", "LPair Semi-Inelastic", n_cost_bins, -1.,1.);
    TH1F *cost3 = new TH1F("cost3", "LPair Inelastic", n_cost_bins, -1.,1.);


    int n_rap_bins = 20;
    float rap_bin_size = 5. / n_rap_bins;
    TH1F *rap1 = new TH1F("rap1", "LPair Elastic", n_rap_bins, -2.5,2.5);
    TH1F *rap2 = new TH1F("rap2", "LPair Semi-Inelastic", n_rap_bins, -2.5,2.5);
    TH1F *rap3 = new TH1F("rap3", "LPair Inelastic", n_rap_bins, -2.5,2.5);



    int n_pt_bins1 = 7;
    Float_t pt_bins1[] = {0., 10., 20., 30., 50., 70., 100., 300., 700. };
    TH1F *pt1 = new TH1F("pt1", "LPair Elastic", 10, 0, 200);
    TH1F *pt2 = new TH1F("pt2", "LPair Semi-Inelastic", 10, 0, 200);
    TH1F *pt3 = new TH1F("pt3", "LPair Inelastic", 10, 0, 200);


    int n_xf_bins1 = 5;
    float xf_bins1[] = {0.,0.04, 0.07, 0.1, 0.2, 0.5};
    TH1F *xf1 = new TH1F("xf1", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *xf2 = new TH1F("xf2", "dy signal", n_xf_bins1,  xf_bins1);

    int n_phi_bins = 20;
    TH1F *phi1 = new TH1F("phi1", "Data", n_phi_bins, -4.,4.);
    TH1F *phi2 = new TH1F("phi2", "Data", n_phi_bins, -4.,4.);

    float m_low = 170.;
    float m_high = 13000;
    bool ss = false;
    bool logy;
    float ratio_min = 0.5;
    float ratio_max = 1.5;


    char plt_file1[100];

    make_m_cost_pt_xf_hist(t1, m1, cost1, pt1, xf1, phi1, rap1, false, type, year, m_low, m_high);
    make_m_cost_pt_xf_hist(t2, m2, cost2, pt2, xf1, phi1, rap2, false, type, year, m_low, m_high);
    make_m_cost_pt_xf_hist(t3, m3, cost3, pt3, xf1, phi1, rap3, false, type, year, m_low, m_high);

    printf("Elastic:\n");
    cost1->Print("range\n");
    print_counting_AFB(cost1);
    printf("Semi-Elastic:\n");
    cost2->Print("range\n");
    print_counting_AFB(cost2);
    printf("Inelastic:\n");
    cost3->Print("range\n");
    print_counting_AFB(cost3);

    m1->SetLineColor(kBlue);
    m2->SetLineColor(kMagenta +3);
    m3->SetLineColor(kRed);

    cost1->SetLineColor(kBlue);
    cost2->SetLineColor(kMagenta +3);
    cost3->SetLineColor(kRed);


    rap1->SetLineColor(kBlue);
    rap2->SetLineColor(kMagenta +3);
    rap3->SetLineColor(kRed);

    pt1->SetLineColor(kBlue);
    pt2->SetLineColor(kMagenta +3);
    pt3->SetLineColor(kRed);


    m1->SetLineWidth(3);
    m2->SetLineWidth(3);
    m3->SetLineWidth(3);

    cost1->SetLineWidth(3);
    cost2->SetLineWidth(3);
    cost3->SetLineWidth(3);


    rap1->SetLineWidth(3);
    rap2->SetLineWidth(3);
    rap3->SetLineWidth(3);

    pt1->SetLineWidth(3);
    pt2->SetLineWidth(3);
    pt3->SetLineWidth(3);

    logy = true;
    TCanvas *c1 = new TCanvas("c1", "", 800,800);
    c1->SetLogy();
    m2->Draw("hist");
    m3->Draw("hist same");
    m1->Draw("hist same");
    m2->GetXaxis()->SetTitle("M (GeV)");
    sprintf(plt_file1, "%s%s_m_cmp.pdf", plot_dir, file_label);
    c1->BuildLegend();
    c1->Print(plt_file1);

    TCanvas *c2 = new TCanvas("c2", "", 800,800);
    cost2->Draw("hist");
    cost3->Draw("hist same");
    cost1->Draw("hist same");
    cost2->GetXaxis()->SetTitle("cos(#theta)");
    sprintf(plt_file1, "%s%s_cost_cmp.pdf", plot_dir, file_label);
    c2->BuildLegend();
    c2->Print(plt_file1);



    TCanvas *c3 = new TCanvas("c3", "", 800,800);
    c3->SetLogy();
    pt2->Draw("hist");
    pt3->Draw("hist same");
    pt1->Draw("hist same");
    pt2->GetXaxis()->SetTitle("dilepton p_T (GeV)");
    sprintf(plt_file1, "%s%s_pt_cmp.pdf", plot_dir, file_label);
    c3->BuildLegend();
    c3->Print(plt_file1);


    TCanvas *c4 = new TCanvas("c4", "", 800,800);
    c4->SetLogy();
    rap2->Draw("hist");
    rap3->Draw("hist same");
    rap1->Draw("hist same");
    rap2->GetXaxis()->SetTitle("dilepton Y");
    sprintf(plt_file1, "%s%s_rap_cmp.pdf", plot_dir, file_label);
    c4->BuildLegend();
    c4->Print(plt_file1);


}

    
    
