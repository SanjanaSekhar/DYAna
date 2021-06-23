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

const int type = FLAG_MUONS;

TFile *f_data, *f_mc, *f_mc_nosig, *f_ttbar, *f_QCD, *f_WJets, *f_WJets_mc, *f_QCD_mc, *f_diboson, *f_wt;
TTree *t_data, *t_mc, *t_mc_nosig, *t_ttbar, *t_QCD, *t_WJets, *t_WJets_mc, *t_QCD_mc, *t_diboson, *t_wt;

void make_hists(TTree *t1, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_rap, int flag1 = FLAG_MUONS, int year = 2018, bool require_filter = false){
    bool is_data = true;
    TempMaker tm(t1, is_data, year);
    if(flag1 == FLAG_MUONS) tm.do_muons = true;
    else tm.do_electrons = true;
    tm.do_RC = true;


    tm.setup();
    int pass_filters = 1;
    tm.t_in->SetBranchAddress("pass_filters", &pass_filters);
    int nEvents=0;

    int n1=0;
    int n2=0;
    int n3=0;


    Long64_t size  =  t1->GetEntries();
    for (int i=0; i<size; i++) {
        tm.getEvent(i);

        bool pass = tm.has_no_bjets && tm.met_pt< met_cut && tm.m >= 170.;
        if(require_filter) pass = pass && pass_filters;
        if(pass){
            nEvents++;
            tm.doCorrections();
            tm.getEvtWeight();

            h_m->Fill(tm.m,tm.evt_weight);
            h_pt->Fill(tm.pt,tm.evt_weight);
            h_cost->Fill(tm.cost, tm.evt_weight);
            h_rap->Fill(tm.cm.Rapidity(), tm.evt_weight);


        }
    }
    printf("Selected %i events \n", nEvents);
    //printf("n1,n2,n3: %i %i %i \n", n1,n2,n3);
    tm.finish();
}


void draw_met_filter_check(){
    int year = 2018;
    char *plot_dir = "Misc_plots/met_filter_check/";
    bool write_out = true;
    int flag1 = FLAG_ELECTRONS;

    //char *file_label = "MuMu18_filter_check";
    char *file_label = "ElEl18_filter_check";
    char *plt_label = "Electrons 2018";

    char *name1 = "No MET Filters";
    char *name2 = "MET Filters Required";

    gROOT->SetBatch(1);

    setTDRStyle();
    TTree *t;
    if(flag1 == FLAG_MUONS){
        TFile *f_mu = TFile::Open("../analyze/output_files/2018/MuMu18_filter_test.root");
        t = (TTree *) f_mu->Get("T_sig");
    }
    else{
        TFile *f_el = TFile::Open("../analyze/output_files/2018/ElEl18_filter_test.root");
        t = (TTree *) f_el->Get("T_sig");
    }



    


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
    TH1F *pt1 = new TH1F("pt1", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *pt2 = new TH1F("pt2", "dy signal", n_pt_bins1, pt_bins1);

    m1->Sumw2();
    cost1->Sumw2();
    pt1->Sumw2();
    rap1->Sumw2();

    m2->Sumw2();
    cost2->Sumw2();
    pt2->Sumw2();
    rap2->Sumw2();


    bool logy;
    float ratio_min = 0.99;
    float ratio_max = 1.01;


    char plt_file1[100];
    bool require_filters = false;


    make_hists(t, m1, cost1, pt1, rap1, flag1, year, require_filters);
    require_filters = true;
    make_hists(t, m2, cost2, pt2, rap2, flag1, year, require_filters);

    //zero errors on other hist to not double count errors
    bool add_err = false;
    setHistError(m1, 0, add_err);
    setHistError(cost1, 0, add_err);
    setHistError(pt1, 0, add_err);
    setHistError(rap1, 0, add_err);

    m1->Print("range");
    m2->Print("range");

    logy = true;
    sprintf(plt_file1, "%s%s_m_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, m1, name1, m2, name2, "Ratio", "M (GeV)", logy, write_out, ratio_min, ratio_max, plt_label);

    logy = false;
    sprintf(plt_file1, "%s%s_cost_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, cost1, name1, cost2, name2, "Ratio", "cos(#theta)", logy, write_out, ratio_min, ratio_max, plt_label);


    logy = true;
    sprintf(plt_file1, "%s%s_pt_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, pt1, name1, pt2, name2, "Ratio", "p_T (GeV)", logy, write_out, ratio_min, ratio_max, plt_label);

    logy = true;
    sprintf(plt_file1, "%s%s_rap_cmp.pdf", plot_dir, file_label);
    make_ratio_plot(plt_file1, rap1, name1, rap2, name2, "Ratio", "Rapidity", logy, write_out, ratio_min, ratio_max, plt_label);


}



