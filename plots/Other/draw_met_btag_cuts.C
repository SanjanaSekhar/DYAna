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

void make_cut_hists(TTree *t1, TH1F *h_both, TH1F* h_metonly, TH1F *h_btagonly , TH1F *h_neither, int flag1 = FLAG_MUONS, int year = 2016){
    bool is_data = false;
    TempMaker tm(t1, is_data, year);
    if(flag1 == FLAG_MUONS) tm.do_muons = true;
    else tm.do_electrons = true;
    tm.do_RC = true;


    tm.setup();
    float nnpdf30_weight = 1;
    int nEvents=0;

    int n1=0;
    int n2=0;
    int n3=0;


    Long64_t size  =  t1->GetEntries();
    for (int i=0; i<size; i++) {
        tm.getEvent(i);

        if(tm.m >= 170.){
            tm.doCorrections();
            tm.getEvtWeight();

            h_neither->Fill(tm.cost, tm.evt_weight);

            if(tm.has_no_bjets && tm.met_pt < met_cut ){
                n1++;
                h_both->Fill(tm.cost, tm.evt_weight);
            }
            if(tm.met_pt < met_cut){
                n2++;
                h_metonly->Fill(tm.cost, tm.evt_weight);
            }
            if(tm.has_no_bjets){
                n3++;
                h_btagonly->Fill(tm.cost, tm.evt_weight);
            }


        }
    }
    //printf("n1,n2,n3: %i %i %i \n", n1,n2,n3);
    tm.finish();
}


void draw_met_btag_cuts(){
    int year = 2018;
    char *plot_dir = "Paper_plots/";
    bool write_out = false;
    int flag1;


    setTDRStyle();
    init(year);


    int n_cos_theta_bins = 20;

    TH1F *h_mumu_mc_both = new TH1F("h_mumu_mc1", "Binned MC", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_mc_metonly = new TH1F("h_mumu_mc2", "Binned MC", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_mc_btagonly = new TH1F("h_mumu_mc3", "Binned MC", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_mc_neither = new TH1F("h_mumu_mc4", "Binned MC", n_cos_theta_bins, -1.,1.);

    TH1F *h_mumu_top_both = new TH1F("h_mumu_top1", "Binned top", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_top_metonly = new TH1F("h_mumu_top2", "Binned top", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_top_btagonly = new TH1F("h_mumu_top3", "Binned top", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_top_neither = new TH1F("h_mumu_top4", "Binned top", n_cos_theta_bins, -1.,1.);

    TH1F *h_mumu_diboson_both = new TH1F("h_mumu_diboson1", "Binned diboson", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_diboson_metonly = new TH1F("h_mumu_diboson2", "Binned diboson", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_diboson_btagonly = new TH1F("h_mumu_diboson3", "Binned diboson", n_cos_theta_bins, -1.,1.);
    TH1F *h_mumu_diboson_neither = new TH1F("h_mumu_diboson4", "Binned diboson", n_cos_theta_bins, -1.,1.);


    TH1F *h_elel_mc_both = new TH1F("h_elel_mc1", "Binned MC", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_mc_metonly = new TH1F("h_elel_mc2", "Binned MC", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_mc_btagonly = new TH1F("h_elel_mc3", "Binned MC", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_mc_neither = new TH1F("h_elel_mc4", "Binned MC", n_cos_theta_bins, -1.,1.);

    TH1F *h_elel_top_both = new TH1F("h_elel_top1", "Binned top", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_top_metonly = new TH1F("h_elel_top2", "Binned top", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_top_btagonly = new TH1F("h_elel_top3", "Binned top", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_top_neither = new TH1F("h_elel_top4", "Binned top", n_cos_theta_bins, -1.,1.);

    TH1F *h_elel_diboson_both = new TH1F("h_elel_diboson1", "Binned diboson", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_diboson_metonly = new TH1F("h_elel_diboson2", "Binned diboson", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_diboson_btagonly = new TH1F("h_elel_diboson3", "Binned diboson", n_cos_theta_bins, -1.,1.);
    TH1F *h_elel_diboson_neither = new TH1F("h_elel_diboson4", "Binned diboson", n_cos_theta_bins, -1.,1.);


    

    flag1 = FLAG_MUONS;

    make_cut_hists(t_mumu_mc, h_mumu_mc_both, h_mumu_mc_metonly, h_mumu_mc_btagonly, h_mumu_mc_neither, flag1, year);
    make_cut_hists(t_mumu_ttbar, h_mumu_top_both, h_mumu_top_metonly, h_mumu_top_btagonly, h_mumu_top_neither, flag1, year);
    make_cut_hists(t_mumu_wt, h_mumu_top_both, h_mumu_top_metonly, h_mumu_top_btagonly, h_mumu_top_neither, flag1, year);
    make_cut_hists(t_mumu_diboson, h_mumu_diboson_both, h_mumu_diboson_metonly, h_mumu_diboson_btagonly, h_mumu_diboson_neither, flag1, year);

    flag1 = FLAG_ELECTRONS;
    make_cut_hists(t_elel_mc, h_elel_mc_both, h_elel_mc_metonly, h_elel_mc_btagonly, h_elel_mc_neither, flag1, year);
    make_cut_hists(t_elel_ttbar, h_elel_top_both, h_elel_top_metonly, h_elel_top_btagonly, h_elel_top_neither, flag1, year);
    make_cut_hists(t_elel_wt, h_elel_top_both, h_elel_top_metonly, h_elel_top_btagonly, h_elel_top_neither, flag1, year);
    make_cut_hists(t_elel_diboson, h_elel_diboson_both, h_elel_diboson_metonly, h_elel_diboson_btagonly, h_elel_diboson_neither, flag1, year);

    

    printf("Year %i (DY mumu, top mumu, diboson mumu, DY ee, top ee, diboson ee): \n", year);
    printf("Neither & %.0f & %.0f & %.0f    & %.0f & %.0f & %.0f  \\\\ \n",    
            h_mumu_mc_neither->Integral(), h_mumu_top_neither->Integral(), h_mumu_diboson_neither->Integral(), 
            h_elel_mc_neither->Integral(), h_elel_top_neither->Integral(), h_elel_diboson_neither->Integral());
    printf("MET Cut & %.0f & %.0f & %.0f    & %.0f & %.0f & %.0f  \\\\ \n",    
            h_mumu_mc_metonly->Integral(), h_mumu_top_metonly->Integral(), h_mumu_diboson_metonly->Integral(),
            h_elel_mc_metonly->Integral(), h_elel_top_metonly->Integral(), h_elel_diboson_metonly->Integral());
    printf("Anti-Btag & %.0f & %.0f & %.0f  & %.0f & %.0f & %.0f  \\\\ \n",  
            h_mumu_mc_btagonly->Integral(), h_mumu_top_btagonly->Integral(), h_mumu_diboson_btagonly->Integral(),
            h_elel_mc_btagonly->Integral(), h_elel_top_btagonly->Integral(), h_elel_diboson_btagonly->Integral());
    printf("Both & %.0f & %.0f & %.0f       & %.0f & %.0f & %.0f  \\\\ \n",       
            h_mumu_mc_both->Integral(), h_mumu_top_both->Integral(), h_mumu_diboson_both->Integral(),
            h_elel_mc_both->Integral(), h_elel_top_both->Integral(), h_elel_diboson_both->Integral());




    char plt_file1[200], plt_file2[200];
    sprintf(plt_file1, "%sMuMu%i_met_cut.pdf", plot_dir, year % 2000);
    sprintf(plt_file2, "%sMuMu%i_btag_cut.pdf", plot_dir, year % 2000);

    float ratio_min = 0.9;
    float ratio_max = 1.1;
    
    bool logy = false;
    char plot_label[100];

    sprintf(plot_label, "Muons %i", year);
    make_ratio_plot(plt_file1, h_mumu_mc_btagonly, "Before",h_mumu_mc_both, "After MET Cut", "Before/After", "cos(#theta_{r})", logy, write_out,ratio_min,ratio_max, plot_label);
    make_ratio_plot(plt_file2, h_mumu_mc_metonly,  "Before",h_mumu_mc_both, "After anti-b-tag Cut", "Before/After", "cos(#theta_{r})", logy, write_out,ratio_min,ratio_max, plot_label);

    sprintf(plot_label, "Electrons %i", year);
    sprintf(plt_file1, "%sElEl%i_met_cut.pdf", plot_dir, year % 2000);
    sprintf(plt_file2, "%sElEl%i_btag_cut.pdf", plot_dir, year % 2000);
    make_ratio_plot(plt_file1, h_elel_mc_btagonly, "Before",h_elel_mc_both, "After MET Cut", "Before/After", "cos(#theta_{r})", logy, write_out,ratio_min,ratio_max, plot_label);
    make_ratio_plot(plt_file2, h_elel_mc_metonly,  "Before",h_elel_mc_both, "After anti-b-tag Cut", "Before/After", "cos(#theta_{r})", logy, write_out,ratio_min,ratio_max, plot_label);

}



