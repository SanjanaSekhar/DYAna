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
#include "../tdrstyle.C"
#include "../../utils/root_files.h"
#include "../../utils/HistUtils.C"
#include "../../utils/ScaleFactors.C"
#include "../../utils/PlotUtils.C"
#include "../../utils/TemplateMaker_systematics.C"



//make templates based on generator level samples (used for assessing impact of
//fiducial cuts on AFB
int fill_acceps(TTree *t_gen, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_rap,
        float m_low, float m_high, bool do_ptrw = false, int year = 2016){
    TH1F *h_m_all = (TH1F *) h_m->Clone("h_m_all");
    TH1F *h_cost_all = (TH1F *)h_cost->Clone("h_cost_all");
    TH1F *h_pt_all = (TH1F *) h_pt->Clone("h_pt_all");
    TH1F *h_rap_all = (TH1F *) h_rap->Clone("h_rap_all");

    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
    float mu_R_up, mu_R_down, mu_F_up, mu_F_down, mu_RF_up, mu_RF_down;
    float evt_weight;
    float pdf_weights[60];
    Bool_t sig_event(1);
    t_gen->SetBranchAddress("gen_p", &gen_lep_p);
    t_gen->SetBranchAddress("gen_m", &gen_lep_m);
    //t_gen->SetBranchAddress("gen_mu_p", &gen_lep_p);
    //t_gen->SetBranchAddress("gen_mu_m", &gen_lep_m);
    t_gen->SetBranchAddress("m", &m);
    t_gen->SetBranchAddress("cost", &cost);
    t_gen->SetBranchAddress("cost_st", &cost_st);
    t_gen->SetBranchAddress("gen_weight", &gen_weight);
    t_gen->SetBranchAddress("sig_event", &sig_event);

    A0_helpers A0_helper; 
    setup_A0_helper(&A0_helper, year);

    ptrw_helper ptrw_SFs; 
    setup_ptrw_helper(&ptrw_SFs, year);

    //float pt_cut = 26.;
    float pt_cut = 30.;


    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);


        if(m >= m_low && m <= m_high){
            nEvents++;
            cm = *gen_lep_p + *gen_lep_m;


            float pt = cm.Pt();
            float rap = cm.Rapidity();
            float evt_weight = gen_weight;

            if(do_ptrw){
                float ptrw = get_ptrw_SF(ptrw_SFs, m, pt, 0); 
                evt_weight *= ptrw;
            }

            h_m_all->Fill(m, evt_weight);
            h_cost_all->Fill(cost_st, evt_weight);
            h_pt_all->Fill(pt, evt_weight);
            h_rap_all->Fill(rap, evt_weight);

            float leading = std::fmax(gen_lep_m->Pt(), gen_lep_p->Pt());
            float subleading = std::fmin(gen_lep_m->Pt(), gen_lep_p->Pt());

            bool pass = true;
            pass = pass && abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4;
            pass = pass && leading > pt_cut && subleading > 15.;
            if(pass){
                h_m->Fill(m, evt_weight);
                h_cost->Fill(cost_st, evt_weight);
                h_pt->Fill(pt, evt_weight);
                h_rap->Fill(rap, evt_weight);
            }
        }

    }
    printf("selected %i events \n", nEvents);

    printf("All:\n");
    h_cost_all->Print("range");
    print_counting_AFB(h_cost_all);
    printf("\n\nFiducial:\n");
    h_cost->Print("range");
    print_counting_AFB(h_cost);


    h_m->Divide(h_m_all);
    h_cost->Divide(h_cost_all);
    h_pt->Divide(h_pt_all);
    h_rap->Divide(h_rap_all);

    return nEvents;

}


void draw_acceptance(){
    const int type = FLAG_MUONS;
    const int year = 2017;
    const bool write_out = true;
    char *plot_dir = "Paper_plots/acceptance/";
    char *file_label = "Accep17";
    char *plot_label = "M > 170 GeV";
    bool do_ptrw = true;


    TFile *f_gen = TFile::Open("../analyze/output_files/DY17_gen_level_nov13.root");

    gROOT->SetBatch(1);
    setTDRStyle();

    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    int n_m_bins = 61;
    float m_bin_size = 30.;
	float m_bin_low = 170.;
	float m_bin_high = m_bin_low + n_m_bins*m_bin_size;
    TH1F *m1 = new TH1F("m1", "", n_m_bins, m_bin_low, m_bin_high);

    int n_cost_bins = 40;
    TH1F *cost1 = new TH1F("cost1", "", n_cost_bins, -1.,1.);

    int n_rap_bins = 20;
    TH1F *rap1 = new TH1F("rap1", "", n_rap_bins, -3.0,3.0);




    int n_pt_bins1 = 7;
    Float_t pt_bins1[] = {0., 10., 20., 30., 50., 70., 100., 300., 700. };
    TH1F *pt1 = new TH1F("pt1", "", n_pt_bins1, pt_bins1);



    float m_low, m_high;

    m_low = 170.;
    m_high = 10000.;

    fill_acceps(t_gen_mu, m1, cost1, pt1, rap1, m_low, m_high, do_ptrw, year);

    cost1->Print("range");

    char plt_file1[100];
    bool logy = false;

    sprintf(plt_file1, "%s%s_m.png", plot_dir, file_label);
    draw_single_plot(plt_file1, m1, "M (GeV)", "Acceptance", plot_label, logy, write_out);

    sprintf(plt_file1, "%s%s_cost.png", plot_dir, file_label);
    draw_single_plot(plt_file1, cost1, "cos(#theta_{*})", "Acceptance", plot_label, logy, write_out);

    sprintf(plt_file1, "%s%s_pt.png", plot_dir, file_label);
    draw_single_plot(plt_file1, pt1, "p_{T} (GeV)", "Acceptance", plot_label, logy, write_out);

    sprintf(plt_file1, "%s%s_rap.png", plot_dir, file_label);
    draw_single_plot(plt_file1, rap1, "Rapidity", "Acceptance", plot_label, logy, write_out);


}

