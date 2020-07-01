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

void print_profile(TProfile *h, int n_bins){
    printf("%s = { ", h->GetName());
    for(int i=1; i < n_bins +1; i++){
        printf("%.3f, ", h->GetBinContent(i));
    }
    printf("}; \n");
}

void fill_RF_hists(TTree *t_dy, TProfile *h_R_up, TProfile *h_R_down, TProfile *h_F_up, TProfile *h_F_down, 
                    TProfile *h_RF_up, TProfile *h_RF_down){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
    float mu_R_up, mu_R_down, mu_F_up, mu_F_down, mu_RF_down, mu_RF_up;

    t_dy->SetBranchAddress("gen_p", &gen_lep_p);
    t_dy->SetBranchAddress("gen_m", &gen_lep_m);
    t_dy->SetBranchAddress("m", &m);
    t_dy->SetBranchAddress("cost", &cost);
    t_dy->SetBranchAddress("cost_st", &cost_st);
    t_dy->SetBranchAddress("mu_R_up", &mu_R_up);
    t_dy->SetBranchAddress("mu_R_down", &mu_R_down);
    t_dy->SetBranchAddress("mu_F_up", &mu_F_up);
    t_dy->SetBranchAddress("mu_F_down", &mu_F_down);
    t_dy->SetBranchAddress("mu_RF_up", &mu_RF_up);
    t_dy->SetBranchAddress("mu_RF_down", &mu_RF_down);
    t_dy->SetBranchAddress("gen_weight", &gen_weight);

    float m_low = 150.;




    int nEvents=0;

    for (int i=0; i<t_dy->GetEntries(); i++) {
        t_dy->GetEntry(i);
        if(m > m_low){
            nEvents++;
            h_R_down->Fill(m, mu_R_down, gen_weight);
            h_F_down->Fill(m, mu_F_down, gen_weight);
            h_RF_down->Fill(m, mu_RF_down, gen_weight);
            h_R_up->Fill(m, mu_R_up, gen_weight);
            h_F_up->Fill(m, mu_F_up, gen_weight);
            h_RF_up->Fill(m, mu_RF_up, gen_weight);
        }


        
    }
    printf("selected %i events \n", nEvents);

}


void get_RF_avg(){
    setTDRStyle();
    TFile *f_gen = TFile::Open("../analyze/output_files/DY18_gen_level_june29.root");

    TTree *t_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_el = (TTree *) f_gen->Get("T_gen_el");


    int n_m_bins = 8;
    double m_bins[] = {150, 171, 200,  250, 320, 510, 700, 1000, 14000};

    TProfile *h_R_up = new TProfile("h_R_up", "h_R_up", n_m_bins, m_bins, 0.8, 1.2);
    TProfile *h_F_up = new TProfile("h_F_up", "h_F_up", n_m_bins, m_bins, 0.8, 1.2);
    TProfile *h_RF_up = new TProfile("h_RF_up", "h_RF_up", n_m_bins, m_bins, 0.8, 1.2);
    TProfile *h_R_down = new TProfile("h_R_down", "h_R_down", n_m_bins, m_bins, 0.8, 1.2);
    TProfile *h_F_down = new TProfile("h_F_down", "h_F_down", n_m_bins, m_bins, 0.8, 1.2);
    TProfile *h_RF_down = new TProfile("h_RF_down", "h_RF_down", n_m_bins, m_bins, 0.8, 1.2);

    fill_RF_hists(t_mu, h_R_up, h_R_down, h_F_up, h_F_down, h_RF_up, h_RF_down);
    fill_RF_hists(t_el, h_R_up, h_R_down, h_F_up, h_F_down, h_RF_up, h_RF_down);

    h_R_up->Print("all");
    print_profile(h_R_up, n_m_bins);
    print_profile(h_R_down, n_m_bins);
    print_profile(h_F_up, n_m_bins);
    print_profile(h_F_down, n_m_bins);
    print_profile(h_RF_up, n_m_bins);
    print_profile(h_RF_down, n_m_bins);


}
