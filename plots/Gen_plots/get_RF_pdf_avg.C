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
#include "../../utils/bins.h"

constexpr int n_pdfs = 60;

void print_profile(TProfile *h, int n_bins){
    printf("%s = { ", h->GetName());
    for(int i=1; i < n_bins +1; i++){
        printf("%.3f, ", h->GetBinContent(i));
    }
    printf("}; \n");
}

void fill_RF_pdf_hists(TTree *t_dy, TProfile *h_R_up, TProfile *h_R_down, TProfile *h_F_up, TProfile *h_F_down, 
                    TProfile *h_RF_up, TProfile *h_RF_down, TProfile *pdfs[60], TProfile *h_nnpdfrw = nullptr){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
    float mu_R_up, mu_R_down, mu_F_up, mu_F_down, mu_RF_down, mu_RF_up, nnpdf30_weight;
    Float_t pdf_weights[60];

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
    t_dy->SetBranchAddress("pdf_weights", &pdf_weights);
    if(h_nnpdfrw != nullptr){
        printf("setting up nnpdfweight\n");
        t_dy->SetBranchAddress("nnpdf30_weight", &nnpdf30_weight);
    }

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
            for(int i=0; i<n_pdfs; i++){
                float pdf_weight = pdf_weights[i];
                pdf_weight = std::max(std::min(pdf_weight, 3.0f), 0.333f);
                pdfs[i]->Fill(m, pdf_weight, gen_weight);
            }
            if(h_nnpdfrw != nullptr) h_nnpdfrw->Fill(m, nnpdf30_weight, gen_weight);

        }


        
    }
    printf("selected %i events \n", nEvents);

}


void get_RF_pdf_avg(){
    setTDRStyle();
    const int year = 2016;
    char *out_file = "../analyze/SFs/2016/RF_pdf_weights.root";
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_april11.root");
    const bool write_out = true;
    const bool print_nnpdf = true;

    TTree *t_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_el = (TTree *) f_gen->Get("T_gen_el");

    TProfile *h_pdfs[60];
    //double m_bins_d[] = {150, 170, 200,  250, 320, 510, 700, 1000, 14000};
    double m_bins_d[] = {100, 200,  400, 500, 700, 800, 1000, 1500, 2000, 3000}; //match mbins of MC samples
    double y_low = 0.5;
    double y_high = 1.5;

    TProfile *h_nnpdfrw = new TProfile("h_nnpdfrw", "h_nnpdfrw", n_m_bins, m_bins_d, y_low, y_high);
    TProfile *h_R_up = new TProfile("h_R_up", "h_R_up", n_m_bins, m_bins_d, y_low, y_high);
    TProfile *h_F_up = new TProfile("h_F_up", "h_F_up", n_m_bins, m_bins_d, y_low, y_high);
    TProfile *h_RF_up = new TProfile("h_RF_up", "h_RF_up", n_m_bins, m_bins_d, y_low, y_high);
    TProfile *h_R_down = new TProfile("h_R_down", "h_R_down", n_m_bins, m_bins_d, y_low, y_high);
    TProfile *h_F_down = new TProfile("h_F_down", "h_F_down", n_m_bins, m_bins_d, y_low, y_high);
    TProfile *h_RF_down = new TProfile("h_RF_down", "h_RF_down", n_m_bins, m_bins_d, y_low, y_high);

    char name[100];
    for(int i=0; i<n_pdfs; i++){
        sprintf(name, "h_pdf%i", i);
        h_pdfs[i] = new TProfile(name, "", n_m_bins, m_bins_d, y_low, y_high);
    }


    if(!print_nnpdf){
        fill_RF_pdf_hists(t_mu, h_R_up, h_R_down, h_F_up, h_F_down, h_RF_up, h_RF_down, h_pdfs);
        fill_RF_pdf_hists(t_el, h_R_up, h_R_down, h_F_up, h_F_down, h_RF_up, h_RF_down, h_pdfs);
    }
    else{
        fill_RF_pdf_hists(t_mu, h_R_up, h_R_down, h_F_up, h_F_down, h_RF_up, h_RF_down, h_pdfs, h_nnpdfrw);
        fill_RF_pdf_hists(t_el, h_R_up, h_R_down, h_F_up, h_F_down, h_RF_up, h_RF_down, h_pdfs, h_nnpdfrw);
    }


    h_R_up->Print("all");
    print_profile(h_R_up, n_m_bins);
    print_profile(h_R_down, n_m_bins);
    print_profile(h_F_up, n_m_bins);
    print_profile(h_F_down, n_m_bins);
    print_profile(h_RF_up, n_m_bins);
    print_profile(h_RF_down, n_m_bins);
    print_profile(h_pdfs[10], n_m_bins);

    if(print_nnpdf){
        h_nnpdfrw->Print("all");
        print_profile(h_nnpdfrw, n_m_bins);
    }

    if(write_out){
        TFile *f_out = TFile::Open(out_file, "RECREATE");
        h_R_up->Write();
        h_F_up->Write();
        h_RF_up->Write();
        h_R_down->Write();
        h_F_down->Write();
        h_RF_down->Write();
        for(int i=0; i<n_pdfs; i++){
            h_pdfs[i]->Write();
        }
    
        f_out->Close();
    }


}
