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
#include "../../utils/root_files.h"
#include "../../utils/HistUtils.C"
#include "../../utils/PlotUtils.C"
#include "../../utils/ScaleFactors.C"

int make_data_hist(TTree *t_gen, int year, TH1F *h_cost,   
        float m_low, float m_high){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
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





    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        if(m >= m_low && m <= m_high && sig_event){
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            float my_cost = cost;
            if(gen_weight >0) nEvents++;
            else  nEvents--;

            h_cost->Fill(my_cost, gen_weight);
            //h_cost_r->Fill(cost, gen_weight);


            

        }
    }
    printf("selected %i events \n", nEvents);

    return nEvents;

}

int make_templates(TTree *t_gen, int year, TH1F *h_sym, TH1F *h_asym, TH1F *h_alpha,
        float m_low, float m_high){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
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



    A0_helpers A0_helper_;
    setup_A0_helper(&A0_helper_, year);


    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        if(m >= m_low && m <= m_high && sig_event){
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            float my_cost = cost;

            if(gen_weight >0) nEvents++;
            else  nEvents--;

            float denom = get_reweighting_denom(A0_helper_, cost_st, m, pt);
            float reweight_a = cost_st/ denom;
            float reweight_s = (1 + cost_st*cost_st)/denom;
            float reweight_alpha = (1 - cost_st*cost_st)/denom;

            h_sym->Fill(my_cost, reweight_s *gen_weight);
            h_sym->Fill(-my_cost, reweight_s *gen_weight);
            h_asym->Fill(my_cost, reweight_a *gen_weight);
            h_asym->Fill(-my_cost, -reweight_a *gen_weight);
            h_alpha->Fill(my_cost, reweight_alpha * gen_weight);
            h_alpha->Fill(-my_cost, reweight_alpha * gen_weight);



            

        }
    }

    return nEvents;

}



void templ_test(){

    int year = 2016;
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_april17.root");

    

    //TFile *f_gen = TFile::Open("../MuMu17_dy_gen.root");
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    char plot_dir[] = "Misc_plots/gen_templ_test/";



    int n_bins = 40;
    float bin_size = 2./n_bins;

    TH1F *h_cost_raw = new TH1F("h_cost1", "", n_bins, -1., 1.);
    TH1F *h_cost_templ = new TH1F("h_cost2", "", n_bins, -1., 1.);
    TH1F *h_cost_sym = new TH1F("h_sym", "", n_bins, -1., 1.);
    TH1F *h_cost_asym = new TH1F("h_asym", "", n_bins, -1., 1.);
    TH1F *h_cost_alpha = new TH1F("h_alpha", "", n_bins, -1., 1.);

    float afb = 0.61;
    float alpha = 0.06;


    float m_low = 150.;
    float m_high = 171.;

    printf("Mass range from %.0f to %.0f \n", m_low, m_high);

    make_data_hist(t_gen_mu,  year, h_cost_raw, m_low, m_high);
    make_data_hist(t_gen_el,  year, h_cost_raw, m_low, m_high);
    make_templates(t_gen_mu,  year, h_cost_sym, h_cost_asym, h_cost_alpha, m_low, m_high);
    make_templates(t_gen_el,  year, h_cost_sym, h_cost_asym, h_cost_alpha, m_low, m_high);

    h_cost_sym->Scale(0.5);
    h_cost_asym->Scale(0.5);
    h_cost_alpha->Scale(0.5);


    float norm = 3./4./(2.+alpha);
    //float norm = 3./8.;


    auto h_pl = *h_cost_sym + *h_cost_asym;
    auto h_mn = *h_cost_sym - *h_cost_asym;
    h_pl.Scale(0.5);
    h_mn.Scale(0.5);


    h_cost_templ->Add(&h_pl, &h_mn, (norm + afb), (norm - afb));
    //h_cost_asym->Scale(1./h_cost_sym->Integral());
    //h_cost_sym->Scale(1./h_cost_sym->Integral());
    //h_cost_templ->Add(h_cost_sym, h_cost_asym, norm, afb);

    h_cost_alpha->Scale(norm * alpha);
    h_cost_templ->Add(h_cost_alpha);

    /*
    TF1 *func = new TF1("func", "3./8.*(1 + x*x + ([1]/2.)*(1-3*x*x)) + [0]*x", -1., 1.);
    //TF1 *func = new TF1("func", "3./4./(2. + [1])*(1 + x*x + [1]*(1-x*x)) + [0]*x", -1., 1.);
    //TF1 *func = new TF1("func", "3./8.*(1 + x*x) + [0]*x", -1., 1.);
    //TF1 *func2 = new TF1("func", "[0]*x", -1., 1.);
    func->SetParameter(0,0.6);
    func->SetParameter(1,0.1);
    func->SetParName(0, "AFB");
    func->SetParName(1, "A0");
    h_cost_templ->Fit(func);
    h_cost_templ->Scale(1./h_cost_templ->Integral() / bin_size);

    gStyle->SetOptFit(1);
    TCanvas *c1 = new TCanvas("c1", "", 1000, 800);
    h_cost_templ->Draw();
    */
    //

    make_ratio_plot("Template Test", h_cost_raw, "Raw", h_cost_templ, "Sum of Templates", "ratio", "cos(#theta)", false);
}
