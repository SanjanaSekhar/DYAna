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
#include "../../utils/ScaleFactors.C"
#include "../../utils/PlotUtils.C"
#include "../../utils/TemplateMaker_systematics.C"

int make_temps(TTree *t_gen, TH1F *h_raw, TH1F *h_sym, TH1F *h_asym, TH1F *h_alpha,  
        float m_low, float m_high, bool do_ptrw = false, int year = 2016){
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


    A0_helpers A0_helper; 
    setup_A0_helper(&A0_helper, year);

    ptrw_helper ptrw_SFs; 
    setup_ptrw_helper(&ptrw_SFs, year);


    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        bool pass = abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4 
            && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > 26. && min(gen_lep_m->Pt(), gen_lep_p->Pt()) > 15.;
        //bool pass = true;
        cm = *gen_lep_p + *gen_lep_m;
        //bool pass = abs(cm.Rapidity()) < 2.4;
        if(m >= m_low && m <= m_high && pass){

            float pt = cm.Pt();
            float rap = abs(cm.Rapidity());
            /*
            float my_cost = get_cost(*gen_lep_p, *gen_lep_m);
            if(cost_st > 0) my_cost = abs(my_cost);
            else my_cost = -abs(my_cost);
            */
            float gen_cost = cost_st;
            if(gen_weight >0) nEvents++;
            else  nEvents--;

            float denom = get_reweighting_denom(A0_helper, gen_cost, m, pt, rap);

            float reweight_a = gen_cost/ denom;
            float reweight_s = (1 + gen_cost*gen_cost)/denom;
            float reweight_alpha = (1 - gen_cost*gen_cost)/denom;


            if(do_ptrw){
                float ptrw = get_ptrw_SF(ptrw_SFs, m, pt, 0); 
                gen_weight *= ptrw;
            }


            h_raw->Fill(gen_cost, gen_weight);

            h_sym->Fill(gen_cost, reweight_s * gen_weight); 
            h_sym->Fill(-gen_cost, reweight_s * gen_weight); 

            h_asym->Fill(gen_cost, reweight_a * gen_weight);
            h_asym->Fill(-gen_cost, -reweight_a * gen_weight);

            h_alpha->Fill(gen_cost, reweight_alpha * gen_weight); 
            h_alpha->Fill(-gen_cost, reweight_alpha * gen_weight); 
             

        }
    }
    printf("selected %i events \n", nEvents);

    //cleanup_template(h_sym);
    //fixup_template_sum(h_sym, h_asym);

    return nEvents;

}




void make_gen_templates(){

    int year = 2016;
    bool do_ptrw = true;
    string fout_name = string("combine/templates/y16_gen_temps.root");
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_nov13.root");
    gROOT->SetBatch(1);

    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");

    TFile * fout = TFile::Open(fout_name.c_str(), "RECREATE");

    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;
    //int i_max = 1;

    float scale_ = 1e6;

    int n_bins = 40;
    TH1F *h_raw = new TH1F("cost_data_obs", "", n_bins, -1., 1.);
    TH1F *h_sym = new TH1F("cost_sym", "", n_bins, -1., 1.);
    TH1F *h_asym = new TH1F("cost_asym", "", n_bins, -1., 1.);
    TH1F *h_pl = new TH1F("cost_pl", "", n_bins, -1., 1.);
    TH1F *h_mn = new TH1F("cost_mn", "", n_bins, -1., 1.);
    TH1F *h_alpha = new TH1F("cost_alpha", "", n_bins, -1., 1.);

    float m_low, m_high;

    

    for(int i=i_start; i<i_max; i++){

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);



        make_temps(t_gen_mu, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year);
        make_temps(t_gen_el, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year);

        Double_t dnB, dnF;

        Double_t nB = h_raw->IntegralAndError(1,n_bins/2, dnB);
        Double_t nF = h_raw->IntegralAndError(n_bins/2 + 1,n_bins, dnF);
        Double_t n_tot = nB + nF;
        
        Double_t AFB = ((nF - nB))/((nF+nB));
        Double_t dAFB_v2 = sqrt( pow(dnB * 2. * nF / (n_tot*n_tot),2) + pow(dnF * 2. * nB / (n_tot*n_tot),2));

        printf("Counting AFB %.3f +/- %.3f \n", AFB, dAFB_v2);
        printf("F, B : %.3e %.3e \n", nF, nB);
        h_raw->Print("range");




        h_sym->Scale(0.5);
        h_asym->Scale(0.5);
        h_alpha->Scale(0.5);

        make_pl_mn_templates(h_sym, h_asym, h_pl, h_mn); 

        h_raw->Scale(scale_);
        h_pl->Scale(scale_);
        h_mn->Scale(scale_);
        h_alpha->Scale(scale_);

        fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

        h_raw->Write();
        h_pl->Write();
        h_mn->Write();
        h_alpha->Write();

        h_raw->Reset();
        h_pl->Reset();
        h_mn->Reset();
        h_alpha->Reset();
        h_asym->Reset();
        h_sym->Reset();


    }


    fout->Close();
    printf("Templates written to %s \n", fout_name.c_str());

}
