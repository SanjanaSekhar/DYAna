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

int make_amc_gen_cost(TTree *t_gen, TH1F *h_cost_st, TH1F *h_cost_r, TH1F *h_pt, TH1F *h_xf,  
        float m_low, float m_high){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    double gen_weight, m, cost, cost_st;
    t_gen->SetBranchAddress("gen_p", &gen_lep_p);
    t_gen->SetBranchAddress("gen_m", &gen_lep_m);
    t_gen->SetBranchAddress("m", &m);
    t_gen->SetBranchAddress("cost", &cost);
    t_gen->SetBranchAddress("cost_st", &cost_st);
    t_gen->SetBranchAddress("gen_weight", &gen_weight);





    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        if(m >= m_low && m <= m_high){
            if(gen_weight >0) nEvents++;
            else  nEvents--;
            cm = *gen_lep_p + *gen_lep_m;

            h_cost_st->Fill(cost_st, gen_weight);
            h_cost_r->Fill(cost, gen_weight);


            double xf = abs(2.*cm.Pz()/13000.);
            double pt = cm.Pt();

            h_pt->Fill(pt, gen_weight);
            h_xf->Fill(xf, gen_weight);

        }
    }
    printf("selected %i events \n", nEvents);

    return nEvents;

}




void fit_amc_gen_cost(){

    TFile *f_gen = TFile::Open("Results/DY_gen_level_nov8.root");
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    int n_bins = 80;
    TH1F *h_cost1 = new TH1F("h_cost1", "", n_bins, -1., 1.);
    TH1F *h_cost2 = new TH1F("h_cost2", "", n_bins, -1., 1.);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, -1., 1.);

    int nEvents = 0;


    for(int m_idx = 0; m_idx < n_m_bins; m_idx++){
        h_cost1->Reset();
        float m_low = m_bins[m_idx];
        float m_high = m_bins[m_idx+1];

        nEvents = make_amc_gen_cost(t_gen_mu,  h_cost1, h_dummy, h_dummy, h_dummy, m_low, m_high);
        nEvents += make_amc_gen_cost(t_gen_el,  h_cost1, h_dummy, h_dummy, h_dummy, m_low, m_high);
        //int nEvents2 = make_gen_cost(t_gen_el,  h_cost2, m_low, m_high);


        //TF1 *func = new TF1("func", "(1 + x*x + [1]*(1-x*x) + (4./3.)*(2. + [1])*[0]*x) /(8./3. + 4.*[1]/3.)", -1., 1.);
        TF1 *func = new TF1("func", "3./8.*(1 + x*x + ([1]/2.)*(1-3*x*x)) + [0]*x", -1., 1.);
        func->SetParameter(0,0.6);
        func->SetParameter(1,0.1);

        
        Double_t nB = h_cost1->Integral(1,n_bins/2);
        Double_t nF = h_cost1->Integral(n_bins/2 + 1,n_bins);
        
        float bin_size = 2./n_bins;
        //Normalize hist & divide by bin size
        h_cost1->Scale(1./h_cost1->Integral() / bin_size);
        h_cost1->Fit(func);
        //h_cost2->Scale(1./h_cost2->Integral() / bin_size);
        //h_cost2->Draw("same");
        //h_cost2->SetLineColor(kBlack);

        Double_t AFB = ((nF - nB))/((nF+nB));
        //Double_t dAFB = sqrt((1-AFB*AFB)/(nEvents));
        Double_t B = (1. - AFB)*nEvents/2.;
        Double_t F = (1. + AFB)*nEvents/2.;
        Double_t dAFB = sqrt(4.*F*B/pow(F+B, 3));
        //Double_t dAFB = sqrt(4.*nF*nB/pow(nF+nB, 3));

        printf("Mass range from %.0f to %.0f \n", m_low, m_high);
        printf("AFB: %.4f +/- %.4f \n", func->GetParameter(0), func->GetParError(0));
        printf("A0: %.3f +/- %.3f \n", func->GetParameter(1), func->GetParError(1));
        printf("Counting: NF %.0f NB %.0f \n", F, B);
        printf("Counting: AFB %.4f +/- %.4f \n", AFB, dAFB);
    }
}
