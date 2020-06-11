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

int make_amc_gen_cost(TTree *t_gen, int FLAG, TH2D *h_2d, float m_low, float m_high, int year){
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


    LQ_rw_helper h_LQ;
    setup_LQ_rw_helper(&h_LQ, year);



    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        //bool pass = abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4 && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > 30.;
        if(m >= m_low && m <= m_high && sig_event){
            nEvents++;
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            /*
            float my_cost = get_cost(*gen_lep_p, *gen_lep_m);
            if(cost_st > 0) my_cost = abs(my_cost);
            else my_cost = -abs(my_cost);
            */

            //gen_weight = gen_weight * 1000. / get_LQ_reweighting_denom(h_LQ, FLAG, m, cost_st);
            gen_weight = gen_weight * 1000. * get_LQ_reweighting_denom(h_LQ, FLAG_MUONS, m, cost_st) / get_LQ_reweighting_denom(h_LQ, FLAG_ELECTRONS, m, cost_st);
            float my_cost = cost_st;
            h_2d->Fill(m, my_cost, gen_weight);

        }
    }
    printf("selected %i events \n", nEvents);

    return nEvents;

}




void test_rw(){

    int year = 2018;
    TFile *f_gen = TFile::Open("../analyze/output_files/DY18_gen_level_april17.root");

    

    //TFile *f_gen = TFile::Open("../MuMu17_dy_gen.root");
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    char plot_dir[] = "Misc_plots/A0_fits/";



    int n_cost_bins = 20;
    float cost_bins[] = {-1.,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};
    //float cost_bins[] = {-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.0};
    //float cost_bins[] = {0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};

    int n_LQ_pt_bins = 1;
    float LQ_pt_bins[] = {0., 10000.};

    int n_LQ_m_bins = 37;
    float m_LQ_bins[] = {350., 375., 400., 425., 450., 475., 500., 525., 550., 575., 600., 625., 650., 675., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 
        1200., 1250., 1300., 1350., 1400., 1450., 1500., 1600., 1700., 1800., 1900.,  2000.,  2500.,   3000., };

    TH2D *h_mu = new TH2D("h_mu", "", n_LQ_m_bins, m_LQ_bins,  n_cost_bins, cost_bins);
    TH2D *h_el = new TH2D("h_el", "", n_LQ_m_bins, m_LQ_bins,  n_cost_bins, cost_bins);

    TH1D *h_1d = new TH1D("h1", "", n_cost_bins, cost_bins);


    //float m_low = m_LQ_bins[0];
    float m_low = 350;
    float m_high = 1225;
    make_amc_gen_cost(t_gen_mu, FLAG_MUONS,  h_mu, m_low, m_high, year);
    make_amc_gen_cost(t_gen_el,  FLAG_ELECTRONS, h_el, m_low, m_high, year);

    //h_el->Print("range");


    TCanvas *c_m = new TCanvas("c_m", "", 1000, 1000);
    TH1D * h_m = h_el->ProjectionX("h_m");
    h_m->GetXaxis()->SetRangeUser(500., 1000.);
    h_m->Print("range");
    
    h_m->Draw();

    TCanvas *c_cost = new TCanvas("c_cost", "", 1000, 1000);
    TH1D *h_cost = h_el->ProjectionY("h_cost");
    h_cost->Print("range");
    
    h_cost->Draw();

}
